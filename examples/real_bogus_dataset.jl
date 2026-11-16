#=
Labelled tracklets for a real/false classifier, from every local
Pan-STARRS1 set (IASC):

    julia --project=. examples/real_bogus_dataset.jl

Declared before running:

Sets. The 13 sets in data/real/iasc, each with 100 synthetic asteroids
injected (inject_movers, V 20.0-23.0, 5-100"/h, seed 2000 plus the set's
place in the list: other seeds than the completeness test's). A set that
analyze_set marks unreliable is left out.

Search. run_pipeline on the frames analyze_set keeps, no stacking, with
the tracklet tests that decide what is real relaxed, so false tracklets
come through: no flux-ratio, background-excess or detector-artifact test,
the straight-line residual allowed up to 3 px. Detection as in production
(6 sigma, single-pixel spikes dropped).

Changed before any result was seen: relaxing detection too (5 sigma, no
sharpness test) left thousands of detections per frame and linking ran
for half an hour on one set without finishing. And the relaxed search was
analyze_set itself at first: it asks the catalogue services about every
unidentified tracklet, hundreds on a dense field, and stalled there.

Labels. A tracklet is real if injection_recovery matches it to an injected
asteroid, or it is one analyze_set, run with its defaults on the same
frames, identifies as a catalogued asteroid; false otherwise. (A catalogued
asteroid only the relaxed search finds is labelled false.)
An uncatalogued real asteroid would be labelled false: on these sets the
production search found 6 unidentified tracklets in 13 sets, so few.

Features, per tracklet, read from the frames' normalized differences
(stack_difference): significance at each point (minimum and median), the
ratio of the largest to the smallest aperture flux, the straight-line
residual, the rate, the fraction of masked pixels around the points, how
point-like the light is (peak over aperture sum, second-moment axis ratio),
the deepest negative within 2 FWHM (a misaligned star leaves a dipole),
and the number of points. Also whether the production search, run on the
same injected frames, keeps the tracklet.

Each set's table is cached in data/real/real_bogus.

The tables behind the validation page were made before pipeline_config set
max_negative: run again, the production search has that test too.
=#
using AsteroidPipeline, FITSIO, Printf, Random, Serialization, Statistics

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const CACHE = joinpath(@__DIR__, "..", "data", "real", "real_bogus")

relaxed(scale) = merge(pipeline_config(scale), (max_flux_ratio=Inf, max_residual=3.0, max_background_excess=Inf,
                                               detector_artifact_radius=0.0, max_negative=Inf))

read_image(p) = FITS(f -> permutedims(Float64.(read(f[1]))), p)

function features(rows, diffs, dmasks, scale)
    pts = sort(rows; by=r -> r.epoch)
    snr = Float64[]; flux = Float64[]; peakfrac = Float64[]; axis = Float64[]; masked = Float64[]; dip = Float64[]
    for p in pts
        d, m = diffs[p.frame], dmasks[p.frame]
        n1, n2 = size(d)
        x, y = round(Int, p.x), round(Int, p.y)
        r = 4
        (r < x <= n2 - r && r < y <= n1 - r) || continue
        box = d[y-r:y+r, x-r:x+r]
        mbox = m[y-r:y+r, x-r:x+r]
        v = ifelse.(mbox, 0.0, box)
        push!(masked, mean(mbox))
        push!(snr, sum(v[r-1:r+3, r-1:r+3]) / 5)
        pos = max.(v, 0.0)
        total = sum(pos)
        push!(flux, total)
        push!(peakfrac, total > 0 ? maximum(pos) / total : 0.0)
        push!(dip, minimum(v))
        if total > 0
            dx = (-r:r)'; dy = -r:r
            mx = sum(pos .* dx) / total; my = sum(pos .* dy) / total
            sxx = sum(pos .* (dx .- mx) .^ 2) / total; syy = sum(pos .* (dy .- my) .^ 2) / total
            sxy = sum(pos .* (dx .- mx) .* (dy .- my)) / total
            t = (sxx + syy) / 2; q = sqrt(max(((sxx - syy) / 2)^2 + sxy^2, 0.0))
            push!(axis, sqrt((t + q) / max(t - q, 1e-6)))
        else
            push!(axis, NaN)
        end
    end
    isempty(snr) && return nothing
    t = [(p.epoch - pts[1].epoch) * 24 for p in pts]
    xs = [(p.ra - pts[1].ra) * cosd(pts[1].dec) * 3600 for p in pts]
    ys = [(p.dec - pts[1].dec) * 3600 for p in pts]
    A = hcat(ones(length(t)), t)
    res = sqrt(mean(abs2, vcat(xs - A * (A \ xs), ys - A * (A \ ys))))
    rate = hypot((A \ xs)[2], (A \ ys)[2])
    pf = filter(>(0), flux)
    return (n=length(pts), snr_min=minimum(snr), snr_med=median(snr),
            flux_ratio=length(pf) >= 2 ? maximum(pf) / minimum(pf) : Inf, residual=res, rate=rate,
            masked=maximum(masked), peak_frac=median(peakfrac), axis_ratio=median_or_nan(filter(isfinite, axis)),
            dipole=minimum(dip))
end
median_or_nan(v) = isempty(v) ? NaN : median(v)

orbits = read_mpcorb(download_mpcorb())
mkpath(CACHE)
for (k, name) in enumerate(SETS)
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && continue
    paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
    work = mktempdir()
    result = try
        injected, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=100, mag_range=(20.0, 23.0),
                                        speed_range=(5.0, 100.0), rng=Xoshiro(2000 + k))
        prod = analyze_set(injected; workdir=mkpath(joinpath(work, "prod")), extract_to=work, observatory="F51",
                           mpcorb=orbits, name=name, stack=false)
        prod.reliable || error("marked unreliable by analyze_set")
        scale = AsteroidPipeline._pixel_scale_arcsec(prod.wcs[1], prod.image_size)
        used = prod.frames_used
        found = collect(run_pipeline(injected[used]; relaxed(scale)...))
        cands = [merge(c, (frame=used[c.frame],)) for c in found]
        rec = injection_recovery(cands, truth)
        injected_ids = Set(rec.tracklet[rec.recovered])
        prod_rec = injection_recovery(prod.candidates, truth)
        imgs = [read_image(p) for p in injected[used]]
        masks = [fill_value_mask(img) for img in imgs]
        d, dm = stack_difference(imgs, prod.wcs[used]; masks)
        diffs = Dict(zip(used, d)); dmasks = Dict(zip(used, dm))
        # a relaxed tracklet is one of the production search's if 3 of its points are within 1" of that one's
        function near(rows, tracklets)
            any(tracklets) do pts
                count(p -> any(q -> q.frame == p.frame &&
                                    AsteroidPipeline._angular_distance_arcsec(p.ra, p.dec, q.ra, q.dec) < 1.0, pts), rows) >= 3
            end
        end
        prod_all = collect(prod.candidates)
        prod_pts = [filter(c -> c.id == id, prod_all) for id in unique(prod.candidates.id)]
        known_pts = [filter(c -> c.id == id, prod_all) for id in keys(prod.known)]
        out = NamedTuple[]
        for id in unique(c.id for c in cands)
            rows = filter(c -> c.id == id, cands)
            f = features(rows, diffs, dmasks, scale)
            f === nothing && continue
            label = id in injected_ids ? :injected : near(rows, known_pts) ? :known : :false
            push!(out, merge((set=name, id=id, label=label, production=near(rows, prod_pts)), f))
        end
        @printf("%s: %d tracklets, %d injected, %d known, %d false; production found %d injected, %d tracklets\n",
                name, length(out), count(o -> o.label == :injected, out), count(o -> o.label == :known, out),
                count(o -> o.label == :false, out), count(prod_rec.recovered), length(unique(prod.candidates.id)))
        flush(stdout)
        (name=name, rows=out)
    catch e
        println(name, ": left out, ", first(sprint(showerror, e), 200))
        (name=name, failed=first(sprint(showerror, e), 200))
    end
    serialize(cache, result)
end
