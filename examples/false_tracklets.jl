#=
What the tracklets analyze_set leaves unidentified are, on one set:

    julia --project=. examples/false_tracklets.jl <set folder> <MPC code>

A diagnostic, not a test: for each unidentified tracklet its points, rate,
magnitude against the limit, whether found only by stacking, and the
nearest Gaia star (G < 21) to its points, with that star's G.
=#
using AsteroidPipeline, FITSIO, Printf, Serialization, Statistics

const DIR = ARGS[1]
paths = sort(filter(f -> occursin(r"\.fits(\.fz)?$"i, f), readdir(DIR; join=true)))
saved = joinpath(tempdir(), "false_tracklets_" * basename(DIR) * ".jls")
r = if isfile(saved)
    deserialize(saved)
else
    work = mktempdir()
    serialize(saved, analyze_set(paths; workdir=work, extract_to=work, observatory=ARGS[2], name=basename(DIR)))
    deserialize(saved)
end
rows = collect(r.candidates)
ra0, dec0 = median(c.ra for c in rows), median(c.dec for c in rows)
radius = maximum(hypot((c.ra - ra0) * cosd(dec0), c.dec - dec0) for c in rows) + 0.02
function patient(f)
    for attempt in 1:12
        try
            return f()
        catch e
            attempt == 12 && rethrow()
            @warn "Gaia unavailable; trying again in 5 min" exception = e
            sleep(300)
        end
    end
end
stars = patient(() -> gaia_reference_stars(ra0, dec0, min(radius, 5.0); epoch_jd=median(c.epoch for c in rows),
                                            mag_range=(0.0, 21.0)))
@printf("%s: L %.2f, %d tracklets, %d unidentified, %d Gaia stars\n", basename(DIR), r.limiting_mag,
        length(unique(r.candidates.id)), length(r.unknown_ids), length(stars.ra))
println("  id  pts  rate\"/h  mag-L  stacked  star\"  starG  on a star  bright\"  brightG")
for id in r.unknown_ids
    t = sort(filter(c -> c.id == id, rows); by=c -> c.epoch)
    dt = (t[end].epoch - t[1].epoch) * 24
    sep(a, b) = hypot((a.ra - b.ra) * cosd(a.dec), a.dec - b.dec) * 3600
    rate = dt > 0 ? sep(t[end], t[1]) / dt : NaN
    d = [minimum(hypot.((stars.ra .- c.ra) .* cosd(c.dec), stars.dec .- c.dec)) * 3600 for c in t]
    j = argmin(hypot.((stars.ra .- t[1].ra) .* cosd(t[1].dec), stars.dec .- t[1].dec))
    mags = filter(isfinite, [c.mag for c in t])
    # points with a Gaia star within 4" as bright as the point or brighter
    on = count(t) do c
        near = hypot.((stars.ra .- c.ra) .* cosd(c.dec), stars.dec .- c.dec) .* 3600 .< 4
        any(near) && isfinite(c.mag) && minimum(stars.gmag[near]) <= c.mag + 0.5
    end
    # the brightest star within 60" of the first point
    near = findall(hypot.((stars.ra .- t[1].ra) .* cosd(t[1].dec), stars.dec .- t[1].dec) .* 3600 .< 60)
    b = isempty(near) ? 0 : near[argmin(stars.gmag[near])]
    bd = b == 0 ? NaN : hypot((stars.ra[b] - t[1].ra) * cosd(t[1].dec), stars.dec[b] - t[1].dec) * 3600
    @printf("%4d  %3d  %7.1f  %5.2f  %7s  %5.1f  %5.1f  %d of %d  %6.1f  %6.1f\n", id, length(t), rate,
            isempty(mags) ? NaN : median(mags) - r.limiting_mag, id in r.stacked, median(d), stars.gmag[j], on, length(t),
            bd, b == 0 ? NaN : stars.gmag[b])
end
