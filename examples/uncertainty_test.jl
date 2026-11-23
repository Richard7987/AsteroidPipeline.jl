#=
Whether the uncertainties analyze_set gives each position and magnitude
(pos_err, mag_err: what ades_psv writes as rmsRA, rmsDec and rmsMag) are the
size of the errors actually made:

    julia --project=. examples/uncertainty_test.jl

Declared before running:

Sets. The 13 sets in data/real/iasc, each with 60 synthetic asteroids
injected (inject_movers, V 19.5-22.0, bright enough to be found often,
5-100"/h, seed 4000 plus the set's place in the list: other seeds than the
other tests'). A set that analyze_set marks unreliable is left out.

Search. analyze_set with its defaults (stacking on), observatory F51, the
MPC's orbits.

Measured. Every point of a tracklet that injection_recovery matches to an
injected asteroid, against that asteroid's true position and magnitude in
the same frame: the error in RA (times cos Dec) and in Dec, in arcsec, and in
magnitude, each divided by the uncertainty given. If the uncertainties are
right these ratios have a standard deviation of 1, 68% of them within 1 and
95% within 2. Reported for all points and by magnitude; the deviation also
robustly (1.4826 MAD), as a few mismatched points would inflate it.

The true positions are those inject_movers placed through each frame's own
astrometric fit, while analyze_set measures through the frames' joint fit:
the two fits' disagreement is part of the error measured, as it is part of
any real one. The magnitudes are both on the frame's Gaia zero point.

Each set's points are cached in data/real/uncertainty_test.

Changed after the first run (seeds 4000+k, 11 sets, 1899 points): positions
were given 0.018" where the errors were 0.073" (z robust deviation 5.5),
magnitudes nearly right (1.19). The astrometric part had been the fit's error
of its mean, rms over sqrt(n); an affine fit leaves the field's distortion in
its residuals, which a target's position carries as a star's does, so it is
now the stars' scatter about the fit, per axis. Checked on new seeds
(5000+k, cache data/real/uncertainty_test2), the first run kept as it was.

Run again after inject_movers came to fit the frames together, as analyze_set
does (2026-11-19). On the sets kept aside (docs/src/validation/holdout.md) its
frame-by-frame fits were found up to 33" off, and on these sets 0.8-1.05"
from the joint fit on 5 frames and missing on 1 (XY72_p11): synthetic
asteroids were placed there up to 1" from where the search measures them.
Same sets, seeds and settings; nothing else changed. Expected if that was
the cause: on synthetic asteroids the position errors over uncertainties
closer to a deviation of 1 than the 1.5 measured, and the mean offset of
0.05" smaller. Cache data/real/uncertainty_test3.

Run again after the candidates' photometry changed (2026-11-23): what lies
still under a candidate is subtracted, measured where it was in the other
frames, and the aperture's noise is measured in empty places of the frame
instead of computed from the sky's pixel scatter. Seen first on DECam
(examples/telescope_test.jl, 590 points): magnitude errors over
uncertainties with a deviation of 3.05 (robust 1.17), from points up to 1.3
mag too bright (a star in the aperture) and, without those, 1.35. Same sets,
seeds and settings as uncertainty_test3, so the same asteroids are measured
both ways. Expected: the magnitudes' deviation and robust deviation both
between 0.8 and 1.25, fewer than 1% of points beyond 5, and the positions
no worse than in uncertainty_test3 (deviation within 0.1 of it). Cache
data/real/uncertainty_test4.
=#
using AsteroidPipeline, Printf, Random, Serialization, Statistics

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const CACHE = joinpath(@__DIR__, "..", "data", "real", "uncertainty_test4")
const PROGRESS = joinpath(homedir(), ".claude", "progress", "a1-uncertainty")

orbits = read_mpcorb(download_mpcorb())
function measure(k, name)
    paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
    work = mktempdir()
    injected, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=60, mag_range=(19.5, 22.0),
                                    speed_range=(5.0, 100.0), rng=Xoshiro(5000 + k))
    r = analyze_set(injected; workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, name=name)
    r.reliable || return (name=name, reliable=false)
    rec = injection_recovery(r.candidates, truth)
    rows = collect(r.candidates)
    pts = NamedTuple[]
    for o in filter(o -> o.recovered, collect(rec))
        for t in filter(t -> t.object == o.object, collect(truth))
            i = findfirst(c -> c.id == o.tracklet && abs(c.epoch - t.epoch) < 1e-6, rows)
            i === nothing && continue
            c = rows[i]
            push!(pts, (mag_true=t.mag, dra=(c.ra - t.ra) * cosd(t.dec) * 3600, ddec=(c.dec - t.dec) * 3600,
                        dmag=c.mag - t.mag, pos_err=c.pos_err, mag_err=c.mag_err,
                        stacked=o.tracklet in r.stacked))
        end
    end
    @printf("%s: %d of 60 recovered, %d points\n", name, count(rec.recovered), length(pts))
    flush(stdout)
    (name=name, reliable=true, points=pts)
end

mkpath(CACHE); mkpath(dirname(PROGRESS))
measured = map(enumerate(SETS)) do (k, name)
    write(PROGRESS, "$(k - 1) $(length(SETS)) incertidumbres\n")
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && return deserialize(cache)
    result = try
        measure(k, name)
    catch e
        (name=name, failed=first(sprint(showerror, e), 200))
    end
    serialize(cache, result)
    result
end
write(PROGRESS, "$(length(SETS)) $(length(SETS)) incertidumbres\n")
for r in measured
    haskey(r, :failed) && println(r.name, ": left out, ", r.failed)
    haskey(r, :reliable) && !r.reliable && println(r.name, ": left out, marked unreliable by analyze_set")
end
pts = reduce(vcat, [r.points for r in measured if haskey(r, :points)])

robust(z) = 1.4826 * median(abs.(z .- median(z)))
function summary(label, sel)
    p = pts[sel]
    zra = [q.dra / q.pos_err for q in p if isfinite(q.pos_err)]
    zdec = [q.ddec / q.pos_err for q in p if isfinite(q.pos_err)]
    zmag = [q.dmag / q.mag_err for q in p if isfinite(q.mag_err) && isfinite(q.dmag)]
    z = vcat(zra, zdec)
    isempty(z) && return
    beyond5 = 100mean(abs.(zmag) .> 5)
    @printf("%-14s %5d | pos: std %5.2f robust %5.2f |z|<1 %3.0f%% <2 %3.0f%% median err %.3f\" given %.3f\" | mag: std %5.2f robust %5.2f |z|<1 %3.0f%% <2 %3.0f%% median err %.3f given %.3f, beyond 5 %.1f%%\n",
            label, length(p), std(z), robust(z), 100mean(abs.(z) .< 1), 100mean(abs.(z) .< 2),
            median(abs.(vcat([q.dra for q in p], [q.ddec for q in p]))), median(filter(isfinite, [q.pos_err for q in p])),
            std(zmag), robust(zmag), 100mean(abs.(zmag) .< 1), 100mean(abs.(zmag) .< 2),
            median(filter(isfinite, abs.([q.dmag for q in p]))), median(filter(isfinite, [q.mag_err for q in p])), beyond5)
end
println()
summary("all", trues(length(pts)))
for (lo, hi) in ((19.5, 20.5), (20.5, 21.5), (21.5, 22.0))
    summary(@sprintf("V %.1f-%.1f", lo, hi), [lo <= q.mag_true < hi for q in pts])
end
summary("frame by frame", [!q.stacked for q in pts])
summary("stacked", [q.stacked for q in pts])
@printf("RA and Dec separately: std %.2f and %.2f; mean error %.3f\" and %.3f\"\n",
        std([q.dra / q.pos_err for q in pts if isfinite(q.pos_err)]), std([q.ddec / q.pos_err for q in pts if isfinite(q.pos_err)]),
        mean(q.dra for q in pts), mean(q.ddec for q in pts))
