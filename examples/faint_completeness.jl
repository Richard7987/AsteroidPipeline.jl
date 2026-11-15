#=
How faint analyze_set finds asteroids, frame by frame and with stacking,
measured on every local Pan-STARRS1 set (IASC):

    julia --project=. examples/faint_completeness.jl

Declared before running:

Sets. The 13 sets in data/real/iasc: XY14_p10, XY15_p01, XY25_p10,
XY26_p01, XY42_p11, XY45_p01, XY45_p10, XY45_p11, XY54_p10, XY72_p01,
XY72_p11, XY73_p00, XY73_p01.

Injected. 100 synthetic asteroids per set with inject_movers, V 21.0 to
23.5 and 5 to 100"/h drawn uniformly, seed 1000 plus the set's place in
that list.

Search. analyze_set with its defaults (stacking on), observatory F51, the
MPC's orbits for identification.

Measured. An injected asteroid is recovered as injection_recovery counts
it. Recovered fractions in bins of 0.25 mag, by the frame-by-frame search
alone (tracklets not in `stacked`) and with stacking; the magnitude where
the fraction falls to 50%, interpolated between bins, for all sets
together and for each; the tracklets found only by stacking that are not
injected.

Changed after the first run stopped at XY72_p11 (inject_movers found no
photometric zero point in one of its frames), before any magnitude bin was
seen: a set whose injection fails, or that analyze_set marks unreliable
(XY14_p10 in that run), is reported and left out. Once inject_movers left such a
frame out instead of stopping, the cache of the two sets it had stopped on
(XY72_p11, XY73_p00) was removed and they were run again. Each set's result is
cached in data/real/faint_completeness.

The question it answers: how far past V 22, where Pan-STARRS's own
pipeline stops and IASC's searchers find their objects (Miller et al.
2024, PASP 136, 024502: 0.30 mag fainter on average), this pipeline finds
asteroids.
=#
using AsteroidPipeline, Printf, Random, Serialization, Statistics

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const EDGES = 21.0:0.25:23.5

const CACHE = joinpath(@__DIR__, "..", "data", "real", "faint_completeness")

orbits = read_mpcorb(download_mpcorb())
function measure(k, name)
    paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
    work = mktempdir()
    injected, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=100, mag_range=(21.0, 23.5),
                                    speed_range=(5.0, 100.0), rng=Xoshiro(1000 + k))
    t0 = time()
    r = analyze_set(injected; workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, name=name)
    all_found = injection_recovery(r.candidates, truth)
    frame = injection_recovery(r.candidates[findall(id -> !(id in r.stacked), r.candidates.id)], truth)
    real = Set(all_found.tracklet[all_found.recovered])
    extra = count(id -> !(id in real), r.stacked)
    @printf("%s: %.0f s, limit V %.2f; recovered %d frame by frame, %d with stacking; stacked not injected: %d\n",
            name, time() - t0, r.limiting_mag, count(frame.recovered), count(all_found.recovered), extra)
    (name=name, mag=all_found.mag, frame=frame.recovered, stacked=all_found.recovered, extra=extra,
     limit=r.limiting_mag, reliable=r.reliable)
end

mkpath(CACHE)
measured = map(enumerate(SETS)) do (k, name)
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
for r in measured
    haskey(r, :failed) && println(r.name, ": left out, ", r.failed)
    haskey(r, :reliable) && !r.reliable && println(r.name, ": left out, marked unreliable by analyze_set")
end
results = filter(r -> !haskey(r, :failed) && r.reliable, measured)

function fractions(mag, found)
    [(lo, count((lo .<= mag .< lo + 0.25)), count(found .& (lo .<= mag .< lo + 0.25))) for lo in EDGES[1:end-1]]
end

function v50(mag, found)
    f = [(lo + 0.125, n == 0 ? NaN : k / n) for (lo, n, k) in fractions(mag, found)]
    for i in 1:length(f)-1
        a, b = f[i], f[i+1]
        a[2] >= 0.5 > b[2] && return a[1] + (a[2] - 0.5) / (a[2] - b[2]) * (b[1] - a[1])
    end
    return NaN
end

mag = reduce(vcat, [r.mag for r in results])
frame = reduce(vcat, [r.frame for r in results])
stacked = reduce(vcat, [r.stacked for r in results])
println("\nV bin          injected  frame by frame  with stacking")
for ((lo, n, k1), (_, _, k2)) in zip(fractions(mag, frame), fractions(mag, stacked))
    @printf("%.2f-%.2f   %8d  %7d (%3.0f%%)  %7d (%3.0f%%)\n", lo, lo + 0.25, n, k1, 100k1 / n, k2, 100k2 / n)
end
@printf("50%% recovered at V %.2f frame by frame, V %.2f with stacking\n", v50(mag, frame), v50(mag, stacked))
println("\nset        limit  V50 frame  V50 stacking  stacked not injected")
for r in results
    @printf("%-9s  %5.2f  %9.2f  %12.2f  %d\n", r.name, r.limit, v50(r.mag, r.frame), v50(r.mag, r.stacked), r.extra)
end
