#=
What the star and faint tests (analyze_set's star_test and faint_test)
change on campaign sets with synthetic asteroids added:

    julia --project=. examples/star_test.jl

Found on other telescopes (docs/src/validation/telescopes.md, run with
examples/false_tracklets.jl on its three most crowded ZTF and Las Cumbres
sets): of 64 unidentified tracklets, all of 3 points on sets of 4-5 frames,
none found by stacking, most had every point within 4" of a Gaia star as
bright as the point, and most of the rest were fainter than the frames'
limit. star_test takes out the frame-by-frame points within their frame's
FWHM of a Gaia star no more than 0.5 mag fainter than the point, and a
tracklet left with fewer than 3; faint_test takes out the frame-by-frame
tracklets fainter than the limit (below_limit), as 3-frame sets already
did. Neither touches a tracklet found by stacking.

Declared before running:

Sets and injections as in dipole_test.jl (the sets of real_bogus_dataset.jl
and its seeds, V 20-23), without the two sets analyze_set marks unreliable.
Each set is searched three times by analyze_set, no stacking: with neither
test, with star_test, and with both (the default). Counted each time:
injected asteroids recovered (injection_recovery), catalogued asteroids
identified, and the other tracklets, which are false but for the few
uncatalogued asteroids.

Expected, all sets together, both tests against neither:
1. injected asteroids recovered: at least 99% as many;
2. catalogued asteroids identified: at most 1 fewer;
3. other tracklets: no more.

Each set's counts are cached in data/real/star_test.
=#
using AsteroidPipeline, FITSIO, Printf, Random, Serialization

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const LEFT_OUT = ("XY14_p10", "XY73_p00")
const TESTS = (neither=(star_test=false, faint_test=false), star=(star_test=true, faint_test=false),
               both=(star_test=true, faint_test=true))
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const CACHE = joinpath(@__DIR__, "..", "data", "real", "star_test")

orbits = read_mpcorb(download_mpcorb())
mkpath(CACHE)
total = Dict(t => zeros(Int, 3) for t in keys(TESTS))
for (k, name) in enumerate(SETS)
    name in LEFT_OUT && continue
    cache = joinpath(CACHE, name * ".jls")
    counts = if isfile(cache)
        deserialize(cache)
    else
        paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
        work = mktempdir()
        injected, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=100, mag_range=(20.0, 23.0),
                                        speed_range=(5.0, 100.0), rng=Xoshiro(2000 + k))
        out = Dict{Symbol,Vector{Int}}()
        for (t, kw) in pairs(TESTS)
            r = analyze_set(injected; workdir=mkpath(joinpath(work, "run_$t")), extract_to=work, observatory="F51",
                            mpcorb=orbits, name=name, stack=false, kw...)
            rec = injection_recovery(r.candidates, truth)
            ids = unique(r.candidates.id)
            mine = Set(rec.tracklet[rec.recovered])
            out[t] = [count(rec.recovered), length(r.known), count(id -> !(id in mine) && !haskey(r.known, id), ids)]
        end
        serialize(cache, out)
        out
    end
    for t in keys(TESTS)
        total[t] .+= counts[t]
    end
    @printf("%s: %s\n", name, join((@sprintf("%s: %d injected, %d known, %d other", t, counts[t]...) for t in keys(TESTS)), " | "))
    flush(stdout)
end
@printf("all: %s\n", join((@sprintf("%s: %d injected, %d known, %d other", t, total[t]...) for t in keys(TESTS)), " | "))
n, b = total[:neither], total[:both]
@printf("1. injected: %d of %d, %.1f%%; %s\n", b[1], n[1], 100b[1] / n[1], b[1] >= 0.99n[1] ? "met" : "not met")
@printf("2. catalogued: %d against %d; %s\n", b[2], n[2], b[2] >= n[2] - 1 ? "met" : "not met")
@printf("3. other: %d against %d; %s\n", b[3], n[3], b[3] <= n[3] ? "met" : "not met")
