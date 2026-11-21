#=
What the dipole test (run_pipeline's max_negative) changes on campaign
sets with synthetic asteroids added:

    julia --project=. examples/dipole_test.jl

Declared before running:

Sets and injections as in real_bogus_dataset.jl (same seeds), without the
two sets analyze_set marks unreliable. Each set is searched three times by
analyze_set, no stacking: without the test, and with max_negative at 8
and at 6 (the detection threshold: a negative as deep as what counts as a
source). Counted each time: injected asteroids recovered
(injection_recovery), catalogued asteroids identified, and the other
tracklets, which are false but for the few uncatalogued asteroids.

Each set's counts are cached in data/real/dipole_test2, the first run's in
data/real/dipole_test.

Run again after inject_movers came to fit the frames together, as analyze_set
does (2026-11-19). On the sets kept aside (docs/src/validation/holdout.md) its
frame-by-frame fits were found up to 33" off, and on these sets 0.8-1.05"
from the joint fit on 5 frames and missing on 1 (XY72_p11): synthetic
asteroids were placed there up to 1" from where the search measures them.
Same sets, seeds and settings; nothing else changed.

The counts on the validation page are from a run made before pipeline_config
set max_negative, when the first search was analyze_set with its defaults.
The configuration is now built here, from the header's pixel scale, so that
the test can be switched off.
=#
using AsteroidPipeline, FITSIO, Printf, Random, Serialization

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const LEFT_OUT = ("XY14_p10", "XY73_p00")
const LIMITS = (Inf, 8.0, 6.0)
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const CACHE = joinpath(@__DIR__, "..", "data", "real", "dipole_test2")

orbits = read_mpcorb(download_mpcorb())
mkpath(CACHE)
total = Dict(l => zeros(Int, 3) for l in LIMITS)
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
        scale = FITS(injected[1]) do f
            AsteroidPipeline._pixel_scale_arcsec(load_wcs(read_header(f[1], String)), reverse(size(f[1])))
        end
        config = pipeline_config(scale; n_frames=length(injected))
        out = Dict{Float64,Vector{Int}}()
        for limit in LIMITS
            r = analyze_set(injected; workdir=mkpath(joinpath(work, "run$limit")), extract_to=work, observatory="F51",
                            mpcorb=orbits, name=name, stack=false, config=merge(config, (max_negative=limit,)))
            rec = injection_recovery(r.candidates, truth)
            ids = unique(r.candidates.id)
            mine = Set(rec.tracklet[rec.recovered])
            out[limit] = [count(rec.recovered), length(r.known), count(id -> !(id in mine) && !haskey(r.known, id), ids)]
        end
        serialize(cache, out)
        out
    end
    for l in LIMITS
        total[l] .+= counts[l]
    end
    @printf("%s: %s\n", name, join((@sprintf("%g: %d injected, %d known, %d other", l, counts[l]...) for l in LIMITS), " | "))
    flush(stdout)
end
@printf("all: %s\n", join((@sprintf("%g: %d injected, %d known, %d other", l, total[l]...) for l in LIMITS), " | "))
