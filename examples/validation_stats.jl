#=
Intervals for the validation pages:

    julia --project=. examples/validation_stats.jl

Included by other scripts for `wilson`, `v50` and `bootstrap_v50`; run on its
own, it prints the 95% intervals quoted on the pages, from the counts there
and, for the magnitudes where half are found, from the cached results of
faint_completeness.jl.

A fraction k of n gets Wilson's score interval (Wilson 1927, JASA 22, 209),
which stays inside 0–1 and is not empty at k = 0 or k = n. A magnitude where
half are found comes from asteroids injected into a few fields, which are
not independent of each other; its interval resamples whole sets with
replacement (bootstrap, 2,000 resamplings, percentiles 2.5 and 97.5).
=#
using Printf, Random, Serialization, Statistics

function wilson(k::Integer, n::Integer; z=1.959964)
    n == 0 && return (NaN, NaN)
    p = k / n
    centre = (p + z^2 / 2n) / (1 + z^2 / n)
    half = z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / 4n^2)
    (max(0.0, centre - half), min(1.0, centre + half))
end

const V50_EDGES = 21.0:0.25:23.5

function v50(mag, found; edges=V50_EDGES)
    f = map(edges[1:end-1]) do lo
        inbin = lo .<= mag .< lo + step(edges)
        n = count(inbin)
        (lo + step(edges) / 2, n == 0 ? NaN : count(found .& inbin) / n)
    end
    for i in 1:length(f)-1
        a, b = f[i], f[i+1]
        a[2] >= 0.5 > b[2] && return a[1] + (a[2] - 0.5) / (a[2] - b[2]) * (b[1] - a[1])
    end
    return NaN
end

"`sets` holds one (mag, found) per set; the sets are resampled, not the asteroids."
function bootstrap_v50(sets; n=2000, rng=Xoshiro(1), edges=V50_EDGES)
    draws = map(1:n) do _
        pick = rand(rng, sets, length(sets))
        v50(reduce(vcat, first.(pick)), reduce(vcat, last.(pick)); edges)
    end
    good = filter(!isnan, draws)
    (quantile(good, 0.025), quantile(good, 0.975), length(good))
end

pct(x) = @sprintf("%.0f", 100x)
interval(k, n) = ((lo, hi) = wilson(k, n); "$k of $n, $(pct(k / n))% ($(pct(lo))–$(pct(hi))%)")

const PUBLISHED = [
    ("pan-starrs1", "catalogued asteroids predicted on the frames, found", 17, 30),
    ("pan-starrs1", "synthetic asteroids V 18–21.5 found", 96, 117),
    ("stacking", "synthetic asteroids V 21–24 found frame by frame", 77, 320),
    ("stacking", "the same, with stacking", 94, 320),
    ("faint", "synthetic asteroids V 22.0–22.5 found frame by frame", 8, 214),
    ("faint", "the same, with stacking", 34, 214),
    ("variable-star", "confirmed variable stars found", 3, 3),
    ("orbits", "predictions 30 days ahead within 2σ", 39, 40),
    ("linking", "objects on 3 nights or more linked whole, none wrong", 47, 47),
    ("itf", "unidentified tracklets already in the ITF", 2, 5),
    ("digest2", "catalogued non-NEOs below 65", 25, 26),
    ("real-bogus", "synthetic asteroids found with the dipole test, over without", 511, 516),
    ("real-bogus", "other tracklets left by the dipole test", 19, 89),
    ("uncertainties", "points of numbered asteroids within 2 uncertainties", 41, 45),
    ("calibration", "catalogued asteroids of BANZAI's frames found here", 14, 18),
    ("ztf-periods", "periods within 1%, of those measured", 2, 6),
    ("period-sample", "period or its half within 1%", 13, 29),
    ("period-sample", "the same, magnitudes shuffled", 3, 29),
]

if abspath(PROGRAM_FILE) == @__FILE__
    for (page, what, k, n) in PUBLISHED
        @printf("%-14s %-62s %s\n", page, what, interval(k, n))
    end
    faint = joinpath(@__DIR__, "..", "data", "real", "faint_completeness")
    sets = filter(r -> !haskey(r, :failed) && r.reliable,
                  [deserialize(joinpath(faint, f)) for f in sort(readdir(faint)) if endswith(f, ".jls")])
    for (label, field) in (("frame by frame", :frame), ("with stacking", :stacked))
        per = [(r.mag, getfield(r, field)) for r in sets]
        lo, hi, kept = bootstrap_v50(per)
        @printf("faint: V50 %s %.2f (%.2f–%.2f, %d of 2000 resamplings cross 50%%), %d sets\n", label,
                v50(reduce(vcat, first.(per)), reduce(vcat, last.(per))), lo, hi, kept, length(per))
    end
end
