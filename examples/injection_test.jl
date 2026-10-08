#=
Measures how complete the IASC configuration of `examples/iasc_demo.jl`
is, by injecting synthetic asteroids of known brightness and speed into
real IASC frames (`inject_movers`) and scoring which ones the pipeline
finds (`injection_recovery`) — the only way to get a completeness curve:
the six local sets contain six real moving objects in all.

Two sweeps per set: a main-belt-like one (5-100"/h) binned by magnitude
and speed, and a fast one (100-1000"/h, near-Earth-object rates) — fast
objects trail during a 45 s exposure, so this is what tells whether
trails survive detection, the sharpness and clustering filters, and
linking. `max_speed` is raised for the fast sweep: the IASC demo's limit
(~210"/h) would reject them outright.

Needs network (Gaia DR3 for each frame's WCS and zero point). Writes the
injected copies under data/real/iasc_injected/.

    julia --project=. examples/injection_test.jl [set ...]   # default: XY54_p10
=#
using AsteroidPipeline
using FITSIO, Statistics, Printf, Random

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "iasc")
const OUT_DIR = joinpath(@__DIR__, "..", "data", "real", "iasc_injected")
const PS1 = 0.2563                                   # arcsec/pixel
const ARCSEC_PER_HOUR = 24 / PS1                     # -> pixels/day, run_pipeline's unit

iasc_config(; max_speed_arcsec_h=210.0) = (
    threshold=6.0, match_radius=2.0 / PS1, max_speed=max_speed_arcsec_h * ARCSEC_PER_HOUR,
    mask_fill_values=true, refine_astrometry=true, difference_stack=true,
    min_sharpness=0.3, min_speed=3.0 * ARCSEC_PER_HOUR, max_flux_ratio=2.5, max_residual=1.0, min_frames=3)

function sweep(paths, label; n, mag_range, speed_range, max_speed_arcsec_h, seed, batches)
    results = []
    for b in 1:batches
        dir = joinpath(OUT_DIR, "$(label)_$b")
        injected, truth = inject_movers(paths, dir; n, mag_range, speed_range, rng=Xoshiro(seed + b))
        candidates = run_pipeline(injected; iasc_config(; max_speed_arcsec_h)...)
        push!(results, injection_recovery(candidates, truth))
        rm(dir; recursive=true)
    end
    rec = reduce(vcat, collect.(results))
    println("\n-- $label: $(count(r -> r.recovered, rec))/$(length(rec)) recovered --")
    return rec
end

function binned(rec, field, edges, fmt)
    for (lo, hi) in zip(edges[1:end-1], edges[2:end])
        sel = filter(r -> lo <= getproperty(r, field) < hi, rec)
        isempty(sel) && continue
        frac = count(r -> r.recovered, sel) / length(sel)
        @printf("  %s %-14s %3d/%-3d %5.0f%%  %s\n", field, Printf.format(Printf.Format(fmt), lo, hi), count(r -> r.recovered, sel),
                length(sel), 100frac, "#"^round(Int, 20frac))
    end
end

sets = isempty(ARGS) ? ["XY54_p10"] : ARGS
for set in sets
    paths = sort(filter(endswith(".fits"), readdir(joinpath(DATA_DIR, set), join=true)))
    println("\n=== $set ===")
    slow = sweep(paths, "$(set)_mb"; n=60, mag_range=(18.0, 23.0), speed_range=(5.0, 100.0),
                 max_speed_arcsec_h=210.0, seed=100, batches=3)
    binned(slow, :mag, 18.0:0.5:23.0, "%.1f-%.1f")
    binned(slow, :speed, [5, 10, 20, 40, 70, 100], "%.0f-%.0f\"/h")
    fast = sweep(paths, "$(set)_neo"; n=40, mag_range=(18.0, 21.0), speed_range=(100.0, 1000.0),
                 max_speed_arcsec_h=1100.0, seed=200, batches=2)
    binned(fast, :speed, [100, 200, 400, 700, 1000], "%.0f-%.0f\"/h")
end
