#=
What stack_search adds to the frame-by-frame search on one image set:

    julia --project=. examples/stacking_test.jl <set folder> [seed]

With a seed, 80 synthetic asteroids (V 21-24, 5-100"/h) are injected first,
and both searches are scored against them. Without one, the set is searched
as it is, and every tracklet found only by stacking is listed to look at.
The observatory is Pan-STARRS1 (F51), as for the IASC sets; change it for others.
=#
using AsteroidPipeline, Random, Statistics, Printf

length(ARGS) in (1, 2) || error("usage: julia --project=. examples/stacking_test.jl <set folder> [seed]")
folder = abspath(ARGS[1])
paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(folder; join=true)))
work = mktempdir()
truth = nothing
if length(ARGS) == 2
    paths, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=80, mag_range=(21.0, 24.0),
                                 speed_range=(5.0, 100.0), rng=Xoshiro(parse(Int, ARGS[2])))
end
orbits = read_mpcorb(download_mpcorb())
t0 = time()
r = analyze_set(paths; workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, name=basename(folder))
@printf("%s: %.0f s; %d tracklets, %d of them found only by stacking\n", r.name, time() - t0,
        length(unique(r.candidates.id)), length(r.stacked))

if truth === nothing
    for id in sort(collect(r.stacked))
        rows = r.candidates[r.candidates.id .== id]
        @printf("  stacked %d: RA %.5f Dec %.5f, V %.1f%s\n", id, rows.ra[1], rows.dec[1], (m = filter(isfinite, rows.mag); isempty(m) ? NaN : median(m)),
                haskey(r.known, id) ? ", " * r.known[id] : "")
    end
else
    all_found = injection_recovery(r.candidates, truth)
    frame_by_frame = injection_recovery(r.candidates[findall(id -> !(id in r.stacked), r.candidates.id)], truth)
    added = setdiff(Set(all_found.object[all_found.recovered]), Set(frame_by_frame.object[frame_by_frame.recovered]))
    real = Set(all_found.tracklet[all_found.recovered])
    # stacked positions against the injected ones, arcsec
    offsets = Float64[]
    for row in r.candidates
        row.id in r.stacked || continue
        d = minimum((3600 * hypot((t.ra - row.ra) * cosd(row.dec), t.dec - row.dec) for t in truth if t.frame == row.frame);
                    init=Inf)
        d < 3 && push!(offsets, d)
    end
    @printf("injected recovered: %d frame by frame, %d with stacking (V %s); stacked tracklets not injected: %d\n",
            count(frame_by_frame.recovered), count(all_found.recovered),
            join(sort([round(m; digits=1) for (o, m) in zip(all_found.object, all_found.mag) if o in added]), " "),
            count(id -> !(id in real), r.stacked))
    isempty(offsets) || @printf("stacked positions from the injected ones: %.2f\" median, %.2f\" max\n",
                                median(offsets), maximum(offsets))
end
