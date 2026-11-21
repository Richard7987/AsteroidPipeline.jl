#=
The dipole test and the depth of the search, measured again on sets kept
aside, which no setting was chosen on:

    julia --project=. examples/holdout_test.jl [run: 1 or 2]

Declared before running, with the list in docs/src/validation/holdout.md:

Sets. The 12 Pan-STARRS1 sets of IASC's October 2026 campaign that no
measurement of a setting used, in data/real/holdout (copied unchanged from
the files IASC sent): XY02_p00, XY04_p01, XY04_p10, XY04_p11, XY05_p00,
XY05_p01, XY55_p00, XY55_p01, XY55_p10, XY55_p11, XY75_p01, XY75_p10. They
were searched before only with analyze_set's defaults, for the ITF and
digest2 pages, and 3 of their candidates were looked at by eye. Nothing was
injected into them before.

Dipole test, as in dipole_test.jl. 100 synthetic asteroids per set, V 20.0
to 23.0 and 5 to 100"/h, seed 3000 plus the set's place in the list;
analyze_set without stacking, with max_negative off and at 6 (the default).
Counted: injected found (injection_recovery), catalogued identified, other
tracklets.

Depth, as in faint_completeness.jl. 100 per set, V 21.0 to 23.5 and 5 to
100"/h, seed 4000 plus the set's place; analyze_set with its defaults.
Counted: found by 0.25 mag, frame by frame and with stacking, and the
magnitude where half are found.

A set analyze_set marks unreliable, or whose injection fails, is reported
and left out. Nothing is changed in the pipeline whatever comes out.

Expected, from the 11 sets the settings were measured on (intervals 95%,
examples/validation_stats.jl):

1. Injected found with the test, over those found without it: 511 of 516
   there, 99.0%. Met if 99.0% is inside this run's Wilson interval.
2. Other tracklets left by the test: 19 of 89 there, 21%. Met if 21% is
   inside this run's interval; with fewer than 10 other tracklets without
   the test, not measurable.
3. Catalogued asteroids: the same identified with and without the test, as
   there (25 of 25).
4. Half found at V 21.67 frame by frame and V 21.79 with stacking there.
   Met if each is inside this run's bootstrap interval (sets resampled).

Each set's counts are cached in data/real/holdout_test.

Second run, declared after the first was seen (docs/src/validation/holdout.md).
The first found inject_movers fitting each frame's WCS on its own: on 6 of
these sets the headers are 4-33" off, its fits failed (no zero point on
XY04_p11 and XY05_p00) or landed about 1" from the frames' joint fit, and
synthetic asteroids were placed where the search does not measure them. It
now fits the frames together, as analyze_set does; the pipeline is
unchanged. The same measurements, seeds and criteria, with the expected
values taken from the 11 sets measured again with the same injector
(dipole_test.jl and faint_completeness.jl, their second runs, read from
their caches). Cache data/real/holdout_test2. Not a clean test of sets kept
aside: their first run was seen, and the injector was fixed on what it
showed.
=#
using AsteroidPipeline, FITSIO, Printf, Random, Serialization
include("validation_stats.jl")

const SETS = ["XY02_p00", "XY04_p01", "XY04_p10", "XY04_p11", "XY05_p00", "XY05_p01",
              "XY55_p00", "XY55_p01", "XY55_p10", "XY55_p11", "XY75_p01", "XY75_p10"]
const DATA = joinpath(@__DIR__, "..", "data", "real", "holdout")
const RUN = parse(Int, get(ARGS, 1, "2"))
const CACHE = joinpath(@__DIR__, "..", "data", "real", RUN == 1 ? "holdout_test" : "holdout_test2")
const TRAINING = joinpath(@__DIR__, "..", "data", "real")
const PROGRESS = joinpath(homedir(), ".claude", "progress", "a4-holdout")

orbits = read_mpcorb(download_mpcorb())

function dipole(k, name, paths, work)
    injected, truth = inject_movers(paths, mkpath(joinpath(work, "dipole")); n=100, mag_range=(20.0, 23.0),
                                    speed_range=(5.0, 100.0), rng=Xoshiro(3000 + k))
    scale = FITS(injected[1]) do f
        AsteroidPipeline._pixel_scale_arcsec(load_wcs(read_header(f[1], String)), reverse(size(f[1])))
    end
    config = pipeline_config(scale; n_frames=length(injected))
    map((Inf, 6.0)) do limit
        r = analyze_set(injected; workdir=mkpath(joinpath(work, "run$limit")), extract_to=work, observatory="F51",
                        mpcorb=orbits, name=name, stack=false, config=merge(config, (max_negative=limit,)))
        rec = injection_recovery(r.candidates, truth)
        mine = Set(rec.tracklet[rec.recovered])
        (reliable=r.reliable, found=rec.recovered, known=sort(collect(values(r.known))),
         other=count(id -> !(id in mine) && !haskey(r.known, id), unique(r.candidates.id)))
    end
end

function depth(k, name, paths, work)
    injected, truth = inject_movers(paths, mkpath(joinpath(work, "depth")); n=100, mag_range=(21.0, 23.5),
                                    speed_range=(5.0, 100.0), rng=Xoshiro(4000 + k))
    r = analyze_set(injected; workdir=mkpath(joinpath(work, "defaults")), extract_to=work, observatory="F51",
                    mpcorb=orbits, name=name)
    all_found = injection_recovery(r.candidates, truth)
    frame = injection_recovery(r.candidates[findall(id -> !(id in r.stacked), r.candidates.id)], truth)
    (reliable=r.reliable, mag=all_found.mag, frame=frame.recovered, stacked=all_found.recovered, limit=r.limiting_mag)
end

mkpath(CACHE)
measured = map(enumerate(SETS)) do (k, name)
    write(PROGRESS, "$(k - 1) $(length(SETS)) holdout$RUN\n")
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && return deserialize(cache)
    paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
    work = mktempdir()
    result = try
        (name=name, dipole=dipole(k, name, paths, work), depth=depth(k, name, paths, work))
    catch e
        (name=name, failed=first(sprint(showerror, e), 200))
    end
    serialize(cache, result)
    rm(work; recursive=true, force=true)
    result
end
write(PROGRESS, "$(length(SETS)) $(length(SETS)) holdout$RUN\n")

for r in measured
    haskey(r, :failed) && println(r.name, ": left out, ", r.failed)
    haskey(r, :failed) || all(d -> d.reliable, r.dipole) && r.depth.reliable ||
        println(r.name, ": left out, marked unreliable by analyze_set")
end
results = filter(r -> !haskey(r, :failed) && all(d -> d.reliable, r.dipole) && r.depth.reliable, measured)

println("\nset        off: injected  known  other | 6σ: injected  known  other | V50 frame  stacking  limit")
for r in results
    a, b = r.dipole
    @printf("%-9s  %13d  %5d  %5d | %12d  %5d  %5d | %9.2f  %8.2f  %5.2f\n", r.name, count(a.found), length(a.known),
            a.other, count(b.found), length(b.known), b.other, v50(r.depth.mag, r.depth.frame),
            v50(r.depth.mag, r.depth.stacked), r.depth.limit)
end

found_off = sum(count(r.dipole[1].found) for r in results)
found_on = sum(count(r.dipole[2].found) for r in results)
both = sum(count(r.dipole[1].found .& r.dipole[2].found) for r in results)
other_off = sum(r.dipole[1].other for r in results)
other_on = sum(r.dipole[2].other for r in results)
same_known = all(r -> r.dipole[1].known == r.dipole[2].known, results)
n_known = sum(length(r.dipole[1].known) for r in results)
inside(x, (lo, hi)) = lo <= x <= hi ? "met" : "not met"

cached(dir) = [deserialize(joinpath(TRAINING, dir, f)) for f in sort(readdir(joinpath(TRAINING, dir))) if endswith(f, ".jls")]
if RUN == 1
    base = (found_off=516, found_on=511, other_off=89, other_on=19, v50_frame=21.67, v50_stacked=21.79)
else
    dip = cached("dipole_test2")
    faint = filter(r -> !haskey(r, :failed) && r.reliable, cached("faint_completeness2"))
    pooled(field) = v50(reduce(vcat, [r.mag for r in faint]), reduce(vcat, [getfield(r, field) for r in faint]))
    base = (found_off=sum(c[Inf][1] for c in dip), found_on=sum(c[6.0][1] for c in dip),
            other_off=sum(c[Inf][3] for c in dip), other_on=sum(c[6.0][3] for c in dip),
            v50_frame=pooled(:frame), v50_stacked=pooled(:stacked))
    @printf("\nexpected, from %d and %d sets measured with the same injector: %d of %d; %d of %d; V %.2f and %.2f\n",
            length(dip), length(faint), base.found_on, base.found_off, base.other_on, base.other_off,
            base.v50_frame, base.v50_stacked)
end

println("\n", length(results), " sets")
w = wilson(found_on, found_off)
@printf("1. injected found with the test, over without: %s; %.1f%% %s\n", interval(found_on, found_off),
        100base.found_on / base.found_off, inside(base.found_on / base.found_off, w))
@printf("   of the %d found without it, %d also with it\n", found_off, both)
if other_off < 10
    @printf("2. other tracklets: %d without the test, %d with it; not measurable\n", other_off, other_on)
else
    @printf("2. other tracklets left by the test: %s; %.0f%% %s\n", interval(other_on, other_off),
            100base.other_on / base.other_off, inside(base.other_on / base.other_off, wilson(other_on, other_off)))
end
@printf("3. catalogued: %d without the test; the same with it: %s\n", n_known, same_known)

mag = reduce(vcat, [r.depth.mag for r in results])
println("\nV bin          injected  frame by frame  with stacking")
for lo in V50_EDGES[1:end-1]
    inbin = lo .<= mag .< lo + 0.25
    n = count(inbin)
    k1 = count(inbin .& reduce(vcat, [r.depth.frame for r in results]))
    k2 = count(inbin .& reduce(vcat, [r.depth.stacked for r in results]))
    @printf("%.2f-%.2f   %8d  %7d (%3.0f%%)  %7d (%3.0f%%)\n", lo, lo + 0.25, n, k1, 100k1 / n, k2, 100k2 / n)
end
for (label, field, before) in (("frame by frame", :frame, base.v50_frame), ("with stacking", :stacked, base.v50_stacked))
    per = [(r.depth.mag, getfield(r.depth, field)) for r in results]
    lo, hi, kept = bootstrap_v50(per)
    @printf("4. half found %s at V %.2f (%.2f–%.2f, %d of 2000 resamplings cross 50%%); V %.2f %s\n", label,
            v50(reduce(vcat, first.(per)), reduce(vcat, last.(per))), lo, hi, kept, before, inside(before, (lo, hi)))
end
