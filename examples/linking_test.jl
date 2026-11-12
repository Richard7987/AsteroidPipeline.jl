#=
Tests link_nights on real astrometry: one station's observations of every
catalogued asteroid in one patch of sky over two weeks, cut into one
tracklet per object and night, mixed, and linked without their names.

    julia --project=. examples/linking_test.jl

Declared before running:

Window. 2025-10-14 0h to 2025-10-29 0h UTC: the new moon of 2025-10-21 plus
and minus 7 days, when surveys observe most.

Objects. Every MPCORB orbit whose two-body position on 2025-10-21.5 UTC is
within 1 degree of the anti-solar point, where asteroids are densest.

Station. Pan-STARRS1 (F51) was declared first; it has no observation of
these objects in the window, so the station with the most CCD observations
of them there is used instead.

Tracklets. Their CCD observations (column 15 "C") from that station in the window,
one tracklet per object and night (nights split at local noon), kept with
2 points or more. Each tracklet gets a new number, so the linker sees no
names.

Linking. link_nights with its defaults then (sigma 0.5", 15 days, 3 sigma),
linkages of 2 nights or more.

Measured. Linkages whose tracklets are all one object (right) or not
(wrong); objects with tracklets on 2 nights or more whose tracklets were
all linked together, in part, or not at all.

Then again with false tracklets that move like asteroids added, 3 per real
one on each night (declared after the first run, before this one ran):
each has the epochs of one real tracklet of that night, the straight-line
motion of another, a start anywhere in the box that night's real points
cover, and 0.3" of noise per point; seeded, so the run repeats.

The first run fetches the observations from the MPC into data/real/linking;
JPL Horizons is asked live.
=#
using AsteroidPipeline, HTTP, JSON, Printf, Random, Statistics, TypedTables
const P = AsteroidPipeline

const START, STOP = 2460962.5, 2460977.5          # 2025-10-14, 2025-10-29 UTC
const MID = 2460970.0                              # 2025-10-21.5 UTC
const RADIUS = 1.0
const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "linking")

desig(o) = (m = match(r"^\((\d+)\)", o.name); m === nothing ? o.name : m[1])

function observations(d)
    path = joinpath(DATA_DIR, replace(d, ' ' => '_') * ".obs80")
    if !isfile(path)
        body = HTTP.request("GET", "https://data.minorplanetcenter.net/api/get-obs",
                            ["Content-Type" => "application/json"],
                            JSON.json(Dict("desigs" => [d], "output_format" => ["OBS80"]))).body
        write(path, JSON.parse(String(body))[1]["OBS80"])
    end
    rows = NamedTuple[]
    for line in eachline(path)
        length(line) == 80 && line[15] == 'C' || continue
        o = P._parse_mpc80_line(line)
        o === nothing || !(START <= o.jd < STOP) || push!(rows, o)
    end
    return rows
end

mkpath(DATA_DIR)
orbits = read_mpcorb(P.default_mpcorb_path())
E = P._observer_heliocentric_tt([MID + P._tt_minus_utc(MID)], nothing)[1]
ra0, dec0 = mod(atand(E[2], E[1]), 360), asind(E[3] / hypot(E...))
sep(ra, dec) = acosd(clamp(sind(dec0) * sind(dec) + cosd(dec0) * cosd(dec) * cosd(ra - ra0), -1, 1))
chosen = [o for o in orbits if (p = P._mpcorb_radec(o, MID, E); sep(p.ra, p.dec) <= RADIUS)]
@printf("anti-solar point RA %.3f Dec %.3f; %d objects within %.1f deg\n", ra0, dec0, length(chosen), RADIUS)

observed = Dict(o => observations(desig(o)) for o in chosen)
counts = Dict{String,Int}()
for obs in values(observed), r in obs
    counts[r.station] = get(counts, r.station, 0) + 1
end
const STATION = first(sort(collect(counts); by=p -> -p[2]))[1]
@printf("station %s, %d observations of them in the window\n", STATION, counts[STATION])
site = only(filter(c -> c.code == STATION, read_obscodes(joinpath(dirname(P.default_mpcorb_path()), "ObsCodes.txt"))))
night_of(jd) = floor(Int, jd + site.lon / 360)     # changes at local noon

tracklets = Dict{Tuple{String,Int},Vector{NamedTuple}}()
for (o, obs) in observed
    d = desig(o)
    for r in obs
        r.station == STATION || continue
        push!(get!(tracklets, (d, night_of(r.jd)), NamedTuple[]), r)
    end
end
filter!(p -> length(p[2]) >= 2, tracklets)

nights = sort(unique(n for (_, n) in keys(tracklets)))
truth = Dict{Tuple{Int,Int},String}()
movers = map(nights) do night
    rows = (id=Int[], ra=Float64[], dec=Float64[], epoch=Float64[])
    for (k, ((d, n), obs)) in enumerate(sort(collect(tracklets)))
        n == night || continue
        truth[(findfirst(==(night), nights), k)] = d
        for r in obs
            push!(rows.id, k); push!(rows.ra, r.ra); push!(rows.dec, r.dec); push!(rows.epoch, r.jd)
        end
    end
    Table(rows)
end
function run(movers, truth, label)
    per_object = Dict{String,Int}()
    for d in values(truth)
        startswith(d, "false") || (per_object[d] = get(per_object, d, 0) + 1)
    end
    linkable = Set(d for (d, c) in per_object if c >= 2)
    println("\n", label)
    @printf("%d tracklets (%d false) of %d objects on %d nights; %d objects on 2 nights or more (%d on 3 or more)\n",
            length(truth), count(startswith("false"), values(truth)), length(per_object), length(movers),
            length(linkable), count(d -> per_object[d] >= 3, linkable))

    elapsed = @elapsed linkages = link_nights(movers, STATION; min_nights=2)
    names = [unique(truth[(t.night, t.id)] for t in l.tracklets) for l in linkages]
    pure = [length(n) == 1 && !startswith(only(n), "false") for n in names]
    for (l, n, p) in zip(linkages, names, pure)
        p || println("wrong: ", join([truth[(t.night, t.id)] for t in l.tracklets], ", "), @sprintf("  rms %.2f\"", l.orbit.rms))
    end
    complete = Set(only(n) for (l, n, p) in zip(linkages, names, pure) if p && length(l.tracklets) == per_object[only(n)])
    partial = setdiff(Set(only(n) for (n, p) in zip(names, pure) if p), complete)
    missed = setdiff(linkable, complete, partial)
    used_false = count(startswith("false"), (truth[(t.night, t.id)] for l in linkages for t in l.tracklets))
    @printf("%d linkages in %.0f s: %d right, %d wrong; %d false tracklets linked\n",
            length(linkages), elapsed, count(pure), count(!, pure), used_false)
    long = [length(l.tracklets) >= 3 for l in linkages]
    @printf("  on 2 nights: %d right, %d wrong; on 3 or more: %d right, %d wrong\n", count(pure .& .!long),
            count(.!pure .& .!long), count(pure .& long), count(.!pure .& long))
    @printf("  %d points left out, in %d linkages\n", sum(l -> length(l.dropped), linkages; init=0), count(l -> !isempty(l.dropped), linkages))
    @printf("objects on 2+ nights: %d all tracklets linked, %d in part, %d not linked\n",
            length(complete), length(partial), length(missed))
    for k in 2:maximum(values(per_object); init=2)
        group = [d for d in linkable if per_object[d] == k]
        isempty(group) && continue
        @printf("  on %d nights: %d objects, %d complete, %d partial\n", k, length(group), count(in(complete), group), count(in(partial), group))
    end
    isempty(partial) || println("in part: ", join(sort(collect(partial)), ", "))
    isempty(missed) || println("not linked: ", join(sort(collect(missed)), ", "))
    rms = [l.orbit.rms for (l, p) in zip(linkages, pure) if p]
    isempty(rms) || @printf("orbit rms of the right linkages: median %.2f\", max %.2f\"\n", median(rms), maximum(rms))
end

run(movers, truth, "Real tracklets only")

# False tracklets moving like asteroids: each night, FALSE_PER_REAL per real
# tracklet, with the epochs of one real tracklet of that night, the motion
# (arcsec/day, a straight line) of another, a start anywhere in the box the
# night's real points cover, and 0.3" of noise on each point.
const FALSE_PER_REAL = 3
truth_false = Dict{Tuple{Int,Int},String}()
rng = Xoshiro(1)
false_movers = map(enumerate(movers)) do (night, rows)
    out = (id=collect(rows.id), ra=collect(rows.ra), dec=collect(rows.dec), epoch=collect(rows.epoch))
    ids = unique(rows.id)
    real = [filter(r -> r.id == id, rows) for id in ids]
    motion(k) = (m = tracklet_motion(k); 24 .* m.rate_arcsec_h .* (sind(m.pa_deg), cosd(m.pa_deg)))
    next = maximum(keys(truth))[2]
    for n in 1:FALSE_PER_REAL * length(ids)
        id = next + 1000 * night + n
        epochs = collect(rand(rng, real).epoch)
        rate = motion(rand(rng, real))
        ra0 = minimum(rows.ra) + rand(rng) * (maximum(rows.ra) - minimum(rows.ra))
        dec0 = minimum(rows.dec) + rand(rng) * (maximum(rows.dec) - minimum(rows.dec))
        for e in epochs
            dt = e - epochs[1]
            push!(out.id, id); push!(out.epoch, e)
            push!(out.dec, dec0 + (rate[2] * dt + 0.3 * randn(rng)) / 3600)
            push!(out.ra, ra0 + (rate[1] * dt + 0.3 * randn(rng)) / 3600 / cosd(dec0))
        end
        truth_false[(night, id)] = "false $night $n"
    end
    Table(out)
end
run(false_movers, merge(truth, truth_false), "With $FALSE_PER_REAL false tracklets per real one")
