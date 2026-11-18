#=
Whether the positions and uncertainties analyze_set gives catalogued
asteroids agree with JPL's ephemerides, a truth independent of injection:

    julia --project=. examples/known_astrometry_test.jl

Declared before running:

Why. On injected asteroids (uncertainty_test.jl) the errors were 1.5 times
the uncertainties given and had a mean of +0.05" in RA and -0.05" in Dec.
The truth there is placed through each frame's own astrometric fit while
analyze_set measures through the joint one, so the offset could be the
test's or the pipeline's. This test has no injection.

Sets. The 13 sets in data/real/iasc, as they are. A set that analyze_set
marks unreliable is left out.

Search. analyze_set with its defaults (stacking on), observatory F51, the
MPC's orbits.

Measured. Every point of every tracklet analyze_set identifies as a
catalogued asteroid, against JPL Horizons' astrometric RA and Dec of that
asteroid from F51 at the point's epoch: the error in RA (times cos Dec) and
in Dec, in arcsec, divided by the uncertainty given; also along and across
the asteroid's motion (an epoch off by a constant shows along it only); the
mean error per set and over all. Horizons' own uncertainty is ignored: for a
numbered asteroid it is milliarcseconds; points of unnumbered ones are
reported apart. If the pipeline has no offset the mean errors are
consistent with 0 (within 2 standard errors over the sets' means).

Each set's identified points are cached in data/real/known_astrometry.
=#
using AsteroidPipeline, Printf, Serialization, Statistics
const P = AsteroidPipeline

const SETS = ["XY14_p10", "XY15_p01", "XY25_p10", "XY26_p01", "XY42_p11", "XY45_p01", "XY45_p10",
              "XY45_p11", "XY54_p10", "XY72_p01", "XY72_p11", "XY73_p00", "XY73_p01"]
const DATA = joinpath(@__DIR__, "..", "data", "real", "iasc")
const CACHE = joinpath(@__DIR__, "..", "data", "real", "known_astrometry")
const PROGRESS = joinpath(homedir(), ".claude", "progress", "a1-known")

function measure(name, orbits)
    paths = sort(filter(f -> endswith(lowercase(f), ".fits"), readdir(joinpath(DATA, name); join=true)))
    work = mktempdir()
    r = analyze_set(paths; workdir=work, extract_to=work, observatory="F51", mpcorb=orbits, name=name)
    r.reliable || return (name=name, reliable=false)
    rows = [(object=r.known[c.id], epoch=c.epoch, ra=c.ra, dec=c.dec, pos_err=c.pos_err)
            for c in r.candidates if haskey(r.known, c.id)]
    return (name=name, reliable=true, rows=rows)
end

function horizons_command(object)
    m = match(r"^\((\d+)\)", object)
    m === nothing || return "'$(m[1]);'"
    occursin(r"^\d{4} ", object) ? "'DES=$(strip(object));'" : "'$(strip(object));'"
end

# Horizons' astrometric RA, Dec (degrees) of `object` from F51 at each UTC Julian Date.
function jpl_radec(object, jds)
    q = Dict("COMMAND" => horizons_command(object), "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'F51@399'",
             "QUANTITIES" => "'1'", "ANG_FORMAT" => "DEG", "TIME_TYPE" => "UT", "CSV_FORMAT" => "YES",
             "OBJ_DATA" => "NO", "EXTRA_PREC" => "YES")
    out = Tuple{Float64,Float64}[]
    for result in P._horizons(q, [@sprintf("%.7f", jd) for jd in jds])
        for line in split(strip(P._horizons_table(result)), '\n')
            f = strip.(split(line, ','))
            push!(out, (parse(Float64, f[4]), parse(Float64, f[5])))
        end
    end
    length(out) == length(jds) || error("Horizons gave $(length(out)) positions for $(length(jds)) times")
    return out
end

mkpath(CACHE); mkpath(dirname(PROGRESS))
orbits = nothing
measured = map(enumerate(SETS)) do (k, name)
    write(PROGRESS, "$(k - 1) $(length(SETS)) conocidos\n")
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && return deserialize(cache)
    global orbits
    orbits === nothing && (orbits = read_mpcorb(download_mpcorb()))
    result = try
        measure(name, orbits)
    catch e
        (name=name, failed=first(sprint(showerror, e), 200))
    end
    serialize(cache, result)
    result
end
write(PROGRESS, "$(length(SETS)) $(length(SETS)) conocidos\n")

pts = NamedTuple[]
for r in measured
    haskey(r, :failed) && (println(r.name, ": left out, ", r.failed); continue)
    r.reliable || (println(r.name, ": left out, marked unreliable by analyze_set"); continue)
    for object in unique(p.object for p in r.rows)
        rows = sort(filter(p -> p.object == object, r.rows); by=p -> p.epoch)
        truth = try
            jpl_radec(object, [p.epoch for p in rows])
        catch e
            println(r.name, " ", object, ": no ephemeris, ", first(sprint(showerror, e), 120))
            continue
        end
        # direction of motion, from the ephemeris
        t1, t2 = first(truth), last(truth)
        mx = (t2[1] - t1[1]) * cosd(t1[2]); my = t2[2] - t1[2]
        u = hypot(mx, my) > 0 ? (mx, my) ./ hypot(mx, my) : (1.0, 0.0)
        for (p, (ra, dec)) in zip(rows, truth)
            dra = (p.ra - ra) * cosd(dec) * 3600; ddec = (p.dec - dec) * 3600
            push!(pts, (set=r.name, object=object, numbered=startswith(object, "("), dra=dra, ddec=ddec,
                        along=dra * u[1] + ddec * u[2], across=-dra * u[2] + ddec * u[1], pos_err=p.pos_err))
        end
    end
end

robust(z) = 1.4826 * median(abs.(z .- median(z)))
function summary(label, p)
    isempty(p) && return
    z = vcat([q.dra / q.pos_err for q in p], [q.ddec / q.pos_err for q in p])
    @printf("%-12s %4d points %3d objects | z std %5.2f robust %5.2f |z|<1 %3.0f%% <2 %3.0f%% | median |err| %.3f\" given %.3f\"\n",
            label, length(p), length(unique((q.set, q.object) for q in p)), std(z), robust(z),
            100mean(abs.(z) .< 1), 100mean(abs.(z) .< 2), median(abs.(vcat([q.dra for q in p], [q.ddec for q in p]))),
            median(q.pos_err for q in p))
    @printf("%-12s mean error RA %+.3f\" Dec %+.3f\" along %+.3f\" across %+.3f\"\n", "",
            mean(q.dra for q in p), mean(q.ddec for q in p), mean(q.along for q in p), mean(q.across for q in p))
end
println()
summary("all", pts)
summary("numbered", filter(q -> q.numbered, pts))
summary("unnumbered", filter(q -> !q.numbered, pts))
println("\nper set: mean RA, Dec error (arcsec)")
means = map(unique(q.set for q in pts)) do s
    p = filter(q -> q.set == s, pts)
    m = (mean(q.dra for q in p), mean(q.ddec for q in p))
    @printf("%-10s %4d points  %+.3f  %+.3f\n", s, length(p), m...)
    m
end
if length(means) >= 2
    se(v) = std(v) / sqrt(length(v))
    @printf("over the sets' means: RA %+.3f ± %.3f\", Dec %+.3f ± %.3f\" (standard error)\n",
            mean(first.(means)), se(first.(means)), mean(last.(means)), se(last.(means)))
end
