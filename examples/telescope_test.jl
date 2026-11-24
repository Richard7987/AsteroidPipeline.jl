#=
Whether analyze_set works on telescopes other than Pan-STARRS1, measured the
same way on each:

    julia --project=. examples/telescope_test.jl <instrument> <MPC code> <seed> [set ...]

Named sets are measured alone, the others read from their cache; runs on
different sets can go side by side.

Declared before running, for every instrument; each instrument's sets are
chosen by the rules in examples/telescope_data.jl, fixed before its images
were looked at.

Sets. Every folder in data/real/telescopes/<instrument>/, as downloaded.

Search. analyze_set with its defaults (stacking on), the instrument's MPC
code (or, given "site", the code nearest each set's headers' site), the MPC's orbits. Nothing is set for an instrument: the settings scale
with the pixel scale the frames' WCS gives. A set analyze_set marks
unreliable, or that fails, is reported and left out of what follows.

Measured, on the frames as they are:
- the 5σ limiting magnitude L analyze_set gives; and, where the header has
  one, the survey's own (ZTF's MAGLIM);
- the catalogued asteroids predicted on every frame, 20 px or more from its
  edges, moving 5-100"/h, predicted V brighter than L - 0.5: found
  (`detected`) or not, in two bins of V - L, [-inf, -1.5) and [-1.5, -0.5);
- the tracklets analyze_set leaves unidentified;
- every point of every catalogued asteroid identified, against JPL
  Horizons' astrometric position from the instrument's site: the errors
  over pos_err.
Then, on the same frames with 100 synthetic asteroids injected
(inject_movers, V uniform in L - 2 to L + 0.5, 5-100"/h, seed <seed> plus
the set's place): found frame by frame and with stacking, as
injection_recovery counts them, by V - L; the V - L at which half are
found (validation_stats.jl's v50 on bins of 0.25 mag, sets resampled for
its interval); their positions' errors over pos_err.

Expected if the search works on this instrument as on Pan-STARRS1, where the
11 sets of the third run of faint.md had half found frame by frame a median
0.89 mag brighter than L (0.76-1.01), the 12 kept aside 0.94; 12
unidentified tracklets on 20 sets; on injected asteroids errors over
uncertainties with deviation 0.95 (98% within 2); against JPL 1.19:
1. half the injected asteroids found frame by frame between L - 1.2 and
   L - 0.6 (the median over the sets);
2. catalogued asteroids found in each bin as often as injected ones in the
   same bin, with stacking (their 95% intervals overlap);
3. at most 2 unidentified tracklets per set, on average;
4. injected asteroids' position errors over pos_err with deviation at most
   1.5, and 90% or more within 2;
5. catalogued asteroids' the same, against JPL, deviation at most 1.5.
Predicted V and the zero point's magnitudes differ by an asteroid's colour
(a few tenths), which the bins of 2 carry.

Each set's result is cached in data/real/telescope_test/<instrument>/.

Changed after the first DECam set, before any result of the test was seen:
inject_movers read a frame's first extension only, where DECam's and Las
Cumbres' files hold the image in the second (tile-compressed), and failed on
it; it now reads such a file as analyze_set does. That set's cache was
removed.

Changed after the first DECam set's result: 0 of 100 injected asteroids
recovered, the search's tracklets 0.07-0.46" from them. injection_recovery
matched points by frame number, inject_movers numbering the frames in the
order given and analyze_set in time order; DECam's files, named by their
checksum, do not sort by time. It now matches by epoch. The DECam cache was
removed; the frames of the other instruments sort by name in time order.

Changed during the second DECam set, with no result of it: a SkyBoT request
waited 3 h for an answer that never came, HTTP.jl's 60 s idle timeout not
ending the TLS read. SkyBoT and JPL requests now end at a deadline of their
own and are retried. Set 1's result, measured before, is kept.

Changed after the first summary of Las Cumbres: 3 of its 708 injected points
and 1 catalogued point had no uncertainty (an aperture that measured no
flux), which made expectations 4 and 5 NaN. Such points are now left out of
both, and counted.

Changed during the Las Cumbres sets: one (K2-25) failed when the ESA Gaia
archive answered 401, and the archive then answered with an HTML page behind
status 200. A copy that refuses, or answers with a page, is now passed over
for the next. ZTF ran with this change; K2-25 was run again with it, its
failed cache removed; the other sets' results are kept.
=#
using AsteroidPipeline, FITSIO, Printf, Random, Serialization, Statistics
const P = AsteroidPipeline
include("validation_stats.jl")

const INSTRUMENT = ARGS[1]
const CODE = ARGS[2]   # "site": each set's, the MPC code nearest its headers' site
const SEED = parse(Int, ARGS[3])
const ONLY = ARGS[4:end]   # sets to measure, if named; the others from their cache, or left out
const DATA = joinpath(@__DIR__, "..", "data", "real", "telescopes", INSTRUMENT)
const CACHE = joinpath(@__DIR__, "..", "data", "real", "telescope_test", INSTRUMENT)
const PROGRESS = joinpath(homedir(), ".claude", "progress", "b1-" * INSTRUMENT)
const REL_EDGES = -2.0:0.25:0.5

fits_in(dir) = sort(filter(f -> occursin(r"\.fits(\.fz)?$"i, f), readdir(dir; join=true)))

function header_value(path, key)
    FITS(path, "r") do f
        for hdu in f
            h = read_header(hdu)
            haskey(h, key) && return Float64(h[key])
        end
        NaN
    end
end

function image_size(path)
    FITS(path, "r") do f
        for hdu in f
            hdu isa ImageHDU && ndims(hdu) == 2 && return size(hdu)
        end
    end
end

function rate(o)
    dt = (o.epochs[end] - o.epochs[1]) * 24
    dt > 0 || return NaN
    hypot((o.ra[end] - o.ra[1]) * cosd(o.dec[1]), o.dec[end] - o.dec[1]) * 3600 / dt
end

function site_code(path)
    h = FITS(f -> read_header(f[1]), plain_frame(path; outdir=mktempdir()), "r")
    first(observatory_candidates(frame_site(h), read_obscodes(download_obscodes()))).code
end

function measure(k, name)
    paths = fits_in(joinpath(DATA, name))
    code = CODE == "site" ? site_code(first(paths)) : CODE
    work = mktempdir()
    t0 = time()
    r = analyze_set(paths; workdir=work, extract_to=work, observatory=code, mpcorb=orbits, name=name)
    seconds = time() - t0
    L = r.limiting_mag
    base = (name=name, code=code, frames=length(paths), reliable=r.reliable, limit=L, seconds=seconds,
            survey_limit=median(header_value.(paths, "MAGLIM")), warnings=r.warnings,
            unknown=length(r.unknown_ids))
    (r.reliable && isfinite(L)) || return base
    nx, ny = image_size(first(paths))
    inside(o) = all(20 .<= o.x .<= nx - 20) && all(20 .<= o.y .<= ny - 20)
    catalogued = [(name=o.name, rel=o.mag - L, rate=rate(o), found=o.detected)
                  for o in r.expected if inside(o) && 5 <= rate(o) <= 100 && o.mag < L - 0.5]
    known_rows = [(code=code, object=r.known[c.id], epoch=c.epoch, ra=c.ra, dec=c.dec, pos_err=c.pos_err)
                  for c in r.candidates if haskey(r.known, c.id)]

    injected, truth = inject_movers(paths, mkpath(joinpath(work, "injected")); n=100, mag_range=(L - 2, L + 0.5),
                                    speed_range=(5.0, 100.0), rng=Xoshiro(SEED + k))
    ri = analyze_set(injected; workdir=mkpath(joinpath(work, "run_injected")), extract_to=work, observatory=code,
                     mpcorb=orbits, name=name)
    rec = injection_recovery(ri.candidates, truth)
    frame = [o.recovered && !(o.tracklet in ri.stacked) for o in rec]
    rows = collect(ri.candidates)
    points = NamedTuple[]
    for o in filter(o -> o.recovered, collect(rec)), t in filter(t -> t.object == o.object, collect(truth))
        i = findfirst(c -> c.id == o.tracklet && abs(c.epoch - t.epoch) < 1e-6, rows)
        i === nothing && continue
        c = rows[i]
        push!(points, (dra=(c.ra - t.ra) * cosd(t.dec) * 3600, ddec=(c.dec - t.dec) * 3600, pos_err=c.pos_err,
                       dmag=c.mag - t.mag, mag_err=c.mag_err))
    end
    merge(base, (catalogued=catalogued, known_rows=known_rows, rel=collect(rec.mag) .- L, frame=frame,
                 stacked=collect(rec.recovered), points=points, injected_reliable=ri.reliable))
end

function jpl_radec(object, jds, code)
    m = match(r"^\((\d+)\)", object)
    command = m !== nothing ? "'$(m[1]);'" : occursin(r"^\d{4} ", object) ? "'DES=$(strip(object));'" : "'$(strip(object));'"
    q = Dict("COMMAND" => command, "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'$code@399'",
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
    out
end

robust(z) = 1.4826 * median(abs.(z .- median(z)))

sets = sort(filter(d -> isdir(joinpath(DATA, d)), readdir(DATA)))
mkpath(CACHE); mkpath(dirname(PROGRESS))
orbits = nothing
measured = map(enumerate(sets)) do (k, name)
    write(PROGRESS, "$(k - 1) $(length(sets)) $INSTRUMENT\n")
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && return deserialize(cache)
    isempty(ONLY) || name in ONLY || return nothing
    global orbits
    orbits === nothing && (orbits = read_mpcorb(download_mpcorb()))
    result = try
        measure(k, name)
    catch e
        (name=name, failed=first(sprint(showerror, e), 300))
    end
    serialize(cache, result)
    result
end
filter!(!isnothing, measured)
write(PROGRESS, "$(length(sets)) $(length(sets)) $INSTRUMENT\n")

println("\n$INSTRUMENT ($CODE), $(length(sets)) sets\n")
println("set                        frames  reliable  limit  survey's  unidentified  half found, V - L: frame  stacking   s")
good = []
for r in measured
    haskey(r, :failed) && (println(rpad(r.name, 27), "failed: ", r.failed); continue)
    ok = r.reliable && haskey(r, :rel)
    ok && push!(good, r)
    @printf("%-27s %6d  %8s  %5.2f  %8.2f  %12d  %24s  %8s  %4.0f\n", r.name, r.frames, r.reliable, r.limit,
            r.survey_limit, r.unknown, ok ? @sprintf("%.2f", v50(r.rel, r.frame; edges=REL_EDGES)) : "",
            ok ? @sprintf("%.2f", v50(r.rel, r.stacked; edges=REL_EDGES)) : "", r.seconds)
    foreach(w -> println("    ", first(w, 160)), r.warnings)
end
isempty(good) && exit()

rel = reduce(vcat, [r.rel for r in good]); frame = reduce(vcat, [r.frame for r in good])
stacked = reduce(vcat, [r.stacked for r in good])
println("\nV - L          injected  frame by frame  with stacking")
for lo in REL_EDGES[1:end-1]
    i = lo .<= rel .< lo + step(REL_EDGES)
    @printf("%+.2f..%+.2f  %8d  %6d (%3.0f%%)  %6d (%3.0f%%)\n", lo, lo + step(REL_EDGES), count(i),
            count(frame .& i), 100 * count(frame .& i) / max(count(i), 1), count(stacked .& i),
            100 * count(stacked .& i) / max(count(i), 1))
end

per_set = [v50(r.rel, r.frame; edges=REL_EDGES) for r in good]
lo, hi, nv = bootstrap_v50([(r.rel, r.frame) for r in good]; edges=REL_EDGES)
m = median(filter(isfinite, per_set))
@printf("\n1. half found frame by frame at V - L = %.2f, median over %d sets (all together %.2f, %.2f to %.2f, %d of 2000 resamplings); %s\n",
        m, count(isfinite, per_set), v50(rel, frame; edges=REL_EDGES), lo, hi, nv, -1.2 <= m <= -0.6 ? "met" : "not met")
lo2, hi2, nv2 = bootstrap_v50([(r.rel, r.stacked) for r in good]; edges=REL_EDGES)
@printf("   with stacking: %.2f (%.2f to %.2f)\n", v50(rel, stacked; edges=REL_EDGES), lo2, hi2)

cat = reduce(vcat, [r.catalogued for r in good])
println("\n2. catalogued against injected, with stacking:")
met2 = true
for (a, b) in ((-Inf, -1.5), (-1.5, -0.5))
    c = filter(o -> a <= o.rel < b, cat)
    i = max(a, -2.0) .<= rel .< b
    kc, kn = count(o -> o.found, c), length(c)
    ki, ni = count(stacked .& i), count(i)
    (cl, ch), (il, ih) = wilson(kc, kn), wilson(ki, ni)
    overlap = kn == 0 || (cl <= ih && il <= ch)
    global met2 &= overlap
    @printf("   V - L %s to %+.1f: catalogued %s; injected %s; %s\n", a == -Inf ? "-inf" : @sprintf("%+.1f", a), b,
            interval(kc, kn), interval(ki, ni), kn == 0 ? "none catalogued" : overlap ? "overlap" : "do not overlap")
    for o in filter(o -> !o.found, c)
        @printf("      not found: %s, V - L %+.2f, %.0f\"/h\n", o.name, o.rel, o.rate)
    end
end
println("   ", met2 ? "met" : "not met")

nu = sum(r.unknown for r in good)
@printf("\n3. unidentified tracklets: %d on %d sets, %.1f per set; %s\n", nu, length(good), nu / length(good),
        nu / length(good) <= 2 ? "met" : "not met")

all_pts = reduce(vcat, [r.points for r in good])
pts = filter(p -> isfinite(p.pos_err) && isfinite(p.mag_err), all_pts)
z = vcat([p.dra / p.pos_err for p in pts], [p.ddec / p.pos_err for p in pts])
zm = [p.dmag / p.mag_err for p in pts]
@printf("\n4. injected: %d points (%d without an uncertainty left out), position error over pos_err deviation %.2f (robust %.2f), %.0f%% within 2; %s\n",
        length(pts), length(all_pts) - length(pts), std(z), robust(z), 100 * mean(abs.(z) .<= 2), std(z) <= 1.5 && mean(abs.(z) .<= 2) >= 0.9 ? "met" : "not met")
@printf("   median error %.3f\", given %.3f\"; magnitudes: deviation %.2f (robust %.2f)\n",
        median(hypot.([p.dra for p in pts], [p.ddec for p in pts])) / sqrt(2), median([p.pos_err for p in pts]),
        std(zm), robust(zm))

jz = Float64[]; jerr = Float64[]
for r in good, object in unique(p.object for p in r.known_rows)
    rows = sort(filter(p -> p.object == object, r.known_rows); by=p -> p.epoch)
    truth = try
        jpl_radec(object, [p.epoch for p in rows], rows[1].code)
    catch e
        println("   ", r.name, " ", object, ": no ephemeris, ", first(sprint(showerror, e), 120))
        continue
    end
    for (p, (ra, dec)) in zip(rows, truth)
        isfinite(p.pos_err) || continue
        dra = (p.ra - ra) * cosd(dec) * 3600; ddec = (p.dec - dec) * 3600
        append!(jz, (dra / p.pos_err, ddec / p.pos_err)); append!(jerr, (dra, ddec))
    end
end
if isempty(jz)
    println("\n5. no catalogued asteroid identified")
else
    @printf("\n5. catalogued against JPL: %d points, error over pos_err deviation %.2f (robust %.2f), %.0f%% within 2; mean error %+.3f\"; %s\n",
            length(jz) ÷ 2, std(jz), robust(jz), 100 * mean(abs.(jz) .<= 2), mean(jerr), std(jz) <= 1.5 ? "met" : "not met")
end
