#=
Tests fit_orbit on real astrometry: the MPC's observations of known
asteroids, cut into short arcs as a new discovery would have, against JPL's
orbits for the same objects.

    julia --project=. examples/orbit_test.jl

Declared before running any of them:

Objects. The numbered asteroids 1000, 2000, ..., 10000 (chosen by number
only, mostly main belt), the near-Earth asteroids (433) Eros, (1036) Ganymed
and (1627) Ivar, and the Jupiter trojans (588) Achilles and (624) Hektor.

Arcs. Their CCD observations (column 15 "C") from 2024-01-01 on, from the
first night: the first 3 nights, every observation within 10 days, and
every one within 30 days, each fitted with sigma = 0.5" per observation.
An arc with fewer than 3 nights is reported as such.

Measured. The fitted osculating a, e, i against JPL Horizons' at the same
epoch, and the position predicted 30 days after the arc's end, seen from the
geocentre, against Horizons' ephemeris there, in arcsec and in units of the
fit's own 1-sigma uncertainty. Every arc is reported.

The first run fetches the observations from the MPC (no authentication) into
data/real/orbits; JPL Horizons is asked live.
=#
using AsteroidPipeline, HTTP, JSON, Printf, Statistics
const P = AsteroidPipeline

const OBJECTS = vcat(string.(1000:1000:10000), ["433", "1036", "1627", "588", "624"])
const START = 2460310.5                      # 2024-01-01 UTC
const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "orbits")

function observations(number)
    path = joinpath(DATA_DIR, "$number.obs80")
    if !isfile(path)
        body = HTTP.request("GET", "https://data.minorplanetcenter.net/api/get-obs",
                            ["Content-Type" => "application/json"],
                            JSON.json(Dict("desigs" => [number], "output_format" => ["OBS80"]))).body
        write(path, JSON.parse(String(body))[1]["OBS80"])
    end
    rows = NamedTuple[]
    for line in eachline(path)
        length(line) == 80 && line[15] == 'C' || continue
        o = P._parse_mpc80_line(line)
        o === nothing || o.jd < START || push!(rows, o)
    end
    return sort(rows; by=o -> o.jd)
end

# Horizons' osculating heliocentric ecliptic a, e, i of `number` at `jd` (TDB, within 2 ms of TT).
function jpl_elements(number, jd)
    q = Dict("COMMAND" => "'$number;'", "EPHEM_TYPE" => "ELEMENTS", "CENTER" => "'500@10'", "REF_PLANE" => "ECLIPTIC",
             "REF_SYSTEM" => "ICRF", "TIME_TYPE" => "TDB", "OUT_UNITS" => "AU-D", "CSV_FORMAT" => "YES", "OBJ_DATA" => "NO")
    result = only(P._horizons(q, [@sprintf("%.6f", jd)]))
    f = strip.(split(strip(match(r"\$\$SOE(.*)\$\$EOE"s, result)[1]), ','))
    return (e=parse(Float64, f[3]), i=parse(Float64, f[5]), a=parse(Float64, f[12]))
end

# Horizons' astrometric RA, Dec (degrees) of `number` from the geocentre at UTC `jd`.
function jpl_radec(number, jd)
    q = Dict("COMMAND" => "'$number;'", "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'500@399'", "QUANTITIES" => "'1'",
             "ANG_FORMAT" => "DEG", "TIME_TYPE" => "UT", "CSV_FORMAT" => "YES", "OBJ_DATA" => "NO")
    result = only(P._horizons(q, [@sprintf("%.6f", jd)]))
    f = strip.(split(strip(match(r"\$\$SOE(.*)\$\$EOE"s, result)[1]), ','))
    return (parse(Float64, f[4]), parse(Float64, f[5]))
end

nights(obs) = cumsum([true; diff([o.jd for o in obs]) .> 0.5])

mkpath(DATA_DIR)
for number in OBJECTS
    obs = try
        observations(number)
    catch e
        @printf("(%s): observations not fetched: %s\n", number, sprint(showerror, e)[1:min(end, 120)])
        continue
    end
    isempty(obs) && (@printf("(%s): no CCD observations since 2024\n", number); continue)
    night = nights(obs)
    for (label, keep) in (("3 nights", night .<= 3), ("10 days", [o.jd - obs[1].jd <= 10 for o in obs]),
                          ("30 days", [o.jd - obs[1].jd <= 30 for o in obs]))
        arc = obs[keep]
        nn = length(unique(night[keep]))
        if nn < 3
            @printf("(%s) %-8s %d observations over %d nights: too few nights\n", number, label, length(arc), nn)
            continue
        end
        try
            f = fit_orbit([o.ra for o in arc], [o.dec for o in arc], [o.jd for o in arc], [o.station for o in arc])
            jpl = jpl_elements(number, f.epoch)
            later = arc[end].jd + 30
            pred = orbit_ephemeris(f, [later], "500")
            ra, dec = jpl_radec(number, later)
            off = 3600 * hypot((pred.ra[1] - ra) * cosd(dec), pred.dec[1] - dec)
            el, s = f.elements, f.sigma_elements
            @printf("(%s) %-8s %3d obs, %2d nights, %5.1f d, rms %.2f\", %d rejected; a %.4f±%.4f (JPL %.4f), e %.4f±%.4f (%.4f), i %.3f±%.3f (%.3f); +30 d: %.1f\" off, σ %.1f\" (%.1fσ)\n",
                    number, label, length(arc), nn, f.arc_days, f.rms, count(!, f.used), el.a, s.a, jpl.a, el.e, s.e, jpl.e,
                    el.i, s.i, jpl.i, off, pred.sigma[1], off / pred.sigma[1])
        catch e
            @printf("(%s) %-8s %3d obs, %2d nights: failed: %s\n", number, label, length(arc), nn, sprint(showerror, e)[1:min(end, 150)])
        end
        flush(stdout)
    end
end
