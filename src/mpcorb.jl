"""
    read_mpcorb(path) -> Vector{NamedTuple}

Every orbit in a Minor Planet Center `MPCORB.DAT` file (the one Astrometrica
downloads into its `Catalogs` directory works): `packed`, `name`, `H`, `G`,
`epoch` (Julian Date, TT), mean anomaly `M`, argument of perihelion `peri`,
ascending node `node`, inclination `incl` (degrees, J2000 ecliptic), `e`,
semi-major axis `a` (AU) and mean daily motion `n` (degrees/day), parsed from
the MPC's documented fixed columns. Header lines and lines too short to hold
an orbit are skipped.
"""
function read_mpcorb(path::AbstractString)
    orbits = NamedTuple{(:packed, :name, :H, :G, :epoch, :M, :peri, :node, :incl, :e, :a, :n),
                        Tuple{String,String,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}[]
    for line in eachline(path)
        length(line) >= 103 || continue
        epoch = _unpack_mpc_epoch(line[21:25])
        epoch === nothing && continue
        num(r) = tryparse(Float64, strip(line[r]))
        vals = (num(9:13), num(15:19), num(27:35), num(38:46), num(49:57), num(60:68), num(71:79), num(81:91), num(93:103))
        any(isnothing, vals[3:end]) && continue
        H = something(vals[1], NaN); G = something(vals[2], 0.15)
        name = length(line) >= 194 ? strip(line[167:194]) : strip(line[1:7])
        push!(orbits, (packed=strip(line[1:7]), name=String(name), H=H, G=G, epoch=epoch, M=vals[3], peri=vals[4],
                       node=vals[5], incl=vals[6], e=vals[7], a=vals[9], n=vals[8]))
    end
    return orbits
end

# MPC packed date ("K2669" = 2026 June 9.0 TT) to a Julian Date; `nothing` if malformed.
function _unpack_mpc_epoch(s::AbstractString)
    length(s) == 5 || return nothing
    century = Dict('I' => 1800, 'J' => 1900, 'K' => 2000)
    haskey(century, s[1]) || return nothing
    digit(c) = isdigit(c) ? c - '0' : 'A' <= c <= 'V' ? c - 'A' + 10 : nothing
    y = tryparse(Int, s[2:3]); m = digit(s[4]); d = digit(s[5])
    (y === nothing || m === nothing || d === nothing) && return nothing
    return datetime2julian(DateTime(century[s[1]] + y, m, d))
end

const _OBLIQUITY_J2000 = 23.4392911                     # degrees
const _AU_PER_DAY_LIGHT = 173.144632674                 # speed of light, AU/day
const _TT_MINUS_UTC = (37 + 32.184) / 86400             # days; leap seconds since 2017

"""
    _heliocentric_equatorial(orbit, jd_tt) -> (x, y, z)

Two-body heliocentric position (AU, J2000 equatorial) of an MPCORB orbit at
`jd_tt`: mean anomaly advanced at the orbit's own mean motion, Kepler's
equation solved by Newton iteration. The same approximation Astrometrica's
known-object overlay relies on — fine over the weeks to months between an
MPCORB epoch and an observation, not for long arcs.
"""
function _heliocentric_equatorial(o, jd_tt::Real)
    M = deg2rad(mod(o.M + o.n * (jd_tt - o.epoch), 360))
    E = o.e < 0.8 ? M : π
    for _ in 1:30
        dE = (E - o.e * sin(E) - M) / (1 - o.e * cos(E))
        E -= dE
        abs(dE) < 1e-12 && break
    end
    xp = o.a * (cos(E) - o.e)
    yp = o.a * sqrt(1 - o.e^2) * sin(E)
    ω, Ω, i = deg2rad(o.peri), deg2rad(o.node), deg2rad(o.incl)
    # orbital plane -> ecliptic
    xe = (cos(ω) * cos(Ω) - sin(ω) * sin(Ω) * cos(i)) * xp + (-sin(ω) * cos(Ω) - cos(ω) * sin(Ω) * cos(i)) * yp
    ye = (cos(ω) * sin(Ω) + sin(ω) * cos(Ω) * cos(i)) * xp + (-sin(ω) * sin(Ω) + cos(ω) * cos(Ω) * cos(i)) * yp
    ze = sin(ω) * sin(i) * xp + cos(ω) * sin(i) * yp
    ε = deg2rad(_OBLIQUITY_J2000)
    return (xe, ye * cos(ε) - ze * sin(ε), ye * sin(ε) + ze * cos(ε))
end

# Horizons' answers (the `result` texts) to `query` at the Julian Dates
# `times`, already formatted, asked 50 at a time: a few hundred make a URL
# longer than it accepts (HTTP 502). If a list of times gets no table back
# (on 2026-11-09 Horizons misread every TLIST), each time is asked alone.
function _horizons(query::AbstractDict, times::AbstractVector{<:AbstractString})
    return map(Iterators.partition(times, 50)) do chunk
        result = _horizons_result(merge(query, Dict("TLIST" => "'$(join(chunk, " "))'", "TLIST_TYPE" => "JD")))
        _horizons_table(result) === nothing || return result
        _horizons_merge([_horizons_result(merge(query, _horizons_span(t))) for t in chunk])
    end
end

function _horizons_result(q::AbstractDict)
    body = _with_retries(_TAP_ATTEMPTS) do
        String(HTTP.get("https://ssd.jpl.nasa.gov/api/horizons.api"; query=merge(q, Dict("format" => "json")),
                        read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
    end
    return JSON.parse(body)["result"]
end

# One time as a span: Horizons lists the start, and a step longer than the span
# nothing after it.
function _horizons_span(time::AbstractString)
    jd = parse(Float64, time)
    return Dict("START_TIME" => "'JD $time'", "STOP_TIME" => "'JD $(@sprintf("%.8f", jd + 0.01))'",
                "STEP_SIZE" => "'1 d'")
end

# The rows between $$SOE and $$EOE, or `nothing` if there are none.
function _horizons_table(result::AbstractString)
    m = match(r"\$\$SOE\n(.*?)\$\$EOE"s, result)
    (m === nothing || isempty(strip(m[1]))) && return nothing
    return m[1]
end

# Answers to one time each, as one answer: the first one's text with every
# answer's rows, in order.
function _horizons_merge(results::AbstractVector{<:AbstractString})
    rows = map(results) do r
        t = _horizons_table(r)
        t === nothing && error("unexpected Horizons response: $(first(r, 300))")
        t
    end
    s = first(results)
    m = match(r"\$\$SOE\n.*?\$\$EOE"s, s)
    return string(SubString(s, 1, prevind(s, m.offset)), "\$\$SOE\n", join(rows), "\$\$EOE",
                  SubString(s, m.offset + ncodeunits(m.match)))
end

"""
    _observer_heliocentric(jds_utc, observatory) -> Vector{NTuple{3,Float64}}

The observer's heliocentric position (AU, ICRF/J2000 equatorial) at each UTC
Julian Date, from JPL Horizons (a live request): the Sun's vector as seen
from the MPC site `observatory`, negated. `observatory=nothing` gives the
Earth's centre instead.
"""
_observer_heliocentric(jds_utc::AbstractVector{<:Real}, observatory) =
    _observer_heliocentric_tt(jds_utc .+ _TT_MINUS_UTC, observatory)

# The same at TT Julian Dates.
function _observer_heliocentric_tt(jds_tt::AbstractVector{<:Real}, observatory)
    center = observatory === nothing ? "500@399" : "$(observatory)@399"
    query = Dict("COMMAND" => "'10'", "EPHEM_TYPE" => "VECTORS", "CENTER" => "'$center'", "TIME_TYPE" => "TT",
                 "REF_PLANE" => "FRAME", "REF_SYSTEM" => "ICRF", "VEC_TABLE" => "1", "OUT_UNITS" => "AU-D",
                 "OBJ_DATA" => "NO")
    positions = NTuple{3,Float64}[]
    for result in _horizons(query, [@sprintf("%.8f", jd) for jd in jds_tt])
        append!(positions, _horizons_vectors(result))
    end
    length(positions) == length(jds_tt) || error("Horizons returned $(length(positions)) vectors for $(length(jds_tt)) epochs")
    return positions
end

# The Sun as seen from the observer, negated, from a vector table in text form:
# on 2026-11-09 Horizons left the numbers out of its CSV vector tables.
function _horizons_vectors(result::AbstractString)
    table = _horizons_table(result)
    table === nothing && error("unexpected Horizons response: $(first(result, 300))")
    return [(-parse(Float64, m[1]), -parse(Float64, m[2]), -parse(Float64, m[3]))
            for m in eachmatch(r"X\s*=\s*(\S+)\s+Y\s*=\s*(\S+)\s+Z\s*=\s*(\S+)", table)]
end

# Topocentric astrometric (RA, Dec) in degrees, distance Δ and heliocentric r (AU)
# of an orbit seen from `obs` at `jd_utc`, with one light-time iteration.
function _mpcorb_radec(o, jd_utc::Real, obs::NTuple{3,Float64})
    t = jd_utc + _TT_MINUS_UTC
    p = _heliocentric_equatorial(o, t)
    Δ = hypot(p[1] - obs[1], p[2] - obs[2], p[3] - obs[3])
    p = _heliocentric_equatorial(o, t - Δ / _AU_PER_DAY_LIGHT)
    v = (p[1] - obs[1], p[2] - obs[2], p[3] - obs[3])
    Δ = hypot(v...)
    return (ra=mod(atand(v[2], v[1]), 360), dec=asind(v[3] / Δ), Δ=Δ, r=hypot(p...))
end

# H-G apparent V magnitude at heliocentric r, geocentric Δ and phase angle.
function _hg_magnitude(H, G, r, Δ, obs_r)
    cosα = clamp((r^2 + Δ^2 - obs_r^2) / (2r * Δ), -1, 1)
    return H + 5log10(r * Δ) + _hg_phase_dimming(acos(cosα), G)
end

# The H-G phase function's dimming (magnitudes) at phase angle α (radians).
function _hg_phase_dimming(α, G)
    φ1 = exp(-3.33 * tan(α / 2)^0.63)
    φ2 = exp(-1.87 * tan(α / 2)^1.22)
    return -2.5log10((1 - G) * φ1 + G * φ2)
end

"""
    _crossmatch_mpcorb(candidates, radius; mpcorb, observatory) -> Table

[`crossmatch_catalog`](@ref)'s `:mpcorb` backend: every row of `candidates`
against the orbits in `mpcorb` (a path to `MPCORB.DAT`, or what
[`read_mpcorb`](@ref) returned) propagated to the row's own epoch, within
`radius` arcsec. Orbits are first screened at one epoch to the ones within a
degree of any candidate, so only those are propagated per row.

Next to SkyBoT because the MPC's current orbits identify recently linked
or poorly observed objects that older orbits place degrees away.
"""
function _crossmatch_mpcorb(candidates, radius::Real; mpcorb, observatory=nothing)
    rows = collect(candidates)
    out = (id=Int[], name=String[], ra=Float64[], dec=Float64[], class=String[], mv=Float64[], distance_arcsec=Float64[])
    isempty(rows) && return Table(; out...)
    orbits = mpcorb isa AbstractString ? read_mpcorb(mpcorb) : mpcorb
    epochs = sort(unique(r.epoch for r in rows))
    observers = Dict(zip(epochs, _observer_heliocentric(epochs, observatory)))

    # screen: orbits within 1 degree of any candidate at the first epoch
    t0 = epochs[1]; obs0 = observers[t0]
    near = filter(orbits) do o
        p = _mpcorb_radec(o, t0, obs0)
        any(r -> _angular_distance_arcsec(p.ra, p.dec, r.ra, r.dec) < 3600, rows)
    end

    for r in rows, o in near
        obs = observers[r.epoch]
        p = _mpcorb_radec(o, r.epoch, obs)
        d = _angular_distance_arcsec(p.ra, p.dec, r.ra, r.dec)
        d <= radius || continue
        push!(out.id, r.id); push!(out.name, o.name); push!(out.ra, p.ra); push!(out.dec, p.dec)
        push!(out.class, ""); push!(out.mv, _hg_magnitude(o.H, o.G, p.r, p.Δ, hypot(obs...)))
        push!(out.distance_arcsec, d)
    end
    return Table(; out...)
end
