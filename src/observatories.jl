#=
Which observatory took a set of frames, from the frames themselves: FITS
headers carry the telescope's name and, often, its site's coordinates — but
not its MPC code, which every report starts with and which places
"known objects" as seen from the right spot on Earth. The site narrows the
MPC's list of observatories down to a few kilometres; the names tell apart
telescopes that share a mountain (Pan-STARRS 1 and 2, F51 and F52, ~10 m
apart; the three Palomar codes).
=#

const OBSCODES_URL = "https://minorplanetcenter.net/iau/lists/ObsCodes.html"

"""
    read_obscodes(path) -> Vector{NamedTuple}

The MPC's list of observatory codes (`ObsCodes.html`, the HTML page or its
text): `code`, east longitude `lon` (degrees), parallax constants `rho_cos`
and `rho_sin` (Earth radii), and `name`, from its fixed columns. Codes
without a site on Earth (spacecraft, `"250 Hubble Space Telescope"`) are
skipped.
"""
function read_obscodes(path::AbstractString)
    codes = NamedTuple{(:code, :lon, :rho_cos, :rho_sin, :name),Tuple{String,Float64,Float64,Float64,String}}[]
    for line in eachline(path)
        (length(line) >= 31 && isascii(line) && occursin(r"^[0-9A-Z]{3} ", line)) || continue
        vals = tryparse.(Float64, strip.((line[5:13], line[14:21], line[22:30])))
        any(isnothing, vals) && continue
        push!(codes, (code=line[1:3], lon=vals[1], rho_cos=vals[2], rho_sin=vals[3], name=String(strip(line[31:end]))))
    end
    return codes
end

"""
    download_obscodes(path=<cache>/ObsCodes.txt; max_age_hours=720, url=OBSCODES_URL) -> String

A copy of the MPC's observatory list at `path` (default next to
[`default_mpcorb_path`](@ref)'s), fetched if missing or older than
`max_age_hours` (30 days: codes are added, rarely changed). On a failed
download an existing copy is kept and returned, with a warning.
"""
function download_obscodes(path::AbstractString=joinpath(dirname(default_mpcorb_path()), "ObsCodes.txt");
                           max_age_hours::Real=720, url::AbstractString=OBSCODES_URL)
    isfile(path) && (time() - mtime(path)) < max_age_hours * 3600 && return path
    part = path * ".part"
    try
        mkpath(dirname(path))
        text = String(HTTP.get(url; read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
        write(part, replace(text, r"</?pre>" => ""))
        n = length(read_obscodes(part))
        n > 1000 || error("downloaded observatory list has only $n sites")
        mv(part, path; force=true)
    catch e
        rm(part; force=true)
        isfile(path) || rethrow()
        @warn "could not refresh the MPC's observatory list; using the existing copy" path exception=e
    end
    return path
end

# Site keywords seen in real headers, and their conventions. Longitudes are
# turned east-positive degrees.
const _LON_KEYS = (("OBSLON", :east), ("LONG-OBS", :east), ("OBS-LONG", :east), ("SITELONG", :east),
                   ("FPA.LONGITUDE", :west_hours))   # Pan-STARRS: west longitude, in hours
const _LAT_KEYS = ("OBSLAT", "LAT-OBS", "OBS-LAT", "SITELAT", "FPA.LATITUDE")
const _ELEV_KEYS = ("OBSALT", "ALT-OBS", "OBS-ELEV", "SITEELEV", "FPA.ELEVATION")

# A number, or a sexagesimal "-116 51 36" / "-116:51:36" string, in its own unit.
function _header_angle(v)
    v isa Real && return Float64(v)
    v isa AbstractString || return nothing
    m = match(r"^\s*([+-]?)(\d+(?:\.\d*)?)(?:[\s:]+(\d+(?:\.\d*)?))?(?:[\s:]+(\d+(?:\.\d*)?))?\s*$", v)
    m === nothing && return nothing
    x = parse(Float64, m[2]) + (m[3] === nothing ? 0.0 : parse(Float64, m[3]) / 60) +
        (m[4] === nothing ? 0.0 : parse(Float64, m[4]) / 3600)
    return m[1] == "-" ? -x : x
end

_header_text(h, key) = haskey(h, key) && h[key] isa AbstractString ? String(strip(h[key])) : ""

"""
    frame_site(header::FITSHeader) -> NamedTuple

What a frame's header says about where it was taken: `telescope`
(`TELESCOP`), `instrument` (`INSTRUME`), `origin` (`ORIGIN`), `observatory`
(`OBSERVAT`), and the site's east longitude `lon` and latitude `lat`
(degrees) and elevation `elev` (metres) — `nothing` for any of the three the
header doesn't give. Reads the usual keywords (`OBSLON`/`OBSLAT`/`OBSALT` as
ZTF writes them, `SITELONG`/`SITELAT`, `LONG-OBS`...) and Pan-STARRS's
`FPA.LONGITUDE` (west, in hours) / `FPA.LATITUDE` / `FPA.ELEVATION`.
"""
function frame_site(h::FITSHeader)
    lon = nothing
    for (key, convention) in _LON_KEYS
        haskey(h, key) || continue
        v = _header_angle(h[key])
        v === nothing && continue
        lon = convention === :west_hours ? -15v : v
        break
    end
    lat = nothing
    for key in _LAT_KEYS
        haskey(h, key) && (lat = _header_angle(h[key])) !== nothing && break
    end
    elev = nothing
    for key in _ELEV_KEYS
        haskey(h, key) && h[key] isa Real && (elev = Float64(h[key]); break)
    end
    lon === nothing || (lon = mod(lon + 180, 360) - 180)
    return (telescope=_header_text(h, "TELESCOP"), instrument=_header_text(h, "INSTRUME"),
            origin=_header_text(h, "ORIGIN"), observatory=_header_text(h, "OBSERVAT"), lon=lon, lat=lat, elev=elev)
end

# Geocentric position (km) of a site from geodetic latitude/elevation (WGS84),
# or of an MPC code from its parallax constants.
const _EARTH_RADIUS_KM = 6378.137
function _site_xyz(lon, lat, elev_m)
    f = 1 / 298.257223563
    C = 1 / sqrt(cosd(lat)^2 + (1 - f)^2 * sind(lat)^2)
    S = (1 - f)^2 * C
    h = something(elev_m, 0.0) / 1000
    rc, rs = (_EARTH_RADIUS_KM * C + h) * cosd(lat), (_EARTH_RADIUS_KM * S + h) * sind(lat)
    return (rc * cosd(lon), rc * sind(lon), rs)
end
_code_xyz(c) = _EARTH_RADIUS_KM .* (c.rho_cos * cosd(c.lon), c.rho_cos * sind(c.lon), c.rho_sin)

const _STOP_WORDS = Set(["the", "of", "and", "de", "del", "la", "el", "observatory", "observatorio", "obs",
                         "telescope", "station", "mount", "mountain", "inch", "m", "cm"])

# Words of a name and their acronym ("Pan-STARRS 1" -> pan, starrs, 1, ps1;
# "Zwicky Transient Facility" -> ztf), lower case, stop words left out.
function _name_tokens(s::AbstractString)
    out = Set{String}()
    for part in split(lowercase(s), r"[,;/()]|--")
        words = [w for w in split(part, r"[^a-z0-9]+") if !isempty(w)]
        union!(out, filter(w -> !(w in _STOP_WORDS) && (length(w) > 1 || isdigit(w[1])), words))
        letters = [all(isdigit, w) ? w : string(w[1]) for w in words if !(w in _STOP_WORDS)]
        length(letters) >= 2 && push!(out, join(letters))
    end
    return out
end

"""
    observatory_candidates(site, codes; max_km=10) -> Vector{NamedTuple}

MPC codes that could be the observatory [`frame_site`](@ref) describes,
best first: `code`, `name`, `distance_km` from the header's site (`NaN` if
it has none) and `score`, the number of words or acronyms the code's name
shares with the header's telescope, instrument, origin and observatory
names. With a site, only codes within `max_km` are considered, ranked by
score, then distance; without one, every code sharing a word, by score.
`codes` is what [`read_obscodes`](@ref) returns.
"""
function observatory_candidates(site, codes; max_km::Real=10)
    words = union(_name_tokens.((site.telescope, site.instrument, site.origin, site.observatory))...)
    has_site = site.lon !== nothing && site.lat !== nothing
    here = has_site ? _site_xyz(site.lon, site.lat, site.elev) : nothing
    out = NamedTuple{(:code, :name, :distance_km, :score),Tuple{String,String,Float64,Int}}[]
    for c in codes
        score = length(intersect(words, _name_tokens(c.name)))
        d = has_site ? hypot((_code_xyz(c) .- here)...) : NaN
        (has_site ? d <= max_km : score > 0) || continue
        push!(out, (code=c.code, name=c.name, distance_km=d, score=score))
    end
    sort!(out; by=c -> (-c.score, isnan(c.distance_km) ? 0.0 : c.distance_km))
    return out
end
