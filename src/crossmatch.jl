const _SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
const _VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

const _SKYBOT_URL = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"

"""
    crossmatch_catalog(candidates, catalog::Symbol; radius::Real)

Cross-match candidate tracklets against a known-object catalog
(`:skybot`, `:vsx`, or `:simbad`) within `radius` arcseconds, to separate
known objects from candidates warranting human verification.

`candidates` is an iterable of rows with `id`, `ra`, `dec` fields (degrees,
J2000), such as the table returned by [`astrometric_calibrate`](@ref).
SkyBoT additionally requires an `epoch` field (Julian Date) on each row,
since it computes solar-system-object ephemerides for a specific instant
rather than querying a static catalog. `:vsx` and `:simbad` are each
queried directly against their own CDS TAP service (ADQL cone search),
**one request per candidate** — not the single batched request the CDS
X-Match service used to offer, since X-Match itself became unreliable
(extended, total outages; see `INVESTIGATION_LOG.md`) while the
underlying SIMBAD and VizieR TAP services stayed up. The trade-off is real
(N candidates means N requests, not 1) and matters for a large candidate
list — batch by querying a shared sky region directly via TAP if that
becomes a bottleneck; not done here since it wasn't yet.

Returns a table with one row per (candidate, catalog match) pair found
within `radius`; `:vsx` and `:simbad` return different columns (VSX
carries variability class/magnitude/period, SIMBAD doesn't), matching
what each catalog actually offers, but both always include `id`, `ra`,
`dec`, and `distance_arcsec`. Candidate `id`s absent from the returned
table have no known counterpart in `catalog`.
"""
function crossmatch_catalog(candidates, catalog::Symbol; radius::Real)
    if catalog === :skybot
        return _crossmatch_skybot(candidates, radius)
    elseif catalog === :simbad
        return _crossmatch_simbad(candidates, radius)
    elseif catalog === :vsx
        return _crossmatch_vsx(candidates, radius)
    else
        throw(ArgumentError("unknown catalog $(repr(catalog)); expected :skybot, :vsx, or :simbad"))
    end
end

"""
    _tap_query(url, query) -> CSV.File

Run `query` (an ADQL string) as a synchronous TAP query against `url`,
returning the CSV response parsed by `CSV.File`. Shared by
`_crossmatch_simbad` and `_crossmatch_vsx`.
"""
function _tap_query(url::AbstractString, query::AbstractString)
    body = HTTP.Form(Dict(
        "REQUEST" => "doQuery", "LANG" => "ADQL", "FORMAT" => "csv", "QUERY" => query))
    response = HTTP.post(url, [], body)
    return CSV.File(response.body)
end

"""
    _cds_cone_query(select, table, ra_col, dec_col, ra, dec, radius_deg) -> String

An ADQL synchronous cone-search query: rows of `table` within
`radius_deg` of `(ra, dec)`, plus a `distance_arcsec` column (via ADQL's
`DISTANCE`, in degrees, converted here) — shared by
`_crossmatch_simbad` and `_crossmatch_vsx`, which differ
only in `select`/`table`/coordinate column names.
"""
function _cds_cone_query(select::AbstractString, table::AbstractString,
                          ra_col::AbstractString, dec_col::AbstractString,
                          ra::Real, dec::Real, radius_deg::Real)
    point = "POINT('ICRS', $ra_col, $dec_col)"
    return "SELECT $select, DISTANCE($point, POINT('ICRS', $ra, $dec)) * 3600 AS distance_arcsec " *
           "FROM $table WHERE CONTAINS($point, CIRCLE('ICRS', $ra, $dec, $radius_deg)) = 1"
end

_nan_if_missing(x) = x === missing ? NaN : Float64(x)

function _crossmatch_simbad(candidates, radius::Real)
    0 < radius <= 180 || throw(ArgumentError("radius must be in (0, 180] arcsec"))
    radius_deg = radius / 3600

    id = Int[]; name = String[]; ra = Float64[]; dec = Float64[]; distance_arcsec = Float64[]
    for c in candidates
        query = _cds_cone_query("main_id, ra, dec", "basic", "ra", "dec", c.ra, c.dec, radius_deg)
        for row in _tap_query(_SIMBAD_TAP_URL, query)
            push!(id, c.id)
            push!(name, String(row.main_id))
            push!(ra, row.ra)
            push!(dec, row.dec)
            push!(distance_arcsec, row.distance_arcsec)
        end
    end

    return Table(; id, name, ra, dec, distance_arcsec)
end

function _crossmatch_vsx(candidates, radius::Real)
    0 < radius <= 180 || throw(ArgumentError("radius must be in (0, 180] arcsec"))
    radius_deg = radius / 3600

    id = Int[]; name = String[]; class = String[]; ra = Float64[]; dec = Float64[]
    mag_max = Float64[]; mag_min = Float64[]; period = Float64[]; distance_arcsec = Float64[]
    for c in candidates
        query = _cds_cone_query("Name, Type, RAJ2000, DEJ2000, max, min, Period", "\"B/vsx/vsx\"",
                                 "RAJ2000", "DEJ2000", c.ra, c.dec, radius_deg)
        for row in _tap_query(_VIZIER_TAP_URL, query)
            push!(id, c.id)
            # VizieR's TAP service returns Name/Type as fixed-width,
            # space-padded strings (e.g. "RS" arrives as "RS" followed by
            # 28 spaces) — found via a real query, not documented anywhere
            # obvious; strip or every string comparison against these
            # silently fails.
            push!(name, String(strip(row.Name)))
            push!(class, String(strip(row.Type)))
            push!(ra, row.RAJ2000)
            push!(dec, row.DEJ2000)
            push!(mag_max, _nan_if_missing(row.max))
            push!(mag_min, _nan_if_missing(row.min))
            push!(period, _nan_if_missing(row.Period))
            push!(distance_arcsec, row.distance_arcsec)
        end
    end

    return Table(; id, name, class, ra, dec, mag_max, mag_min, period, distance_arcsec)
end

function _crossmatch_skybot(candidates, radius::Real)
    radius_deg = radius / 3600

    id = Int[]
    name = String[]
    ra = Float64[]
    dec = Float64[]
    class = String[]
    mv = Float64[]
    distance_arcsec = Float64[]

    for c in candidates
        query = Dict(
            "-ra" => string(c.ra), "-dec" => string(c.dec), "-rd" => string(radius_deg),
            # Julian Dates are ~2.4e6, which Julia's default Float64 printing
            # renders in scientific notation (e.g. "2.4592886174e6"); SkyBoT
            # rejects that outright as an empty/null epoch, so every match
            # silently came back empty. @sprintf forces fixed-point.
            "-ep" => @sprintf("%.6f", c.epoch), "-mime" => "text", "-output" => "object",
        )
        response = HTTP.get(_SKYBOT_URL; query=query)
        for match in _parse_skybot(String(response.body))
            push!(id, c.id)
            push!(name, match.name)
            push!(ra, match.ra)
            push!(dec, match.dec)
            push!(class, match.class)
            push!(mv, match.mv)
            push!(distance_arcsec, match.distance_arcsec)
        end
    end

    return Table(; id, name, ra, dec, class, mv, distance_arcsec)
end

function _parse_skybot(text::AbstractString)
    matches = NamedTuple[]
    for line in eachline(IOBuffer(text))
        (isempty(line) || startswith(line, '#')) && continue
        fields = strip.(split(line, '|'))
        length(fields) < 8 && continue
        push!(matches, (
            name=String(fields[2]),
            ra=_skybot_hms_to_deg(fields[3]),
            dec=_skybot_dms_to_deg(fields[4]),
            class=String(fields[5]),
            mv=parse(Float64, fields[6]),
            distance_arcsec=parse(Float64, fields[8]),
        ))
    end
    return matches
end

function _skybot_hms_to_deg(s::AbstractString)
    h, m, sec = parse.(Float64, split(s))
    return 15 * (h + m / 60 + sec / 3600)
end

function _skybot_dms_to_deg(s::AbstractString)
    sign = startswith(s, '-') ? -1.0 : 1.0
    d, m, sec = parse.(Float64, split(lstrip(s, ['+', '-'])))
    return sign * (d + m / 60 + sec / 3600)
end
