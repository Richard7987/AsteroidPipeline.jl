const _CDS_XMATCH_URL = "https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
const _CDS_CATALOG_NAMES = Dict(:vsx => "vizier:B/vsx/vsx", :simbad => "simbad")

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
rather than querying a static catalog; `:vsx` and `:simbad` are queried in
a single batched request via the CDS X-Match service.

Returns a table with one row per (candidate, catalog match) pair found
within `radius`. Candidate `id`s absent from the returned table have no
known counterpart in `catalog`.
"""
function crossmatch_catalog(candidates, catalog::Symbol; radius::Real)
    if catalog === :skybot
        return _crossmatch_skybot(candidates, radius)
    elseif haskey(_CDS_CATALOG_NAMES, catalog)
        return _crossmatch_cds(candidates, _CDS_CATALOG_NAMES[catalog], radius)
    else
        throw(ArgumentError("unknown catalog $(repr(catalog)); expected :skybot, :vsx, or :simbad"))
    end
end

function _crossmatch_cds(candidates, cds_name::AbstractString, radius::Real)
    0 < radius <= 180 ||
        throw(ArgumentError("radius must be in (0, 180] arcsec for the CDS X-Match service"))

    upload = IOBuffer()
    println(upload, "id,ra,dec")
    for c in candidates
        println(upload, "$(c.id),$(c.ra),$(c.dec)")
    end
    seekstart(upload)

    body = HTTP.Form(Dict(
        "request" => "xmatch",
        "distMaxArcsec" => string(radius),
        "RESPONSEFORMAT" => "csv",
        "cat1" => HTTP.Multipart("candidates.csv", upload, "text/csv"),
        "colRA1" => "ra",
        "colDec1" => "dec",
        "cat2" => cds_name,
    ))
    response = HTTP.post(_CDS_XMATCH_URL, [], body)

    return Table(CSV.File(response.body))
end

function _crossmatch_skybot(candidates, radius::Real)
    radius_deg = radius / 3600
    rows = NamedTuple[]

    for c in candidates
        query = Dict(
            "-ra" => string(c.ra), "-dec" => string(c.dec), "-rd" => string(radius_deg),
            "-ep" => string(c.epoch), "-mime" => "text", "-output" => "object",
        )
        response = HTTP.get(_SKYBOT_URL; query=query)
        for match in _parse_skybot(String(response.body))
            push!(rows, merge((id=c.id,), match))
        end
    end

    return Table(rows)
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
