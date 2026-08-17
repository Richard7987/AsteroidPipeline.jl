const _SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
const _VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

const _SKYBOT_URL = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"

# Candidates per TAP request for :vsx/:simbad. Conservative, not measured
# at higher N on these free, shared, anonymous-use CDS services — kept
# well short of any observed limit rather than pushed to one.
const _CDS_BATCH_SIZE = 50

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
batched `_CDS_BATCH_SIZE` (50) candidates per request instead of the
single batched request the old CDS X-Match service used to offer
(replaced since X-Match itself became unreliable — extended, total
outages; see the [Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#Batching-crossmatch_catalog...;-:vsx/...;-:simbad) — while the underlying SIMBAD and
VizieR TAP services stayed up). A real TAP `UPLOAD`-based batch (one
request for *any* N) was tried first and rejected: verified directly
against both services that it doesn't work simply here — SIMBAD accepts
an uploaded CSV but fails to resolve its columns, VizieR rejects anything
that isn't a full VOTable document — and ADQL `UNION` (the other
obvious batching route) is rejected outright by SIMBAD's parser. What
*does* work, confirmed against a real multi-candidate query: `OR`-chaining
one `CONTAINS(...)=1` clause per candidate in a single request, then
resolving which candidate each returned row matches client-side (no
`UNION`/per-row tag available in one `OR`-chained query, so each
returned row's distance is checked against every candidate in that
batch). Cuts N requests to `ceil(N / 50)`.

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

Retries once on `HTTP.ParseError` ("unexpected EOF while reading HTTP/1
data"): confirmed, via repeated direct `curl` requests against the exact
query that triggered it, to be a connection-reuse quirk on our end, not a
real SIMBAD/VizieR outage — the service itself answered the same query
successfully every time it was tried directly. A fresh connection (a new
request, not a retried read of the same one) resolves it in practice.
Any other exception, or a second `HTTP.ParseError`, still propagates —
this is a targeted retry for one confirmed-transient failure mode, not a
general-purpose retry loop.
"""
function _tap_query(url::AbstractString, query::AbstractString)
    make_body() = HTTP.Form(Dict(
        "REQUEST" => "doQuery", "LANG" => "ADQL", "FORMAT" => "csv", "QUERY" => query))
    try
        response = HTTP.post(url, [], make_body())
        return CSV.File(response.body)
    catch e
        e isa HTTP.ParseError || rethrow()
        response = HTTP.post(url, [], make_body())
        return CSV.File(response.body)
    end
end

"""
    _cds_batch_query(select, table, ra_col, dec_col, batch, radius_deg) -> String

An ADQL synchronous query matching rows of `table` against *any* member
of `batch` (an iterable of rows with `ra`/`dec` fields): one
`CONTAINS(...) = 1` cone-search clause per candidate, `OR`-chained into a
single request — TAP's `UNION` is unsupported here (confirmed against a
real query) and table `UPLOAD` doesn't work simply on these services
either (confirmed too — see [`crossmatch_catalog`](@ref)'s docstring),
so this is the batching approach that's actually verified to work.
Shared by `_crossmatch_simbad` and `_crossmatch_vsx`, which differ only
in `select`/`table`/coordinate column names. Since a single `OR`-chained
query can't tag which candidate a row matched, that's resolved
separately, client-side, after the query runs.
"""
function _cds_batch_query(select::AbstractString, table::AbstractString,
                           ra_col::AbstractString, dec_col::AbstractString,
                           batch, radius_deg::Real)
    point = "POINT('ICRS', $ra_col, $dec_col)"
    conditions = ["CONTAINS($point, CIRCLE('ICRS', $(c.ra), $(c.dec), $radius_deg)) = 1" for c in batch]
    return "SELECT $select FROM $table WHERE " * join(conditions, " OR ")
end

"""
    _angular_distance_arcsec(ra1, dec1, ra2, dec2) -> Float64

Great-circle angular distance (haversine formula — accurate at any
declination, including near the poles, unlike a flat small-angle
approximation) between two `(ra, dec)` points in degrees, returned in
arcseconds. Used to resolve which candidate in a
[`_cds_batch_query`](@ref) batch each returned row actually matches,
since ADQL's own `DISTANCE` needs one fixed reference point per query and
a batch has several.
"""
function _angular_distance_arcsec(ra1::Real, dec1::Real, ra2::Real, dec2::Real)
    deg2rad = pi / 180
    phi1, phi2 = dec1 * deg2rad, dec2 * deg2rad
    dphi = (dec2 - dec1) * deg2rad
    dlambda = (ra2 - ra1) * deg2rad
    a = sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
    return 2 * asin(min(1.0, sqrt(a))) / deg2rad * 3600
end

_nan_if_missing(x) = x === missing ? NaN : Float64(x)

function _crossmatch_simbad(candidates, radius::Real)
    0 < radius <= 180 || throw(ArgumentError("radius must be in (0, 180] arcsec"))
    radius_deg = radius / 3600

    id = Int[]; name = String[]; ra = Float64[]; dec = Float64[]; distance_arcsec = Float64[]
    for batch in Iterators.partition(collect(candidates), _CDS_BATCH_SIZE)
        query = _cds_batch_query("main_id, ra, dec", "basic", "ra", "dec", batch, radius_deg)
        for row in _tap_query(_SIMBAD_TAP_URL, query)
            for c in batch
                d = _angular_distance_arcsec(c.ra, c.dec, row.ra, row.dec)
                d > radius && continue
                push!(id, c.id)
                push!(name, String(row.main_id))
                push!(ra, row.ra)
                push!(dec, row.dec)
                push!(distance_arcsec, d)
            end
        end
    end

    return Table(; id, name, ra, dec, distance_arcsec)
end

function _crossmatch_vsx(candidates, radius::Real)
    0 < radius <= 180 || throw(ArgumentError("radius must be in (0, 180] arcsec"))
    radius_deg = radius / 3600

    id = Int[]; name = String[]; class = String[]; ra = Float64[]; dec = Float64[]
    mag_max = Float64[]; mag_min = Float64[]; period = Float64[]; distance_arcsec = Float64[]
    for batch in Iterators.partition(collect(candidates), _CDS_BATCH_SIZE)
        query = _cds_batch_query("Name, Type, RAJ2000, DEJ2000, max, min, Period", "\"B/vsx/vsx\"",
                                  "RAJ2000", "DEJ2000", batch, radius_deg)
        for row in _tap_query(_VIZIER_TAP_URL, query)
            # VizieR's TAP service returns Name/Type as fixed-width,
            # space-padded strings (e.g. "RS" arrives as "RS" followed by
            # 28 spaces) — found via a real query, not documented anywhere
            # obvious; strip or every string comparison against these
            # silently fails.
            row_ra, row_dec = row.RAJ2000, row.DEJ2000
            for c in batch
                d = _angular_distance_arcsec(c.ra, c.dec, row_ra, row_dec)
                d > radius && continue
                push!(id, c.id)
                push!(name, String(strip(row.Name)))
                push!(class, String(strip(row.Type)))
                push!(ra, row_ra)
                push!(dec, row_dec)
                push!(mag_max, _nan_if_missing(row.max))
                push!(mag_min, _nan_if_missing(row.min))
                push!(period, _nan_if_missing(row.Period))
                push!(distance_arcsec, d)
            end
        end
    end

    return Table(; id, name, class, ra, dec, mag_max, mag_min, period, distance_arcsec)
end

"""
    _skybot_matches(c, radius_deg) -> Vector{NamedTuple}

One SkyBoT cone-search request for a single candidate `c`. Factored out
of [`_crossmatch_skybot`](@ref) so it can be run concurrently, one task
per candidate — see that function's docstring for why.

Retries once on any `HTTP.HTTPError` (covers `HTTP.ConnectError`,
`HTTP.ParseError`, etc.): confirmed real on a long real crossmatch run
(thousands of candidates, several thousand real SkyBoT requests) — a
"tls write failed: connection is closed" error killed the whole run
partway through a large candidate list, after running cleanly for over
two hours. A fresh connection on retry is enough in practice; a second
failure still propagates.
"""
function _skybot_matches(c, radius_deg::Real)
    query = Dict(
        "-ra" => string(c.ra), "-dec" => string(c.dec), "-rd" => string(radius_deg),
        # Julian Dates are ~2.4e6, which Julia's default Float64 printing
        # renders in scientific notation (e.g. "2.4592886174e6"); SkyBoT
        # rejects that outright as an empty/null epoch, so every match
        # silently came back empty. @sprintf forces fixed-point.
        "-ep" => @sprintf("%.6f", c.epoch), "-mime" => "text", "-output" => "object",
    )
    try
        response = HTTP.get(_SKYBOT_URL; query=query)
        return _parse_skybot(String(response.body))
    catch e
        e isa HTTP.HTTPError || rethrow()
        response = HTTP.get(_SKYBOT_URL; query=query)
        return _parse_skybot(String(response.body))
    end
end

# Concurrent SkyBoT requests per crossmatch_catalog call. SkyBoT (unlike
# :vsx/:simbad's CDS TAP services, batched in a single request — see
# crossmatch_catalog's docstring) has no batch mode, only a one-candidate
# cone search; measured directly, real candidate lists (500+ tracklets
# from one IASC practice field) took ~9 minutes fully sequential, almost
# entirely spent waiting on network round trips, not computing anything.
# Benchmarked directly against the live service (20 real, identical
# requests): concurrency=8 gave a real 2.6x speedup (12.5s vs 32.8s
# sequential, same results); concurrency up to 60 ran clean with zero
# errors, though per-wave latency stopped shrinking much past ~20-40
# (IMCCE's own server-side queueing, not our bottleneck by then). Settled
# on 20 — comfortably inside the tested-clean range, not pushed to it, in
# the same spirit as `_CDS_BATCH_SIZE`'s conservative choice.
const _SKYBOT_CONCURRENCY = 20

"""
    _crossmatch_skybot(candidates, radius) -> Table

Query SkyBoT once per candidate, concurrently (up to
`_SKYBOT_CONCURRENCY` requests in flight at a time via Julia
`Task`s, not extra threads — this is a network-latency-bound workload,
not a CPU-bound one, so cooperative concurrency on however many threads
Julia was started with is enough to see the full speedup). Order of
results does not depend on request completion order: each candidate's
matches are collected into their own slot and the final table is built
from those slots in `candidates`' original order, not arrival order.
"""
function _crossmatch_skybot(candidates, radius::Real)
    radius_deg = radius / 3600
    candidates_vec = collect(candidates)
    per_candidate = Vector{Vector{NamedTuple}}(undef, length(candidates_vec))

    semaphore = Base.Semaphore(_SKYBOT_CONCURRENCY)
    @sync for (i, c) in enumerate(candidates_vec)
        @async begin
            Base.acquire(semaphore)
            try
                per_candidate[i] = _skybot_matches(c, radius_deg)
            finally
                Base.release(semaphore)
            end
        end
    end

    id = Int[]
    name = String[]
    ra = Float64[]
    dec = Float64[]
    class = String[]
    mv = Float64[]
    distance_arcsec = Float64[]
    for (c, matches) in zip(candidates_vec, per_candidate)
        for match in matches
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
