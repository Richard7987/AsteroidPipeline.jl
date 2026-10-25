const _SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
const _VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

const _SKYBOT_URL = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"

# Candidates per TAP request for :vsx/:simbad: kept small on these shared,
# anonymous CDS services, not measured against their limit.
const _CDS_BATCH_SIZE = 50

"""
    crossmatch_catalog(candidates, catalog::Symbol; radius::Real,
                       observatory::Union{Nothing,AbstractString}=nothing, mpcorb=nothing)

Cross-match candidate tracklets against a known-object catalog
(`:skybot`, `:mpcorb`, `:vsx` or `:simbad`) within `radius` arcseconds, to
separate known objects from candidates warranting human verification.

`candidates` is an iterable of rows with `id`, `ra`, `dec` fields (degrees,
J2000), such as the table returned by [`astrometric_calibrate`](@ref);
`:skybot` and `:mpcorb` also need each row's `epoch` (Julian Date).

- `:skybot` asks IMCCE's SkyBoT for the solar-system objects at each row's
  epoch.
- `:mpcorb` propagates the Minor Planet Center's orbits in an `MPCORB.DAT`
  file (`mpcorb`: its path, or [`read_mpcorb`](@ref)'s result) to each
  row's epoch (see [`_crossmatch_mpcorb`](@ref)); same columns as `:skybot`.
- `:vsx` and `:simbad` query each CDS TAP service with one `OR`-chained
  ADQL cone per candidate, `_CDS_BATCH_SIZE` (50) per request, and assign
  each returned row to its candidates client-side. (TAP `UPLOAD` and ADQL
  `UNION` don't work on these services.) VSX rows carry variability class,
  magnitudes and period.

Every result has `id`, `ra`, `dec` and `distance_arcsec`, one row per
(candidate, match) pair; a candidate `id` absent from it has no known
counterpart in `catalog`.

For `:skybot` and `:mpcorb`, `observatory` (an MPC code, e.g. `"F51"` for
Pan-STARRS1, `"I41"` for ZTF) computes each position from that site rather
than the Earth's centre. Leave it unset only for distant objects: parallax
moves a near one by many arcseconds (~15" for a Mars-crosser seen from F51).

Each row is matched at its own `epoch`, so passing *every* row of a
tracklet and counting matches per `(id, name)` checks that the tracklet
follows the object frame by frame. One matched point is a weak test: a
static star near a known asteroid's track matches once.
"""
function crossmatch_catalog(candidates, catalog::Symbol; radius::Real,
                            observatory::Union{Nothing,AbstractString}=nothing, mpcorb=nothing)
    observatory === nothing || catalog in (:skybot, :mpcorb) ||
        throw(ArgumentError("observatory only applies to :skybot and :mpcorb"))
    observatory === nothing || length(observatory) == 3 ||
        throw(ArgumentError("observatory must be a 3-character MPC observatory code"))
    if catalog === :skybot
        return _crossmatch_skybot(candidates, radius; observatory)
    elseif catalog === :mpcorb
        mpcorb === nothing && throw(ArgumentError(":mpcorb needs `mpcorb`, a path to MPCORB.DAT or read_mpcorb's result"))
        return _crossmatch_mpcorb(candidates, radius; mpcorb, observatory)
    elseif catalog === :simbad
        return _crossmatch_simbad(candidates, radius)
    elseif catalog === :vsx
        return _crossmatch_vsx(candidates, radius)
    else
        throw(ArgumentError("unknown catalog $(repr(catalog)); expected :skybot, :mpcorb, :vsx, or :simbad"))
    end
end

"""
    _tap_query(url, query) -> CSV.File

Run `query` (an ADQL string) as a synchronous TAP query against `url`,
returning the CSV response parsed by `CSV.File`. Shared by
`_crossmatch_simbad`, `_crossmatch_vsx` and [`gaia_reference_stars`](@ref).

A request that stops receiving data for `_CATALOG_IDLE_TIMEOUT` seconds
is abandoned (`HTTP.TimeoutError`) instead of waiting forever. Retried
`_TAP_ATTEMPTS` tries in all on a failure of the service itself (see
[`_is_service_failure`](@ref)), never on a 4xx answer.
"""
function _tap_query(url::AbstractString, query::AbstractString)
    make_body() = HTTP.Form(Dict(
        "REQUEST" => "doQuery", "LANG" => "ADQL", "FORMAT" => "csv", "QUERY" => query))
    return _with_retries(_TAP_ATTEMPTS) do
        CSV.File(HTTP.post(url, [], make_body(); read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
    end
end

# `f()`, tried up to `attempts` times while `retryable(e)`, waiting `backoff`
# seconds and doubling it each time; the last failure propagates.
function _with_retries(f, attempts::Integer; retryable=_is_service_failure, backoff::Real=2.0)
    for attempt in 1:attempts
        try
            return f()
        catch e
            (retryable(e) && attempt < attempts) || rethrow()
            sleep(backoff * 2.0^(attempt - 1))
        end
    end
end

const _TAP_ATTEMPTS = 3

"""
    _is_service_failure(e) -> Bool

Whether `e` is the catalog *service* failing — unreachable, reset, TLS
or protocol trouble, a timeout, a 5xx status — which a retry (or another
copy of the catalog) can get past, as opposed to a 4xx status, which
means the query itself is wrong and would fail anywhere.
"""
_is_service_failure(e) = e isa Base.IOError ||
                         (e isa HTTP.HTTPError && !(e isa HTTP.StatusError && e.status < 500))

# Seconds a catalog request may go without receiving data: an idle limit,
# so a large result that keeps streaming is never cut off.
const _CATALOG_IDLE_TIMEOUT = 60

"""
    _cds_batch_query(select, table, ra_col, dec_col, batch, radius_deg) -> String

An ADQL synchronous query matching rows of `table` against *any* member
of `batch` (an iterable of rows with `ra`/`dec` fields): one
`CONTAINS(...) = 1` cone-search clause per candidate, `OR`-chained into a
single request (see [`crossmatch_catalog`](@ref)). Shared by
`_crossmatch_simbad` and `_crossmatch_vsx`; which candidate a row matched
is resolved client-side.
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
            # VizieR pads Name/Type with spaces; unstripped, every comparison fails.
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

One SkyBoT cone-search request for a single candidate `c`, run one task
per candidate by [`_crossmatch_skybot`](@ref).

Up to `_SKYBOT_ATTEMPTS` tries on any `HTTP.HTTPError` or
[`SkyBoTServiceError`](@ref) (SkyBoT failing server-side behind HTTP 200,
which clears on a resend). `fetch` is the HTTP call itself, injectable so
the retries can be tested without the network.

The cone asked of SkyBoT is never smaller than `_SKYBOT_MIN_CONE` (2′),
and matches beyond `radius_deg` are dropped here: asked for 5″, SkyBoT
missed a known asteroid 0.94″ from the centre in two of four ZTF frames.
"""
function _skybot_matches(c, radius_deg::Real;
                         observatory::Union{Nothing,AbstractString}=nothing,
                         fetch=query -> String(HTTP.get(_SKYBOT_URL; query=query,
                                                       read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body),
                         backoff_seconds::Real=1.0)
    query = Dict(
        "-ra" => string(c.ra), "-dec" => string(c.dec), "-rd" => string(max(radius_deg, _SKYBOT_MIN_CONE)),
        # SkyBoT reads a JD in scientific notation ("2.4592886174e6") as no epoch.
        "-ep" => @sprintf("%.6f", c.epoch), "-mime" => "text", "-output" => "object",
    )
    observatory === nothing || (query["-loc"] = observatory)
    retryable(e) = e isa HTTP.HTTPError || e isa SkyBoTServiceError
    return _with_retries(_SKYBOT_ATTEMPTS; retryable, backoff=backoff_seconds) do
        filter(m -> m.distance_arcsec <= radius_deg * 3600, _parse_skybot(fetch(query)))
    end
end

# The smallest cone _skybot_matches asks for, in degrees (see its docstring).
const _SKYBOT_MIN_CONE = 2 / 60

# A request to an online service that failed: callers leave a warning instead.
_is_unavailable(e) = e isa HTTP.HTTPError || e isa Base.IOError || e isa SkyBoTServiceError

# The exception at the root of a failed @sync/@async block.
_task_root(e) = e isa CompositeException && !isempty(e.exceptions) ? _task_root(first(e.exceptions)) :
                e isa TaskFailedException ? _task_root(e.task.exception) : e

# Tries per candidate in _skybot_matches (1/2/4 s backoff): enough for a short
# burst of server errors, without stalling long on a real outage.
const _SKYBOT_ATTEMPTS = 4

"""
    SkyBoTServiceError(message)

SkyBoT answered, but with a server-side failure instead of a result — a
`# Flag:` status other than `1` (objects found) or `0` (none found), or
no status line at all. Thrown by [`crossmatch_catalog`](@ref)`(...;
:skybot)` once retries are exhausted, rather than reporting the
candidate as having no known counterpart: a silently empty match would
make every known object look like a new discovery.
"""
struct SkyBoTServiceError <: Exception
    message::String
end

Base.showerror(io::IO, e::SkyBoTServiceError) = print(io, "SkyBoTServiceError: ", e.message)

# Concurrent SkyBoT requests: it has no batch mode, and requests are mostly
# waiting on the network. Up to 60 ran clean; gains flatten past ~20.
const _SKYBOT_CONCURRENCY = 20

"""
    _crossmatch_skybot(candidates, radius; observatory=nothing) -> Table

Query SkyBoT once per candidate, concurrently (up to
`_SKYBOT_CONCURRENCY` requests in flight at a time via Julia
`Task`s, not extra threads — this is a network-latency-bound workload,
not a CPU-bound one, so cooperative concurrency on however many threads
Julia was started with is enough to see the full speedup). Order of
results does not depend on request completion order: each candidate's
matches are collected into their own slot and the final table is built
from those slots in `candidates`' original order, not arrival order.
"""
function _crossmatch_skybot(candidates, radius::Real;
                            observatory::Union{Nothing,AbstractString}=nothing)
    radius_deg = radius / 3600
    candidates_vec = collect(candidates)
    per_candidate = Vector{Vector{NamedTuple}}(undef, length(candidates_vec))

    semaphore = Base.Semaphore(_SKYBOT_CONCURRENCY)
    try
        @sync for (i, c) in enumerate(candidates_vec)
            @async begin
                Base.acquire(semaphore)
                try
                    per_candidate[i] = _skybot_matches(c, radius_deg; observatory)
                finally
                    Base.release(semaphore)
                end
            end
        end
    catch e
        # Callers recognise the request's own error, not @sync's wrapping of it.
        throw(_task_root(e))
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

"""
    _parse_skybot(text) -> Vector{NamedTuple}

Parse one SkyBoT text-mode cone-search response. Its first line is a
status, `# Flag: <n>`: `1` means a result table follows, `0` means no
object in the field (verified against the live service: a short "No
solar system object was found" message, still HTTP 200). Anything else
— `-1`, or no status line — is SkyBoT failing server-side, and throws a
[`SkyBoTServiceError`](@ref): before this check, a `-1` response's body
(a crash backtrace, no `|`-delimited rows) parsed as zero matches, so an
outage looked exactly like an empty field.
"""
function _parse_skybot(text::AbstractString)
    m = match(r"^#\s*Flag:\s*(-?\d+)", lstrip(text))
    flag = m === nothing ? nothing : parse(Int, m[1])
    if flag != 0 && flag != 1
        detail = first(strip(text), 200)
        throw(SkyBoTServiceError(flag === nothing ?
            "response has no '# Flag:' status line: $(repr(detail))" :
            "server reported Flag $flag: $(repr(detail))"))
    end

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
