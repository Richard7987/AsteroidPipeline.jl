const _ASTROMETRY_LOGIN_URL = "https://nova.astrometry.net/api/login"
const _ASTROMETRY_UPLOAD_URL = "https://nova.astrometry.net/api/upload"
const _ASTROMETRY_SUBMISSIONS_URL = "https://nova.astrometry.net/api/submissions"
const _ASTROMETRY_JOBS_URL = "https://nova.astrometry.net/api/jobs"
const _ASTROMETRY_WCS_URL = "https://nova.astrometry.net/wcs_file"

# nova.astrometry.net requires it for programmatic file downloads.
const _ASTROMETRY_REFERER_HEADER = ["Referer" => _ASTROMETRY_LOGIN_URL]

"""
    plate_solve(image_path; api_key, timeout=300, poll_interval=5) -> WCSTransform

Solve for a FITS image's astrometric (WCS) solution via the
nova.astrometry.net web API, for frames with none already in the header
that [`load_wcs`](@ref) can parse — this requires an `api_key` rather
than any local plate-solving dependency.

Runs the full submit/poll/fetch cycle: login with `api_key` for a session
token, upload `image_path`, poll for a submission to produce a job, poll
that job until it succeeds or fails (both loops give up after `timeout`
seconds total, checking every `poll_interval` seconds), then download and
parse the resulting WCS with [`load_wcs`](@ref).

Errors (with the service's own message where given) on a login failure,
an upload failure, a job that fails to solve, or either poll loop timing
out — a frame astrometry.net cannot solve is a real failure, not
something to silently pass through.

The live round trip is tested when `ENV["ASTROMETRY_API_KEY"]` is set; the
individual steps are always tested against fabricated JSON.
"""
function plate_solve(image_path::AbstractString; api_key::AbstractString,
                      timeout::Real=300, poll_interval::Real=5)
    session = _astrometry_login(api_key)
    subid = _astrometry_upload(session, image_path)
    job_id = _astrometry_poll(
        () -> _astrometry_submission_job(subid),
        timeout, poll_interval,
        "astrometry.net submission $subid to produce a job",
    )
    _astrometry_poll(
        () -> _astrometry_job_done(job_id),
        timeout, poll_interval,
        "astrometry.net job $job_id to complete",
    )
    return _astrometry_fetch_wcs(job_id)
end

"""
    _astrometry_poll(check, timeout, poll_interval, description) -> Any

Call `check()` (which returns `nothing` if not ready yet, or the result
to return once it's ready) every `poll_interval` seconds until it returns
non-`nothing` or `timeout` seconds have elapsed, in which case this errors
mentioning `description`.
"""
function _astrometry_poll(check, timeout::Real, poll_interval::Real, description::AbstractString)
    deadline = time() + timeout
    while true
        result = check()
        result === nothing || return result
        time() >= deadline && error("timed out waiting for $description")
        sleep(poll_interval)
    end
end

_astrometry_login_body(api_key::AbstractString) =
    "request-json=" * HTTP.escapeuri(JSON.json(Dict("apikey" => api_key)))

function _parse_astrometry_response(text::AbstractString, action::AbstractString)
    d = JSON.parse(text)
    d["status"] == "success" ||
        error("astrometry.net $action failed: $(get(d, "errormessage", d["status"]))")
    return d
end

function _astrometry_login(api_key::AbstractString)
    response = HTTP.post(_ASTROMETRY_LOGIN_URL, ["Content-Type" => "application/x-www-form-urlencoded"],
                          _astrometry_login_body(api_key))
    return _parse_astrometry_response(String(response.body), "login")["session"]
end

function _astrometry_upload(session::AbstractString, image_path::AbstractString)
    request_json = JSON.json(Dict("session" => session, "publicly_visible" => "n"))
    return open(image_path) do io
        body = HTTP.Form(Dict(
            "request-json" => request_json,
            "file" => HTTP.Multipart(basename(image_path), io, "application/octet-stream"),
        ))
        response = HTTP.post(_ASTROMETRY_UPLOAD_URL, [], body)
        _parse_astrometry_response(String(response.body), "upload")["subid"]
    end
end

"""
    _astrometry_submission_job(subid) -> Union{Nothing,Integer}

The submission's first job id, once it has one; `nothing` if it doesn't
yet (still processing — the API can report a reserved-but-unassigned job
slot as a JSON `null`, i.e. `nothing` once parsed, as well as an empty
`jobs` array; both mean "not ready").
"""
function _astrometry_submission_job(subid)
    response = HTTP.get("$_ASTROMETRY_SUBMISSIONS_URL/$subid", _ASTROMETRY_REFERER_HEADER)
    d = JSON.parse(String(response.body))
    jobs = get(d, "jobs", [])
    isempty(jobs) && return nothing
    return jobs[1]
end

"""
    _astrometry_job_done(job_id) -> Union{Nothing,Bool}

`true` once `job_id` has solved successfully; errors immediately on a
reported failure (no point continuing to poll); `nothing` while still
processing.
"""
function _astrometry_job_done(job_id)
    response = HTTP.get("$_ASTROMETRY_JOBS_URL/$job_id", _ASTROMETRY_REFERER_HEADER)
    d = JSON.parse(String(response.body))
    status = d["status"]
    status == "success" && return true
    status == "failure" && error("astrometry.net job $job_id failed to solve")
    return nothing
end

function _astrometry_fetch_wcs(job_id)
    response = HTTP.get("$_ASTROMETRY_WCS_URL/$job_id", _ASTROMETRY_REFERER_HEADER)
    path, io = mktemp()
    write(io, response.body)
    close(io)
    header = FITS(path, "r") do f
        read_header(f[1], String)
    end
    rm(path)
    return load_wcs(header)
end

"""
    frame_pointing(header) -> Union{Nothing,Tuple{Float64,Float64}}

Where the telescope pointed, `(ra, dec)` in degrees, from the first pair of
the usual header keywords: `OBJCTRA`/`OBJCTDEC`, `RA`/`DEC`, `TELRA`/`TELDEC`,
`RA_OBJ`/`DEC_OBJ`, `OBJRA`/`OBJDEC`. A number is degrees; text is
sexagesimal, hours for RA. `nothing` if none is there.
"""
function frame_pointing(header)
    for (a, d) in (("OBJCTRA", "OBJCTDEC"), ("RA", "DEC"), ("TELRA", "TELDEC"),
                   ("RA_OBJ", "DEC_OBJ"), ("OBJRA", "OBJDEC"))
        haskey(header, a) && haskey(header, d) || continue
        ra, dec = _angle(header[a], 15.0), _angle(header[d], 1.0)
        ra === nothing || dec === nothing || return (mod(ra, 360.0), dec)
    end
    return nothing
end

function _angle(v, unit)
    v isa Real && return Float64(v)
    v isa AbstractString || return nothing
    s = strip(v)
    parts = split(s, r"[:\s]+"; keepempty=false)
    vals = tryparse.(Float64, parts)
    any(isnothing, vals) && return nothing
    length(vals) == 1 && return vals[1]
    x = sum(abs(p) / 60.0^(i - 1) for (i, p) in enumerate(vals))
    return unit * (startswith(s, "-") ? -x : x)
end

"""
    solve_wcs(image, ra, dec, stars; scale_range=(0.2, 3.0), search_radius=600.0,
              threshold=10.0, mask=nothing, n_bright=30, min_matches=6)
        -> (wcs, refined, n_matches, rms_arcsec, offset_px, n_detections, matches, scale)

A frame's astrometric solution when its header has none, found here rather
than by a web service: only the approximate pointing `(ra, dec)` (degrees,
within `search_radius` arcseconds; [`frame_pointing`](@ref) reads it from a
header) is needed. The pixel scale, anywhere in `scale_range` (″/px), the
orientation and a mirrored axis are found from the stars themselves.

`image` is in [`detect_sources`](@ref)'s `(y, x)` order, and `stars` has
`ra`, `dec` and `gmag` for every star the field could hold, e.g. from
[`gaia_reference_stars`](@ref) around the pointing. Triangles of
neighbouring stars are compared by the ratios of their sides, which no
shift, rotation, scale or mirror changes: the `n_bright` brightest
detections against as many catalogue stars per field area, for each scale
tried. Each matching pair of triangles gives a transformation, kept if
the most detections land on catalogue stars, and [`refine_wcs`](@ref)
then fits the solution to every star. `scale` is the solution's pixel
scale; `refined` is `false`, and `wcs` `nothing`, when no transformation
put `min_matches` detections on stars.
"""
function solve_wcs(image::AbstractMatrix{<:Real}, ra::Real, dec::Real, stars;
                   scale_range=(0.2, 3.0), search_radius::Real=600.0, threshold::Real=10.0,
                   mask::Union{Nothing,AbstractMatrix{Bool}}=nothing, n_bright::Integer=30,
                   min_matches::Integer=6)
    det = detect_sources(image; threshold, mask)
    order = sortperm(det.flux; rev=true)[1:min(n_bright, length(det))]
    D = [(det.xcen[i], det.ycen[i]) for i in order]
    failed = (wcs=nothing, refined=false, n_matches=0, rms_arcsec=NaN, offset_px=(0.0, 0.0),
              n_detections=length(det), matches=nothing, scale=NaN)
    length(D) >= 3 || return failed

    ny, nx = size(image)
    half_diagonal = hypot(nx, ny) / 2
    tangent(r, d) = _tangent_arcsec(r, d, ra, dec)
    S = [tangent(r, d) for (r, d) in zip(stars.ra, stars.dec)]
    by_mag = sortperm(collect(stars.gmag))
    image_triangles = _neighbour_triangles(D, 20.0)

    best = (n=0, s=0.0, a=0.0, b=0.0, t=(0.0, 0.0))
    lo, hi = scale_range
    s = lo
    while s <= hi * 1.2
        region = search_radius + half_diagonal * s
        inside = [i for i in by_mag if hypot(S[i]...) <= region]
        # as many stars per area as there are bright detections per frame area
        n = clamp(round(Int, length(D) * π * region^2 / (nx * ny * s^2)), length(D), 600)
        C = S[inside[1:min(n, length(inside))]]
        catalogue_triangles = sort!(_neighbour_triangles(C, 20.0 * s); by=first)
        firsts = first.(catalogue_triangles)
        for (u, v, idx) in image_triangles
            for k in searchsortedfirst(firsts, u - 0.01):searchsortedlast(firsts, u + 0.01)
                _, v2, cidx = catalogue_triangles[k]
                abs(v - v2) < 0.01 || continue
                fit = _similarity(D[idx], C[cidx])
                lo / 1.2 <= fit.s <= hi * 1.2 && fit.max_residual < 2.0 * fit.s || continue
                hits = count(D) do p
                    q = _apply_similarity(fit, p)
                    any(c -> hypot(q[1] - c[1], q[2] - c[2]) < 2.0 * fit.s, C)
                end
                hits > best.n && (best = (n=hits, fit=fit))
            end
        end
        s *= 1.2
    end
    best.n >= min_matches || return failed

    M = _similarity_matrix(best.fit)
    t = best.fit.t
    # the pixel the pointing falls on: M * p + t = 0
    crpix = -(M \ [t[1], t[2]])
    header = WCSTransform(2; crpix=crpix, crval=[Float64(ra), Float64(dec)], cdelt=[1.0, 1.0],
                          pc=M ./ 3600, ctype=["RA---TAN", "DEC--TAN"])
    r = refine_wcs(image, header, stars; threshold, mask, max_offset=10.0, min_matches)
    return (; r..., wcs=r.refined ? r.wcs : nothing, scale=best.fit.s)
end

# Gnomonic projection about (ra0, dec0), arcseconds: ξ toward east, η north.
function _tangent_arcsec(ra, dec, ra0, dec0)
    c = sind(dec0) * sind(dec) + cosd(dec0) * cosd(dec) * cosd(ra - ra0)
    ξ = cosd(dec) * sind(ra - ra0) / c
    η = (cosd(dec0) * sind(dec) - sind(dec0) * cosd(dec) * cosd(ra - ra0)) / c
    return (rad2deg(ξ) * 3600, rad2deg(η) * 3600)
end

# Each point with every pair of its 5 nearest neighbours, as
# (shortest / longest side, middle / longest, vertices ordered by opposite side).
function _neighbour_triangles(P, min_side)
    seen = Set{NTuple{3,Int}}()
    out = Tuple{Float64,Float64,Vector{Int}}[]
    dist(i, j) = hypot(P[i][1] - P[j][1], P[i][2] - P[j][2])
    for i in eachindex(P)
        nn = sort([j for j in eachindex(P) if j != i]; by=j -> dist(i, j))[1:min(5, length(P) - 1)]
        for a in 1:length(nn), b in a+1:length(nn)
            tri = (i, nn[a], nn[b])
            key = Tuple(sort([tri...]))
            key in seen && continue
            push!(seen, key)
            d = [dist(tri[2], tri[3]), dist(tri[1], tri[3]), dist(tri[1], tri[2])]
            p = sortperm(d)
            d[p[3]] >= min_side || continue
            push!(out, (d[p[1]] / d[p[3]], d[p[2]] / d[p[3]], [tri[k] for k in p]))
        end
    end
    return out
end

# Least-squares similarity P -> Q, with or without a mirrored x axis:
# q = s * R(θ) * (mirror ? (-x, y) : (x, y)) + t, as a = s cos θ, b = s sin θ.
function _similarity(P, Q)
    fits = map((false, true)) do mirror
        P2 = [(mirror ? -p[1] : p[1], p[2]) for p in P]
        mp = (sum(first, P2) / length(P2), sum(last, P2) / length(P2))
        mq = (sum(first, Q) / length(Q), sum(last, Q) / length(Q))
        num_a = num_b = den = 0.0
        for (p, q) in zip(P2, Q)
            x, y, u, v = p[1] - mp[1], p[2] - mp[2], q[1] - mq[1], q[2] - mq[2]
            num_a += x * u + y * v
            num_b += x * v - y * u
            den += x^2 + y^2
        end
        a, b = num_a / den, num_b / den
        t = (mq[1] - (a * mp[1] - b * mp[2]), mq[2] - (b * mp[1] + a * mp[2]))
        fit = (a=a, b=b, mirror=mirror, t=t, s=hypot(a, b))
        residual = maximum(hypot((_apply_similarity(fit, p) .- q)...) for (p, q) in zip(P, Q))
        (; fit..., max_residual=residual)
    end
    return argmin(f -> f.max_residual, fits)
end

_apply_similarity(f, p) = Tuple(_similarity_matrix(f) * [p[1], p[2]] .+ [f.t[1], f.t[2]])

_similarity_matrix(f) = [f.a -f.b; f.b f.a] * [(f.mirror ? -1.0 : 1.0) 0.0; 0.0 1.0]

# Scale ranges tried in turn, with the faintest Gaia magnitude each needs:
# the wider the field, the larger the sky to query and the brighter its limit.
const _SOLVE_LADDER = ((0.2, 0.6, 19.0), (0.6, 1.8, 18.0), (1.8, 3.0, 16.5))

# The fields solved so far: each one's pointing, Gaia stars, scale, and the
# frames solved there (epoch, pointing, WCS). Safe to share between threads:
# a WCS goes in and out as a copy (wcslib isn't safe on one shared between
# threads), and `query` lets one thread at a time fetch a new field's stars.
struct _SolveCache
    fields::Vector{Any}
    lock::ReentrantLock
    query::ReentrantLock
end
_SolveCache() = _SolveCache(Any[], ReentrantLock(), ReentrantLock())

_solved_field(ra, dec, stars, scale, epoch, wcs) =
    (ra=ra, dec=dec, stars=stars, scale=scale, solved=Any[(epoch=epoch, ra=ra, dec=dec, wcs=deepcopy(wcs))])

"""
    _header_or_solved_wcs(image, header_string, header, epoch, cache) -> (wcs, solved)

The frame's WCS from its header ([`load_wcs`](@ref)) or, when it has none,
from [`solve_wcs`](@ref) at its pointing ([`frame_pointing`](@ref)). Gaia
stars are fetched for the first frame of a field and kept in `cache` (a
`_SolveCache`) with the scale found, so the other frames of the field are
solved from them, or, where too few of their stars match, from the WCS of
the frame solved there nearest in time, refined ([`refine_wcs`](@ref)). Every field is kept: a
telescope whose header gives the pointing in two equinoxes, frame by frame,
queries Gaia twice, not at every change; and a frame of a known field that
none of this solves is not queried for again. With `reference_stars`, a
new field is solved from them before Gaia is queried. `solved` says which.
"""
function _header_or_solved_wcs(image, header_string, header, epoch, cache::_SolveCache; reference_stars=nothing)
    try
        return load_wcs(header_string), false
    catch e
        pointing = frame_pointing(header)
        pointing === nothing && rethrow(e)
        ra, dec = pointing
        known = lock(() -> length(cache.fields), cache.lock)
        w = _solve_from_cache(image, ra, dec, epoch, cache, 1:known)
        w === nothing || return w, true
        lock(cache.query)
        try
            # another thread may have fetched this field meanwhile
            now = lock(() -> length(cache.fields), cache.lock)
            w = _solve_from_cache(image, ra, dec, epoch, cache, known+1:now)
            w === nothing || return w, true
            if !any(c -> _angular_distance_arcsec(ra, dec, c.ra, c.dec) < 600, lock(() -> copy(cache.fields), cache.lock))
                if reference_stars !== nothing
                    r = solve_wcs(image, ra, dec, reference_stars)
                    if r.wcs !== nothing
                        field = _solved_field(ra, dec, reference_stars, r.scale, epoch, r.wcs)
                        lock(() -> push!(cache.fields, field), cache.lock)
                        return r.wcs, true
                    end
                end
                ny, nx = size(image)
                for (lo, hi, faintest) in _SOLVE_LADDER
                    radius = (600 + hypot(nx, ny) / 2 * hi) / 3600
                    stars = gaia_reference_stars(ra, dec, radius; epoch_jd=epoch, mag_range=(6.0, faintest))
                    r = solve_wcs(image, ra, dec, stars; scale_range=(lo, hi))
                    r.wcs === nothing && continue
                    lock(() -> push!(cache.fields, _solved_field(ra, dec, stars, r.scale, epoch, r.wcs)), cache.lock)
                    return r.wcs, true
                end
            end
        finally
            unlock(cache.query)
        end
        error("no WCS in the header, and none found around its pointing ($ra, $dec)")
    end
end

# The frame solved from the cached fields `range` within 10′ of its pointing,
# nearest first; `nothing` if none solves it.
function _solve_from_cache(image, ra, dec, epoch, cache::_SolveCache, range)
    distance(c) = _angular_distance_arcsec(ra, dec, c.ra, c.dec)
    fields = lock(() -> cache.fields[range], cache.lock)
    for c in sort(filter(c -> distance(c) < 600, fields); by=distance)
        r = solve_wcs(image, ra, dec, c.stars; scale_range=(0.95, 1.05) .* c.scale)
        if r.wcs === nothing
            # a sparse field can leave too few triangles; the WCS of the frame
            # nearest in time, moved by the change in pointing, refines instead
            near = lock(cache.lock) do
                p = argmin(q -> abs(q.epoch - epoch), c.solved)
                (ra=p.ra, dec=p.dec, wcs=deepcopy(p.wcs))
            end
            w = near.wcs
            w.crval = w.crval .+ [ra - near.ra, dec - near.dec]
            r = refine_wcs(image, w, c.stars)
            r = r.refined ? r : (wcs=nothing,)
        end
        r.wcs === nothing && continue
        lock(() -> push!(c.solved, (epoch=epoch, ra=ra, dec=dec, wcs=deepcopy(r.wcs))), cache.lock)
        return r.wcs
    end
    return nothing
end

# A copy of the frame in `dir` with the WCS as header keywords: the original
# is never written to, even through a symbolic link. A frame already in `dir`
# (a cleaned copy) is written in place.
function _copy_with_wcs(path, w::WCSTransform, dir)
    out = joinpath(dir, basename(path))
    if abspath(out) != abspath(path)
        mkpath(dir)
        cp(path, out; force=true, follow_symlinks=true)
        chmod(out, 0o644)
    end
    # WCS.jl's `pc` field reads transposed from the FITS PCi_j it was built with
    cd = permutedims(w.pc) .* w.cdelt
    FITS(out, "r+") do f
        for (k, v) in (("CTYPE1", w.ctype[1]), ("CTYPE2", w.ctype[2]), ("CRVAL1", w.crval[1]),
                       ("CRVAL2", w.crval[2]), ("CRPIX1", w.crpix[1]), ("CRPIX2", w.crpix[2]),
                       ("CD1_1", cd[1, 1]), ("CD1_2", cd[1, 2]), ("CD2_1", cd[2, 1]), ("CD2_2", cd[2, 2]))
            write_key(f[1], k, v)
        end
    end
    return out
end
