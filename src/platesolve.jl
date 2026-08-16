const _ASTROMETRY_LOGIN_URL = "https://nova.astrometry.net/api/login"
const _ASTROMETRY_UPLOAD_URL = "https://nova.astrometry.net/api/upload"
const _ASTROMETRY_SUBMISSIONS_URL = "https://nova.astrometry.net/api/submissions"
const _ASTROMETRY_JOBS_URL = "https://nova.astrometry.net/api/jobs"
const _ASTROMETRY_WCS_URL = "https://nova.astrometry.net/wcs_file"

# Required by nova.astrometry.net's own API docs for programmatic file
# downloads, as an anti-scraper-bot check.
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

Validated end-to-end against the live service: a real ZTF frame (with
its existing WCS ignored) solved successfully, recovering sky coordinates
consistent with that frame's known field centre; a synthetic image with
no genuine star pattern correctly failed to solve rather than silently
returning a wrong answer. `test/runtests.jl`'s live round-trip test runs
whenever `ENV["ASTROMETRY_API_KEY"]` is set; the individual
request/response steps are always tested directly against fabricated
JSON, without a live connection.
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
