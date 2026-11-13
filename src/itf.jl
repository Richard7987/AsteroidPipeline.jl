# The MPC's Isolated Tracklet File: observations reported to the MPC and not
# yet linked to any object, Pan-STARRS's own unlinked detections among them.

const ITF_URL = "https://www.minorplanetcenter.net/iau/ITF/itf.txt.gz"

"""
    default_itf_path() -> String

Where [`download_itf`](@ref) keeps its copy, next to [`default_mpcorb_path`](@ref):
`\$XDG_CACHE_HOME` (default `~/.cache`) `/AsteroidPipeline/itf.txt.gz`.
"""
default_itf_path() = joinpath(dirname(default_mpcorb_path()), "itf.txt.gz")

"""
    download_itf(path=default_itf_path(); max_age_hours=24, url=ITF_URL) -> String

A current copy of the MPC's Isolated Tracklet File at `path`, kept gzipped
(~130 MB, ten million observations): downloaded if missing or older than
`max_age_hours`. If the download fails and an older copy exists, that copy is
returned with a warning; with no copy at all, the error propagates.
"""
function download_itf(path::AbstractString=default_itf_path(); max_age_hours::Real=24,
                      url::AbstractString=ITF_URL)
    isfile(path) && (time() - mtime(path)) < max_age_hours * 3600 && return path
    part = path * ".part"
    try
        mkpath(dirname(path))
        write(part, HTTP.get(url; read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body)
        # a truncated or HTML error body must not replace a good copy
        first_lines = open(io -> collect(Iterators.take(eachline(GzipDecompressorStream(io)), 100)), part)
        count(l -> _parse_mpc80_line(l) !== nothing, first_lines) >= 90 ||
            error("downloaded ITF does not hold 80-column observations")
        mv(part, path; force=true)
    catch e
        rm(part; force=true)
        isfile(path) || rethrow()
        @warn "could not refresh the ITF; using the existing copy" path exception=e
    end
    return path
end

"""
    read_itf(path, epochs; margin_days=1) -> Vector{NamedTuple}

The observations in the Isolated Tracklet File at `path` (gzipped, as
[`download_itf`](@ref) keeps it) taken within `margin_days` of `epochs` (JD
UTC): `trksub`, `jd`, `ra`, `dec` (degrees) and `station`. Lines of other
dates are skipped unparsed, so a night's observations are read in seconds.
"""
function read_itf(path::AbstractString, epochs; margin_days::Real=1)
    days = floor(Int, minimum(epochs) - 0.5 - margin_days):ceil(Int, maximum(epochs) - 0.5 + margin_days)
    dates = Set(Dates.format(Dates.julian2datetime(d + 0.5), "yyyy mm dd") for d in days)
    obs = NamedTuple[]
    open(path) do io
        for line in eachline(GzipDecompressorStream(io))
            length(line) >= 25 && line[16:25] in dates || continue
            o = _parse_mpc80_line(line)
            o === nothing && continue
            push!(obs, (trksub=String(o.provisional), jd=o.jd, ra=o.ra, dec=o.dec, station=String(o.station)))
        end
    end
    return obs
end

"""
    itf_check(candidates, itf; radius_arcsec=2.0, max_gap_hours=0.5) -> Dict{Int,Vector}

The observations of `itf` ([`read_itf`](@ref)'s) that each tracklet of
`candidates` (`id`, `ra`, `dec`, `epoch`) passes through: within
`radius_arcsec` of the tracklet's straight-line motion, at most
`max_gap_hours` before its first point or after its last. A tracklet with
such observations was most likely reported already, unlinked: by
Pan-STARRS's own pipeline, for an IASC set, which IASC rejects as a
duplicate. Only tracklets with matches are keys; each match is the ITF
observation with its `distance_arcsec` from the motion.
"""
function itf_check(candidates, itf; radius_arcsec::Real=2.0, max_gap_hours::Real=0.5)
    checks = Dict{Int,Vector{NamedTuple}}()
    isempty(itf) && return checks
    rows = collect(candidates)
    for id in unique(r.id for r in rows)
        pts = filter(r -> r.id == id, rows)
        length(pts) >= 2 || continue
        ra0, dec0, t0 = pts[1].ra, pts[1].dec, pts[1].epoch
        offset(ra, dec) = (3600 * (mod(ra - ra0 + 180, 360) - 180) * cosd(dec0), 3600 * (dec - dec0))
        t = [p.epoch - t0 for p in pts]
        xy = [offset(p.ra, p.dec) for p in pts]
        fit(v) = (A = hcat(ones(length(t)), t); A \ v)
        cx, cy = fit(first.(xy)), fit(last.(xy))
        lo, hi = minimum(t) - max_gap_hours / 24, maximum(t) + max_gap_hours / 24
        matches = NamedTuple[]
        for o in itf
            dt = o.jd - t0
            lo <= dt <= hi || continue
            x, y = offset(o.ra, o.dec)
            d = hypot(x - (cx[1] + cx[2] * dt), y - (cy[1] + cy[2] * dt))
            d <= radius_arcsec && push!(matches, merge(o, (distance_arcsec=d,)))
        end
        isempty(matches) || (checks[id] = matches)
    end
    return checks
end
