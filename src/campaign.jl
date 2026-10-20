#=
One image set end to end — a same-field sequence of frames, as a .zip, a
folder or a list of FITS files — to the MPC report that is submitted for it.
Built and validated on IASC (Pan-STARRS1) campaign sets; everything here was
written after real sets from a live campaign (2026-10) went wrong in ways the
detection code alone didn't cover: a report with no `COD` line, a known
asteroid measured by mistake into a "no moving objects" report, another
set's lines left in Astrometrica's report file, a set with one frame Astrometrica
couldn't reduce, and an orbit SkyBoT didn't know.
=#

"""
    pipeline_config(pixel_scale_arcsec) -> NamedTuple

[`run_pipeline`](@ref) keywords for a same-field sequence of frames at this
pixel scale (arcsec/pixel): the configuration validated against every local
IASC set (see the IASC validation page of the docs), with its tolerances kept
in arcseconds and converted to this scale's pixels — a 2" linking radius,
speeds of 3"/h to ~210"/h. Detection on each frame minus the sequence's
median at 6σ, fill-value gaps masked, single-pixel spikes, hot pixels and
bright-star light dropped, each frame's WCS refit to Gaia DR3, tracklets held
to straight-line, constant-brightness tests, and 3 frames enough for a
tracklet.
"""
pipeline_config(pixel_scale_arcsec::Real) = (
    threshold=6.0, match_radius=2.0 / pixel_scale_arcsec, max_speed=5050.0 / pixel_scale_arcsec,
    mask_fill_values=true, refine_astrometry=true, difference_stack=true, min_sharpness=0.3,
    min_speed=3.0 * 24 / pixel_scale_arcsec, max_flux_ratio=2.5, max_residual=1.0, min_frames=3,
    max_background_excess=1.0, detector_artifact_radius=2.0)

"""
    IASC_PS1_CONFIG

[`pipeline_config`](@ref) for Pan-STARRS1's IASC frames (0.2563"/pixel).
"""
const IASC_PS1_CONFIG = pipeline_config(0.2563)

# Arcsec per pixel at the image's centre, through `wcs`.
function _pixel_scale_arcsec(wcs, image_size)
    ny, nx = image_size
    a = pix_to_world(wcs, Float64[nx / 2, ny / 2])
    b = pix_to_world(wcs, Float64[nx / 2 + 1, ny / 2])
    return _angular_distance_arcsec(a..., b...)
end

const _FITS_EXTENSIONS = (".fits", ".fit", ".fts")

"""
    extract_set(zip_path, dest_dir) -> (name, paths)

Unpack an image set (an IASC one: `ps1-20261005_9_XY45_p01.zip`) into
`dest_dir/<set name>/` and return the set's name (the folder inside the
archive, `XY45_p01`, or the archive's own name if it has none) and its FITS
paths, sorted. An already-extracted set is reused, not overwritten.

The archive is unpacked into a scratch directory first, and only its
regular FITS files are moved out, by base name — so a crafted entry name or
a symbolic link inside it can't place anything else in `dest_dir`.
"""
function extract_set(zip_path::AbstractString, dest_dir::AbstractString)
    isfile(zip_path) || throw(ArgumentError("no such archive: $zip_path"))
    mkpath(dest_dir)
    scratch = mktempdir(dest_dir; prefix=".extract-")
    try
        run(pipeline(`$(p7zip_jll.p7zip()) x -y -bd -o$scratch $zip_path`; stdout=devnull, stderr=devnull))
        fits = String[]
        for (root, _, files) in walkdir(scratch), f in files
            p = joinpath(root, f)
            (!islink(p) && isfile(p) && any(e -> endswith(lowercase(f), e), _FITS_EXTENSIONS)) && push!(fits, p)
        end
        isempty(fits) && throw(ArgumentError("no FITS files in $zip_path"))
        parent = dirname(first(fits))
        name = parent == scratch ? first(splitext(basename(zip_path))) : basename(parent)
        target = joinpath(dest_dir, name)
        existing = isdir(target) ? _fits_in(target) : String[]
        if isempty(existing)
            mkpath(target)
            for f in fits
                mv(f, joinpath(target, basename(f)); force=true)
            end
        end
        return (name=name, paths=_fits_in(target))
    finally
        rm(scratch; recursive=true, force=true)
    end
end

_fits_in(dir) = sort([joinpath(dir, f) for f in readdir(dir)
                      if any(e -> endswith(lowercase(f), e), _FITS_EXTENSIONS) && isfile(joinpath(dir, f))])

"""
    clean_blank_pixels(path; outdir) -> String

Write a copy of the FITS file at `path` with `BLANK`-sentinel pixels replaced
by the frame's own valid-region median, into `outdir`, header untouched (WCS
and epochs stay exact), and return the copy's path — or `path` itself if the
frame has no `BLANK` pixels. An existing copy is reused.

PS1's IASC frames mark invalid pixels with the `BLANK` keyword (scaled
through `BZERO`/`BSCALE`), which FITSIO.jl doesn't turn into `NaN`: left in,
one real frame read 158,443 sentinel pixels as bright flux and gave 176,165
spurious detections at 8σ (176 once cleaned).
"""
function clean_blank_pixels(path::AbstractString; outdir::AbstractString)
    out = joinpath(outdir, basename(path))
    isfile(out) && return out
    cleaned = FITS(path, "r") do fin
        hdu = fin[1]
        h = read_header(hdu)
        haskey(h, "BLANK") || return false
        bzero = haskey(h, "BZERO") ? h["BZERO"] : 0.0
        bscale = haskey(h, "BSCALE") ? h["BSCALE"] : 1.0
        sentinel = h["BLANK"] * bscale + bzero
        img = Float64.(read(hdu))
        is_blank = img .== sentinel
        any(is_blank) || return false
        img[is_blank] .= median(img[.!is_blank])
        mkpath(outdir)
        FITS(out, "w") do fout
            write(fout, img; header=h)
        end
        return true
    end
    return cleaned ? out : path
end

const MPCORB_URL = "https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz"

"""
    default_mpcorb_path() -> String

Where [`download_mpcorb`](@ref) keeps its copy: `\$XDG_CACHE_HOME`
(default `~/.cache`) `/AsteroidPipeline/MPCORB.DAT`.
"""
default_mpcorb_path() = joinpath(get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache")),
                                 "AsteroidPipeline", "MPCORB.DAT")

"""
    download_mpcorb(path=default_mpcorb_path(); max_age_hours=24, url=MPCORB_URL) -> String

A current copy of the MPC's orbit file at `path`: downloaded (gzipped, ~70 MB)
and unpacked if the file is missing or older than `max_age_hours`, kept as
is otherwise. If the download fails and an older copy exists, that copy is
returned with a warning; with no copy at all, the error propagates.

Why: on a 2026-10-04 set, SkyBoT placed 2019 UA185 far from where it was,
while the MPC's orbit file of that same day put it within ~4" — the file is
only as useful as it is recent. Astrometrica's own copy is refreshed only
when someone updates it from its menu.
"""
function download_mpcorb(path::AbstractString=default_mpcorb_path(); max_age_hours::Real=24,
                         url::AbstractString=MPCORB_URL)
    isfile(path) && (time() - mtime(path)) < max_age_hours * 3600 && return path
    part = path * ".part"
    try
        mkpath(dirname(path))
        body = HTTP.get(url; read_idle_timeout=_CATALOG_IDLE_TIMEOUT).body
        open(part, "w") do io
            write(io, GzipDecompressorStream(IOBuffer(body)))
        end
        # a truncated or HTML error body must not replace a good copy
        n = countlines(part)
        n > 100_000 || error("downloaded MPCORB.DAT has only $n lines")
        mv(part, path; force=true)
    catch e
        rm(part; force=true)
        isfile(path) || rethrow()
        @warn "could not refresh MPCORB.DAT; using the existing copy" path exception=e
    end
    return path
end

"""
    astrometrica_settings(cfg_path) -> NamedTuple

The observing-site fields of an Astrometrica settings file (`Settings/PS1.cfg`
in IASC's installer): `code` (MPC observatory code), `observers`,
`measurers`, `telescope`, `contact`, and `net` — the reference catalog's name
as Astrometrica writes it in a report's `NET` line, for the catalog codes
seen in real reports (`10`, "Gaia DR1"), else `""`.
"""
function astrometrica_settings(cfg_path::AbstractString)
    kv = Dict{String,String}()
    for line in eachline(cfg_path)
        line = strip(rstrip(line, '\r'))
        (isempty(line) || startswith(line, '[')) && continue
        k, sep, v = partition_first(line, '=')
        sep && (kv[k] = v)
    end
    nets = Dict("10" => "Gaia DR1")
    return (code=get(kv, "MPCCode", ""), observers=get(kv, "Observer", ""), measurers=get(kv, "Measurer", ""),
            telescope=get(kv, "Telescope", ""), contact=get(kv, "Contact1", ""),
            net=get(nets, get(kv, "Catalog", ""), ""))
end

function partition_first(s::AbstractString, c::Char)
    i = findfirst(c, s)
    i === nothing && return (String(s), false, "")
    return (String(strip(s[1:prevind(s, i)])), true, String(strip(s[nextind(s, i):end])))
end

"""
    mpc_report_header(code; observers="", measurers="", telescope="", net="",
                      contact="", ack_time=Dates.now()) -> Vector{String}

The header lines of an MPC 80-column report, in the order Astrometrica writes
them: `COD` first (IASC's portal rejects a report that doesn't start with
it), then `CON`, `OBS`, `MEA`, `TEL`, `ACK`, `NET` — each only if given.
"""
function mpc_report_header(code::AbstractString; observers::AbstractString="", measurers::AbstractString="",
                           telescope::AbstractString="", net::AbstractString="", contact::AbstractString="",
                           ack_time::DateTime=Dates.now())
    length(code) == 3 || throw(ArgumentError("code must be a 3-character MPC observatory code"))
    lines = ["COD $code"]
    for (key, value) in (("CON", contact), ("OBS", observers), ("MEA", measurers), ("TEL", telescope))
        isempty(value) || push!(lines, "$key $value")
    end
    push!(lines, "ACK MPCReport file updated " * Dates.format(ack_time, "yyyy.mm.dd HH:MM:SS"))
    isempty(net) || push!(lines, "NET $net")
    return lines
end

const _MPC_REPORT_END = "----- end -----"

"""
    blank_mpc_report(code; header...) -> String

The report for a set with no moving objects: header and end line, no
observations — what Astrometrica writes after a reduction once the image has
been clicked, and what IASC's portal expects for "no moving objects
detected". Keywords as for [`mpc_report_header`](@ref).
"""
blank_mpc_report(code::AbstractString; header...) =
    join([mpc_report_header(code; header...); _MPC_REPORT_END], "\n") * "\n"

"""
    mpc_report(candidates, code; trksub_prefix="", designations=nothing, header...) -> String

[`mpc80_report`](@ref)'s observation lines for `candidates`, between a
[`mpc_report_header`](@ref) and the end line: a whole report, in the form
Astrometrica produces.
"""
function mpc_report(candidates, code::AbstractString; trksub_prefix::AbstractString="", designations=nothing,
                    header...)
    obs = isempty(candidates) ? String[] : split(chomp(mpc80_report(candidates, code; trksub_prefix, designations)), '\n')
    return join([mpc_report_header(code; header...); obs; _MPC_REPORT_END], "\n") * "\n"
end

const _BASE62 = ['0':'9'; 'A':'Z'; 'a':'z']

"""
    pack_designation(name) -> Union{String,Nothing}

The MPC packed form of a minor-planet designation, as columns 1-12 of an
80-column report carry it: a number (`"(656932) 2016 FA12"` or `"656932"`
-> `"~09bg"`, `"433"` -> `"00433"`), else a provisional designation
(`"2016 SB115"` -> `"K16SB5B"`) or a survey one (`"4311 T-1"` ->
`"T1S4311"`). `nothing` for a name it can't pack (a bare proper name).
"""
function pack_designation(name::AbstractString)
    s = strip(name)
    m = match(r"^\((\d+)\)", s)
    m === nothing && (m = match(r"^(\d+)$", s))
    m === nothing || return _pack_number(parse(Int, m[1]))
    m = match(r"^(\d{4}) ([A-HJ-Y])([A-HJ-Z])(\d*)$", s)
    if m !== nothing
        year = parse(Int, m[1]); cycle = isempty(m[4]) ? 0 : parse(Int, m[4])
        century = get(Dict(18 => 'I', 19 => 'J', 20 => 'K'), year ÷ 100, nothing)
        (century === nothing || cycle >= 620) && return nothing
        c = cycle < 100 ? lpad(cycle, 2, '0') : string(_BASE62[cycle ÷ 10 + 1], cycle % 10)
        return string(century, lpad(year % 100, 2, '0'), m[2], c, m[3])
    end
    m = match(r"^(\d{4}) (P-L|T-1|T-2|T-3)$", s)
    m === nothing || return string(Dict("P-L" => "PL", "T-1" => "T1", "T-2" => "T2", "T-3" => "T3")[m[2]], 'S', m[1])
    return nothing
end

function _pack_number(n::Integer)
    n < 100_000 && return lpad(n, 5, '0')
    n < 620_000 && return string(_BASE62[n ÷ 10_000 + 1], lpad(n % 10_000, 4, '0'))
    k = n - 620_000
    k < 62^4 || return nothing
    return "~" * join(_BASE62[(k ÷ 62^p) % 62 + 1] for p in 3:-1:0)
end

# One 80-column observation line -> designation fields, epoch (JD UTC), RA/Dec
# (degrees) and station; `nothing` if it isn't one.
function _parse_mpc80_line(line::AbstractString)
    line = rstrip(line, ['\r', '\n'])
    (length(line) == 80 && isascii(line)) || return nothing
    m = match(r"^(\d{4}) (\d\d) (\d\d\.\d+)\s*$", line[16:32])
    m === nothing && return nothing
    t = match(r"^(\d\d) (\d\d) (\d\d\.?\d*)\s*$", line[33:44])
    d = match(r"^([+-])(\d\d) (\d\d) (\d\d\.?\d*)\s*$", line[45:56])
    (t === nothing || d === nothing) && return nothing
    jd = datetime2julian(DateTime(parse(Int, m[1]), parse(Int, m[2]), 1)) + parse(Float64, m[3]) - 1
    ra = 15 * (parse(Float64, t[1]) + parse(Float64, t[2]) / 60 + parse(Float64, t[3]) / 3600)
    dec = (d[1] == "-" ? -1 : 1) * (parse(Float64, d[2]) + parse(Float64, d[3]) / 60 + parse(Float64, d[4]) / 3600)
    return (number=strip(line[1:5]), provisional=strip(line[6:12]), jd=jd, ra=ra, dec=dec, station=line[78:80])
end

const _MPC_HEADER_KEYS = ("COD", "CON", "OBS", "MEA", "TEL", "NET", "ACK", "AC2", "COM", "NUM", "BND")

"""
    check_mpc_report(text; observatory, epochs=nothing, field=nothing, known=nothing,
                     epoch_tolerance=120.0, position_tolerance=10.0)
        -> (ok, issues, cleaned)

Check an MPC 80-column report — normally Astrometrica's `MPCReport.txt` —
before it's submitted, for the mistakes made on a real campaign, and return a
corrected version.

`issues` is a vector of `(code, line, detail)` named tuples (`line` is
the 1-based line number, `0` for the report as a whole); `ok` is `false` if
any of them is an error. Codes:

- `:empty` — nothing in it. *Error.*
- `:missing_cod` — no `COD` line: IASC's portal answers "No Telescope Code
  found at start of MPC Report", as it did for a report copied before the
  image was clicked once in Astrometrica. *Error.*
- `:cod_not_first` — a `COD` line, but not first. *Error.*
- `:wrong_code` — `COD` isn't `observatory` (`detail`: the code found; `XXX`
  is Astrometrica's default settings, not the PS1 ones). *Error.*
- `:station_mismatch` — an observation's columns 78-80 aren't
  `observatory`. *Error.*
- `:malformed_line` — neither a header line, an 80-column observation, nor
  the end line. *Error.*
- `:other_set` — an observation more than `epoch_tolerance` seconds from
  every one of `epochs` (the set's frame epochs, JD UTC mid-exposure), or
  off the set's frames (`field(ra, dec)` false): another set's line still in
  Astrometrica's report file. Both tests are needed — IASC's sets from one
  night come from the same exposures (`XY45_p01` and `XY45_p10` share their
  epochs to the second), so only the position tells them apart. *Error.*
- `:known_object` — an observation of a known object (`detail`: its name):
  its packed designation is one of `known`'s, or it lies within
  `position_tolerance` arcsec of one's position at that epoch. Clicking on a
  known asteroid by mistake put such a line into a "no moving objects"
  report. *Warning* — a known object may be reported on purpose.
- `:single_observation` — an object with one observation only; the MPC
  needs at least two. *Warning.*
- `:missing_end` — no `----- end -----` line. *Warning.*

`known` is a vector of named tuples with `name` and optionally `packed`, and
`epochs`/`ra`/`dec` vectors (as [`known_objects_in_field`](@ref) returns).
`cleaned` is the report with every observation flagged `:other_set` or
`:known_object` and every malformed line removed, `COD` (`observatory`)
first and the end line last — a blank report if nothing is left.
"""
function check_mpc_report(text::AbstractString; observatory::AbstractString, epochs=nothing,
                          field=nothing, known=nothing, epoch_tolerance::Real=120.0, position_tolerance::Real=10.0)
    issues = NamedTuple{(:code, :line, :detail),Tuple{Symbol,Int,String}}[]
    flag(code, line=0, detail="") = push!(issues, (code=code, line=line, detail=String(detail)))
    lines = [rstrip(l, ['\r']) for l in split(text, '\n')]
    while !isempty(lines) && isempty(strip(lines[end]))
        pop!(lines)
    end
    if all(l -> isempty(strip(l)), lines)
        flag(:empty)
        return (ok=false, issues=issues, cleaned=blank_mpc_report(observatory))
    end

    header = String[]; kept = String[]; obs_by_object = Dict{String,Vector{Int}}()
    cod_line = 0; has_end = false; first_content = 0
    known_keys = Dict{String,String}()
    if known !== nothing
        for k in known
            for key in (hasproperty(k, :packed) ? k.packed : "", something(pack_designation(k.name), ""))
                isempty(key) || (known_keys[key] = k.name)
            end
        end
    end
    for (i, raw) in enumerate(lines)
        l = String(raw)
        isempty(strip(l)) && continue
        first_content == 0 && (first_content = i)
        if strip(l) == _MPC_REPORT_END
            has_end = true
            continue
        end
        key = length(l) >= 3 ? l[1:3] : ""
        if key in _MPC_HEADER_KEYS && (length(l) == 3 || l[4] == ' ')
            if key == "COD"
                cod_line = i
                code = strip(l[4:end])
                code == observatory || flag(:wrong_code, i, code)
                i == first_content || flag(:cod_not_first, i)
            else
                push!(header, l)
            end
            continue
        end
        obs = _parse_mpc80_line(l)
        if obs === nothing
            flag(:malformed_line, i, l)
            continue
        end
        obs.station == observatory || flag(:station_mismatch, i, obs.station)
        if (epochs !== nothing && !isempty(epochs) &&
            minimum(abs(obs.jd - e) for e in epochs) * 86400 > epoch_tolerance) ||
           (field !== nothing && !field(obs.ra, obs.dec))
            flag(:other_set, i, isempty(obs.number) ? obs.provisional : obs.number)
            continue
        end
        name = get(known_keys, obs.number, get(known_keys, obs.provisional, nothing))
        if name === nothing && known !== nothing
            for k in known
                hasproperty(k, :epochs) || continue
                j = findfirst(e -> abs(e - obs.jd) * 86400 <= epoch_tolerance, k.epochs)
                j === nothing && continue
                _angular_distance_arcsec(obs.ra, obs.dec, k.ra[j], k.dec[j]) <= position_tolerance && (name = k.name; break)
            end
        end
        if name !== nothing
            flag(:known_object, i, name)
            continue
        end
        push!(kept, l)
        push!(get!(obs_by_object, obs.number * "|" * obs.provisional, Int[]), i)
    end
    cod_line == 0 && flag(:missing_cod)
    for (object, idx) in obs_by_object
        length(idx) == 1 && flag(:single_observation, idx[1], strip(replace(object, '|' => ' ')))
    end
    has_end || flag(:missing_end)
    sort!(issues; by=x -> x.line)

    errors = (:empty, :missing_cod, :cod_not_first, :wrong_code, :station_mismatch, :malformed_line, :other_set)
    cleaned = join(["COD $observatory"; header; kept; _MPC_REPORT_END], "\n") * "\n"
    return (ok=!any(x -> x.code in errors, issues), issues=issues, cleaned=cleaned)
end

# "(656932) 2016 FA12" and "2016 FA12", "(1234) Hasslacher" and "Hasslacher":
# the same object as MPCORB and SkyBoT name it.
_bare_name(name::AbstractString) = String(strip(replace(name, r"^\(\d+\)\s*" => "")))

"""
    known_objects_in_field(wcs_per_frame, epochs, image_size; mpcorb=nothing,
                           observatory, skybot=true) -> (objects, warnings)

Every catalogued minor planet whose predicted position falls on the frames at
any of `epochs` (JD UTC, one per frame), whether or not the pipeline detected
it: what Astrometrica's known-object overlay shows, for checking a set by eye
and for not measuring a known asteroid by mistake.

`image_size` is `(ny, nx)`. Positions come from `mpcorb` (what
[`read_mpcorb`](@ref) returns, propagated as in [`crossmatch_catalog`](@ref)'s
`:mpcorb`; observer positions from JPL Horizons) and, with `skybot`, from a
SkyBoT cone search at each epoch; an object both know keeps MPCORB's
position, the more current one. Each object is a named tuple: `name`,
`packed` (MPCORB's packed designation, `""` if SkyBoT only), `mag`
(predicted V), `source`, and per-frame vectors `epochs`, `ra`, `dec`, `x`,
`y` (pixels of each frame's own `wcs`). A network failure leaves that
source out and adds a line to `warnings` instead of failing.
"""
function known_objects_in_field(wcs_per_frame::AbstractVector, epochs::AbstractVector{<:Real}, image_size;
                                mpcorb=nothing, observatory::AbstractString, skybot::Bool=true)
    ny, nx = image_size
    w1 = wcs_per_frame[1]
    centre = pix_to_world(w1, Float64[nx / 2, ny / 2])
    radius = maximum(_angular_distance_arcsec(centre..., pix_to_world(w1, Float64[x, y])...)
                     for (x, y) in ((1, 1), (nx, 1), (1, ny), (nx, ny))) / 3600
    onframe(x, y) = 1 <= x <= nx && 1 <= y <= ny
    objects = NamedTuple[]
    warnings = String[]
    function add!(name, packed, mag, source, ras, decs)
        px = [world_to_pix(w, Float64[r, d]) for (w, r, d) in zip(wcs_per_frame, ras, decs)]
        any(p -> onframe(p...), px) || return
        push!(objects, (name=String(name), packed=String(packed), mag=Float64(mag), source=source,
                        epochs=collect(Float64, epochs), ra=collect(Float64, ras), dec=collect(Float64, decs),
                        x=first.(px), y=last.(px)))
    end

    if mpcorb !== nothing
        try
            orbits = mpcorb isa AbstractString ? read_mpcorb(mpcorb) : mpcorb
            observers = _observer_heliocentric(collect(Float64, epochs), observatory)
            for o in orbits
                p = _mpcorb_radec(o, epochs[1], observers[1])
                _angular_distance_arcsec(p.ra, p.dec, centre...) <= (radius + 0.1) * 3600 || continue
                ps = [_mpcorb_radec(o, e, obs) for (e, obs) in zip(epochs, observers)]
                mag = _hg_magnitude(o.H, o.G, ps[1].r, ps[1].Δ, hypot(observers[1]...))
                add!(o.name, o.packed, mag, "MPCORB", [q.ra for q in ps], [q.dec for q in ps])
            end
        catch e
            _is_service_failure(e) || e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            push!(warnings, "MPCORB positions unavailable (JPL Horizons): $(sprint(showerror, e))")
        end
    end

    if skybot
        try
            per_epoch = [Dict(m.name => m for m in _skybot_matches((ra=centre[1], dec=centre[2], epoch=e),
                                                                  radius + 0.02; observatory))
                         for e in epochs]
            have = Set(_bare_name(o.name) for o in objects)
            for name in keys(per_epoch[1])
                _bare_name(name) in have && continue
                all(d -> haskey(d, name), per_epoch) || continue
                add!(name, "", per_epoch[1][name].mv, "SkyBoT", [d[name].ra for d in per_epoch],
                     [d[name].dec for d in per_epoch])
            end
        catch e
            (e isa HTTP.Exceptions.HTTPError || e isa SkyBoTServiceError || e isa Base.IOError) || rethrow()
            push!(warnings, "SkyBoT unavailable: $(sprint(showerror, e))")
        end
    end
    sort!(objects; by=o -> o.mag)
    return (objects=objects, warnings=warnings)
end

"""
    analyze_set(source; workdir, observatory, mpcorb=nothing, header=(;), name=nothing,
                extract_to=dirname(source), config=nothing,
                progress=(stage, fraction) -> nothing) -> NamedTuple

Everything [`run_pipeline`](@ref) and the checks around it do for one image
set — a same-field sequence of at least 3 frames — `source` being its `.zip`,
its directory of FITS frames, or a vector of FITS paths. `observatory` is the
MPC code of the site that took them (`"F51"` for Pan-STARRS1).

1. unpack a `.zip` ([`extract_set`](@ref)) into `extract_to`, next to it by
   default, where Astrometrica can open it too; frames are taken in time
   order, whatever order they were given in;
2. clean `BLANK` pixels ([`clean_blank_pixels`](@ref), copies under
   `workdir`);
3. refine every frame's WCS against Gaia DR3, the frames together
   ([`refine_wcs_set`](@ref)), and
   measure its quality — sources at 10σ per reference star, as
   [`set_quality`](@ref) does. A frame above the limit (30) is left out if at
   least 3 others remain — on a real 2026-10-04 set, one frame Astrometrica
   couldn't reduce made the whole set look unusable; with fewer than 3 good
   frames the set is `reliable = false` and not searched;
4. search the remaining frames with `config` — by default
[`pipeline_config`](@ref) at the frames' own pixel scale;
5. identify tracklets with SkyBoT (5") and, given `mpcorb`, the MPC's orbits
   (15"), confirmed in at least 3 frames (all, with 3 frames). `mpcorb` is
   used only for images within a year of its orbits' epoch (it is renewed
   about twice a year; 130 days out, it placed 2019 UA185 within 4"):
   its orbits, propagated two-body, drift far from it — 155" for a 2019 ZTF
   image and a 2026 file — and SkyBoT's integrated ones are used instead;
6. check every candidate left against the known objects within 10' of it
   ([`nearby_known_objects`](@ref)): one there moving the same way is
   most likely a known object with a poor orbit — and, for recent images,
   against the MPC's NEO Confirmation Page ([`neocp_check`](@ref));
7. list the known objects expected in the field
   ([`known_objects_in_field`](@ref)), marking which were detected;
8. write the report: the unidentified tracklets' observations (for vetting in
   Astrometrica — never to submit unchecked), or, with none, the blank report
   for "no moving objects detected" ([`blank_mpc_report`](@ref)). `header`
   holds [`mpc_report_header`](@ref)'s keywords (e.g. from
   [`astrometrica_settings`](@ref)).

`progress(stage, fraction)` is called as it goes (`stage` a `Symbol`:
`:extract`, `:clean`, `:refine`, `:search`, `:identify`, `:nearby`, `:field`,
`:report`, `:done`). `name` names a set given as files (default: their
folder's name). Returns a named tuple: `name`, `paths` (the original
frames), `epochs`, `wcs` (refined, per frame), `image_size`, `ratios`,
`frames_used`, `dropped`, `reliable`, `candidates` (an
[`astrometric_calibrate`](@ref) table with `frame` indexing `paths`, and a
`mag` column), `known` (`Dict` tracklet id => name), `unknown_ids`, `nearby`
(`nearby_known_objects`'s check of each unknown tracklet), `neocp`
(`neocp_check`'s),
`expected` (from `known_objects_in_field`, each with a `detected` flag),
`report`, `report_kind` (`:blank`, `:candidates` or `:unreliable`) and
`warnings`.
"""
function analyze_set(source; workdir::AbstractString, observatory::AbstractString, mpcorb=nothing,
                     header=(;), name=nothing,
                     extract_to::AbstractString=dirname(abspath(source isa AbstractString ? source : first(source))),
                     config=nothing, progress=(stage, fraction) -> nothing)
    warnings = String[]
    progress(:extract, 0.0)
    set_name, paths = if source isa AbstractVector
        (basename(dirname(abspath(first(source)))), collect(String, source))
    elseif isdir(source)
        (basename(rstrip(abspath(source), '/')), _fits_in(source))
    else
        extract_set(source, extract_to)
    end
    name = something(name, set_name)
    length(paths) >= 3 || throw(ArgumentError("$name: $(length(paths)) frame(s); a set needs at least 3"))
    paths = paths[sortperm([FITS(p -> frame_epoch(read_header(p[1])), q, "r") for q in paths])]

    progress(:clean, 0.05)
    cleaned = [clean_blank_pixels(p; outdir=joinpath(workdir, "cleaned", name)) for p in paths]

    header_wcs = WCSTransform[]; epochs = Float64[]; image_size = (0, 0)
    images = Matrix{Float64}[]; masks = BitMatrix[]
    stars = nothing
    for (k, p) in enumerate(cleaned)
        progress(:refine, 0.05 + 0.2 * (k - 1) / length(cleaned))
        image, hstr, hdr = FITS(p, "r") do f
            permutedims(Float64.(read(f[1]))), read_header(f[1], String), read_header(f[1])
        end
        image_size = size(image)
        epoch = frame_epoch(hdr)
        w = load_wcs(hstr)
        stars === nothing && (stars = _field_reference_stars(image, w, epoch))
        push!(images, image); push!(masks, fill_value_mask(image))
        push!(header_wcs, w); push!(epochs, epoch)
    end
    progress(:refine, 0.25)
    refined = refine_wcs_set(images, header_wcs, stars; masks)
    images = masks = nothing
    for (k, r) in enumerate(refined)
        r.shift_only && push!(warnings, "frame $k: too few Gaia stars; only the header WCS's offset corrected")
        r.refined || r.shift_only || push!(warnings, "frame $k: astrometric refinement failed; header WCS kept")
    end
    wcs = [r.wcs for r in refined]
    ratios = [r.n_detections / max(length(stars), 1) for r in refined]
    mpcorb = _mpcorb_for(mpcorb, epochs, warnings)
    good = findall(<=(_MAX_DETECTIONS_PER_REFERENCE_STAR), ratios)
    reliable = length(good) >= 3
    used = reliable ? good : collect(eachindex(paths))
    dropped = setdiff(eachindex(paths), used)

    candidates = Table(id=Int[], frame=Int[], x=Int[], y=Int[], ra=Float64[], dec=Float64[], epoch=Float64[],
                       mag=Float64[])
    known = Dict{Int,String}()
    if reliable
        progress(:search, 0.3)
        cfg = something(config, pipeline_config(_pixel_scale_arcsec(wcs[used[1]], image_size)))
found = run_pipeline(cleaned[used]; cfg..., reference_stars=stars)
        if !isempty(found)
            mags = candidate_magnitudes(cleaned[used], found; reference_stars=stars)
            candidates = Table(id=collect(found.id), frame=used[found.frame], x=collect(found.x), y=collect(found.y),
                               ra=collect(found.ra), dec=collect(found.dec), epoch=collect(found.epoch), mag=mags)
            progress(:identify, 0.75)
            known = _identify_tracklets(candidates, mpcorb, observatory, min(3, length(used)), warnings)
        end
    end

    # every candidate left: is a known object near it, moving the same way?
    nearby = Dict{Int,Any}()
    neocp = Dict{Int,Any}()
    unknown_rows = filter(r -> !haskey(known, r.id), collect(candidates))
    if !isempty(unknown_rows)
        progress(:nearby, 0.8)
        check = nearby_known_objects(unknown_rows; observatory, mpcorb)
        nearby = check.checks
        append!(warnings, check.warnings)
        on_neocp = neocp_check(unknown_rows; observatory)
        neocp = on_neocp.checks
        append!(warnings, on_neocp.warnings)
    end

    progress(:field, 0.85)
    field = known_objects_in_field(wcs, epochs, image_size; mpcorb, observatory)
    append!(warnings, field.warnings)
    detected_names = Set(_bare_name(n) for n in values(known))
    expected = [merge(o, (detected=_bare_name(o.name) in detected_names,)) for o in field.objects]

    progress(:report, 0.95)
    unknown_ids = sort(setdiff(unique(candidates.id), keys(known)))
    report, report_kind = if !reliable
        blank_mpc_report(observatory; header...), :unreliable
    elseif isempty(unknown_ids)
        blank_mpc_report(observatory; header...), :blank
    else
        mpc_report(Table(filter(r -> r.id in unknown_ids, collect(candidates))), observatory;
                   trksub_prefix="P", header...), :candidates
    end
    progress(:done, 1.0)
    return (name=name, paths=paths, epochs=epochs, wcs=wcs, image_size=image_size, ratios=ratios,
            frames_used=used, dropped=dropped, reliable=reliable, candidates=candidates, known=known,
            unknown_ids=unknown_ids, nearby=nearby, neocp=neocp, expected=expected, report=report, report_kind=report_kind,
            warnings=warnings)
end

# Tracklet id => known object's name, for tracklets SkyBoT (5") or MPCORB
# (15") match in at least `min_frames` of their rows — a static star near a
# known asteroid's track matches in one frame, not three. MPCORB's name wins
# (it carries the number) when both catalogs identify a tracklet.
# MPCORB's orbits osculate at one recent epoch; propagated two-body far from
# it they drift (155" for a 2019 ZTF image and a 2026 file), so beyond
# _MPCORB_MAX_DAYS the orbit file is set aside, with a warning.
const _MPCORB_MAX_DAYS = 365

function _mpcorb_for(mpcorb, epochs, warnings)
    mpcorb === nothing && return nothing
    orbits = mpcorb isa AbstractString ? read_mpcorb(mpcorb) : mpcorb
    isempty(orbits) && return orbits
    epoch = median(o.epoch for o in Iterators.take(orbits, 1000))
    days = maximum(abs(e - epoch) for e in epochs)
    days <= _MPCORB_MAX_DAYS && return orbits
    push!(warnings, "MPC orbit file not used: its orbits are for $(round(Int, days)) days from these " *
                    "images, too far for two-body propagation; known objects from SkyBoT")
    return nothing
end

function _identify_tracklets(candidates, mpcorb, observatory, min_frames, warnings)
    rows = collect(candidates)
    sources = Any[(:skybot, (; radius=5.0, observatory))]
    mpcorb === nothing || pushfirst!(sources, (:mpcorb, (; radius=15.0, mpcorb, observatory)))
    known = Dict{Int,String}()
    for (catalog, kw) in sources
        matches = try
            crossmatch_catalog(rows, catalog; kw...)
        catch e
            (_is_service_failure(e) || e isa HTTP.Exceptions.HTTPError || e isa SkyBoTServiceError ||
             e isa Base.IOError) || rethrow()
            push!(warnings, "$catalog unavailable: $(sprint(showerror, e))")
            continue
        end
        hits = Dict{Tuple{Int,String},Int}()
        for m in matches                        # one match per row and object
            hits[(m.id, m.name)] = get(hits, (m.id, m.name), 0) + 1
        end
        for ((id, name), n) in sort(collect(hits); by=first)
            n >= min_frames && !haskey(known, id) && (known[id] = name)
        end
    end
    return known
end

"""
    frame_footprint(wcs_per_frame, image_size; margin=0) -> Function

`(ra, dec) -> Bool`: whether a sky position falls on any of the frames
(`image_size = (ny, nx)`, `margin` pixels of slack) — the `field` test
[`check_mpc_report`](@ref) takes.
"""
function frame_footprint(wcs_per_frame::AbstractVector, image_size; margin::Real=0)
    ny, nx = image_size
    return function (ra, dec)
        any(wcs_per_frame) do w
            x, y = world_to_pix(w, Float64[ra, dec])
            1 - margin <= x <= nx + margin && 1 - margin <= y <= ny + margin
        end
    end
end
