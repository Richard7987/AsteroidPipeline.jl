#=
The interface's analysis engine: a long-lived Julia process that loads
AsteroidPipeline once and then answers commands, one JSON object per line on
stdin, with events, one JSON object per line on stdout. Keeping it alive
means paying Julia's start-up and compilation once per session, not once per
image set.

Commands (every one carries an `id`, echoed in its events):

    {"cmd": "inspect", "id": 0, "source": "/path/set.zip"}
        -> {"event": "inspection", "telescope": "PS1", "instrument": "gpc1", "lon": -156.26, ...,
            "candidates": [{"code": "F51", "name": "Pan-STARRS 1, Haleakala", ...}, ...]}
    {"cmd": "analyze", "id": 1, "source": "/path/set.zip", "out": "/cache/set",
     "observatory": "F51", "name": "set", "header": {"observers": "...", ...}}
        -> {"event": "progress", "stage": "refine", "fraction": 0.12}, ...
        -> {"event": "result", "path": "/cache/set/result.json"}
    {"cmd": "build_report", "id": 2, "result": "/cache/set/result.json",
     "tracklets": [3, 7], "header": {...}, "designations": {"3": "NHU0002"}}
        -> {"event": "report", "text": "COD F51..."}
    {"cmd": "check_report", "id": 3, "text": "COD F51...", "observatory": "F51",
     "result": "/cache/set/result.json"}
        -> {"event": "report_check", "ok": false, "issues": [...], "cleaned": "COD F51..."}

"source" is a .zip, a folder of FITS frames or a list of FITS paths; "name",
"header" (observers, measurers, telescope, contact — the report's header) and
"result" are optional.

Any command can end in {"event": "error", "message": "..."}. Unprompted:
{"event": "ready"} once loaded, and {"event": "status", "message": "..."}.
Everything else the libraries print goes to stderr, which the interface
keeps in its log file.
=#
const PROTOCOL = stdout
redirect_stdout(stderr)

using AsteroidPipeline
using AsteroidPipeline: _angular_distance_arcsec
using FITSIO, WCS, JSON, Dates, Statistics
using TypedTables: Table

# Side of the pixel cut-outs kept around every tracklet, for training a
# classifier later: odd, so the object sits on the centre pixel.
const STAMP_SIZE = 63

function emit(; kw...)
    println(PROTOCOL, JSON.json(Dict(string(k) => v for (k, v) in kw)))
    flush(PROTOCOL)
end

# The MPC's orbit file, loaded on first use: today's copy if it can be
# fetched, else the last one fetched (download_mpcorb's own fallback), else
# none — known objects then come from SkyBoT alone.
const ORBITS = Ref{Any}(nothing)
const ORBITS_INFO = Ref{Dict{String,Any}}(Dict{String,Any}())
function orbits(id)
    ORBITS[] === nothing || return ORBITS[] === false ? nothing : ORBITS[]
    try
        emit(; id, event="status", message="mpcorb_download")
        path = download_mpcorb()
        emit(; id, event="status", message="mpcorb_load")
        ORBITS[] = read_mpcorb(path)
        ORBITS_INFO[] = Dict{String,Any}("path" => path, "updated" => string(unix2datetime(mtime(path))))
    catch e
        @warn "no MPCORB.DAT available; known objects from SkyBoT only" exception = e
        ORBITS[] = false
        ORBITS_INFO[] = Dict{String,Any}("path" => nothing)
    end
    return ORBITS[] === false ? nothing : ORBITS[]
end

# The report header's fields from the interface's settings; the reference
# catalog is always Gaia DR3, the one positions are calibrated against.
function report_header(h)
    h === nothing && return (net="Gaia DR3",)
    get_(k) = String(something(get(h, k, ""), ""))
    return (observers=get_("observers"), measurers=get_("measurers"), telescope=get_("telescope"),
            contact=get_("contact"), net="Gaia DR3")
end

# Frames for display, all on frame 1's pixel grid (integer shifts: the
# telescope's dithers between exposures), stretched to 8 bits, top row first
# (FITS row 1 is the bottom), as binary PGM files, which Pillow reads as is.
function write_display_frames(images, wcs, out)
    w1 = wcs[1]
    shifts = Tuple{Int,Int}[]
    for (k, (img, w)) in enumerate(zip(images, wcs))
        nx, ny = size(img)
        c = Float64[nx / 2, ny / 2]
        q = world_to_pix(w1, pix_to_world(w, c))
        dx, dy = round(Int, q[1] - c[1]), round(Int, q[2] - c[2])
        push!(shifts, (dx, dy))
        sample = filter(isfinite, vec(img)[1:7:end])
        med = median(sample)
        σ = 1.4826 * median(abs.(sample .- med))
        lo, hi = med - 1.5σ, med + 25σ
        stretch(v) = isfinite(v) ? round(UInt8, 255 * clamp(asinh(10 * (v - lo) / (hi - lo)) / asinh(10), 0, 1)) : 0x00
        buf = zeros(UInt8, nx, ny)
        for y in 1:ny, x in 1:nx
            xs, ys = x - dx, y - dy                               # display pixel -> this frame's pixel
            (1 <= xs <= nx && 1 <= ys <= ny) && (buf[x, y] = stretch(img[xs, ys]))
        end
        open(joinpath(out, "frame_$k.pgm"), "w") do io
            write(io, "P5\n$nx $ny\n255\n")
            for y in ny:-1:1
                write(io, @view buf[:, y])
            end
        end
    end
    return shifts
end

# One FITS cube per tracklet, STAMP_SIZE x STAMP_SIZE x frames, of the cleaned
# frames' own pixel values (no stretch, no alignment: each frame centred on
# the object's position in it through that frame's WCS). Frames the tracklet
# has no point in are centred on its position interpolated in time, so every
# cube has all frames — what a classifier needs to see the object move.
function write_stamps(images, wcs, epochs, tracklet_points, path)
    h = STAMP_SIZE ÷ 2
    cube = fill(NaN32, STAMP_SIZE, STAMP_SIZE, length(images))
    ts = [p.epoch for p in tracklet_points]
    for (k, (img, w)) in enumerate(zip(images, wcs))
        # nearest two points in time, linear in RA/Dec (fine over an hour)
        j = clamp(searchsortedfirst(ts, epochs[k]), 2, max(2, length(ts)))
        a, b = tracklet_points[max(1, j - 1)], tracklet_points[min(j, end)]
        f = b.epoch == a.epoch ? 0.0 : (epochs[k] - a.epoch) / (b.epoch - a.epoch)
        ra, dec = a.ra + f * (b.ra - a.ra), a.dec + f * (b.dec - a.dec)
        x, y = round.(Int, world_to_pix(w, Float64[ra, dec]))
        nx, ny = size(img)
        for dy in -h:h, dx in -h:h
            xs, ys = x + dx, y + dy
            (1 <= xs <= nx && 1 <= ys <= ny) && (cube[dx + h + 1, dy + h + 1, k] = img[xs, ys])
        end
    end
    mkpath(dirname(path))
    FITS(path, "w") do f
        write(f, cube)
    end
end

label(id) = "P" * uppercase(string(id; base=36))
utc(jd) = Dates.format(julian2datetime(jd), "yyyy-mm-dd HH:MM:SS")

function analyze(id, msg)
    source = msg["source"] isa AbstractVector ? collect(String, msg["source"]) : String(msg["source"])
    out = msg["out"]
    observatory = msg["observatory"]
    mkpath(out)
    mp = orbits(id)
    r = analyze_set(source; workdir=joinpath(out, "work"), observatory, name=get(msg, "name", nothing),
                    mpcorb=mp, header=report_header(get(msg, "header", nothing)),
                    progress=(stage, fraction) -> emit(; id, event="progress", stage=string(stage), fraction))
    emit(; id, event="progress", stage="images", fraction=1.0)
    cleaned = [AsteroidPipeline.clean_blank_pixels(p; outdir=joinpath(out, "work", "cleaned", r.name)) for p in r.paths]
    images = [FITS(f -> Float32.(read(f[1])), p, "r") for p in cleaned]   # (x, y)
    shifts = write_display_frames(images, r.wcs, out)
    w1 = r.wcs[1]
    px(ra, dec) = world_to_pix(w1, Float64[ra, dec])

    rows = collect(r.candidates)
    tracklets = Dict{String,Any}[]
    for tid in sort(unique(r.candidates.id))
        pts = sort(filter(x -> x.id == tid, rows); by=x -> x.epoch)
        motion = tracklet_motion(pts)
        mags = filter(isfinite, [p.mag for p in pts])
        stamp = joinpath("stamps", "tracklet_$tid.fits")
        write_stamps(images, r.wcs, r.epochs, pts, joinpath(out, stamp))
        check = get(r.nearby, tid, nothing)
        neocp = get(r.neocp, tid, nothing)
        push!(tracklets, Dict{String,Any}(
            "id" => tid, "label" => label(tid), "known" => get(r.known, tid, nothing),
            "speed_arcsec_h" => motion.rate_arcsec_h, "pa_deg" => motion.pa_deg, "rms_arcsec" => motion.rms_arcsec,
            "mag" => isempty(mags) ? nothing : median(mags), "stamp" => stamp,
            "check" => check === nothing ? nothing : Dict(
                "verdict" => string(check.verdict),
                "objects" => [Dict("name" => o.name, "source" => o.source, "separation_arcsec" => o.separation_arcsec,
                                   "mag" => o.mag, "rate_arcsec_h" => o.rate_arcsec_h, "pa_deg" => o.pa_deg,
                                   "agreement" => string(o.agreement)) for o in check.objects]),
            "neocp" => neocp === nothing ? nothing : Dict(
                "verdict" => string(neocp.verdict),
                "objects" => [Dict("designation" => o.designation, "score" => o.score,
                                   "separation_arcsec" => o.separation_arcsec, "rate_arcsec_h" => o.rate_arcsec_h,
                                   "pa_deg" => o.pa_deg, "uncertain" => o.uncertain,
                                   "agreement" => string(o.agreement)) for o in neocp.objects]),
            "points" => [Dict("frame" => p.frame, "ra" => p.ra, "dec" => p.dec, "epoch" => p.epoch,
                              "mag" => isfinite(p.mag) ? p.mag : nothing,
                              "x" => px(p.ra, p.dec)[1], "y" => px(p.ra, p.dec)[2]) for p in pts]))
    end
    bare(n) = strip(replace(n, r"^\(\d+\)\s*" => ""))
    known = [Dict{String,Any}(
        "name" => o.name, "packed" => o.packed, "mag" => o.mag, "source" => o.source, "detected" => o.detected,
        "tracklet" => findfirst(n -> bare(n) == bare(o.name), r.known),
        "points" => [Dict("frame" => k, "epoch" => o.epochs[k], "ra" => o.ra[k], "dec" => o.dec[k],
                          "x" => px(o.ra[k], o.dec[k])[1], "y" => px(o.ra[k], o.dec[k])[2]) for k in eachindex(o.ra)])
        for o in r.expected]

    ny, nx = r.image_size
    result = Dict{String,Any}(
        "name" => r.name, "observatory" => observatory, "set_dir" => dirname(r.paths[1]),
        "source" => source isa AbstractVector ? abspath.(source) : abspath(source),
        "analyzed_at" => string(now()), "image_size" => [nx, ny], "reliable" => r.reliable,
        "report_kind" => string(r.report_kind), "report" => r.report, "warnings" => r.warnings,
        "mpcorb" => ORBITS_INFO[], "stamp_size" => STAMP_SIZE,
        "frames" => [Dict("index" => k, "file" => basename(p), "epoch" => r.epochs[k], "time_utc" => utc(r.epochs[k]),
                          "ratio" => r.ratios[k], "used" => k in r.frames_used, "image" => "frame_$k.pgm",
                          "shift" => collect(shifts[k])) for (k, p) in enumerate(r.paths)],
        "tracklets" => tracklets, "known" => known,
        "wcs_headers" => [WCS.to_header(w) for w in r.wcs])
    path = joinpath(out, "result.json")
    write(path * ".part", JSON.json(result))
    mv(path * ".part", path; force=true)
    emit(; id, event="result", path)
end

# The report for the tracklets a person confirmed as real: their
# observations, or the blank report if there are none.
function build_report(id, msg)
    res = JSON.parse(read(msg["result"], String))
    wanted = Set(Int.(msg["tracklets"]))
    rows = [(id=Int(t["id"]), frame=Int(p["frame"]), ra=Float64(p["ra"]), dec=Float64(p["dec"]),
             epoch=Float64(p["epoch"]), mag=p["mag"] === nothing ? NaN : Float64(p["mag"]))
            for t in res["tracklets"] if Int(t["id"]) in wanted for p in t["points"]]
    header = report_header(get(msg, "header", nothing))
    given = get(msg, "designations", nothing)       # {"<tracklet id>": "NHU0002", ...}
    designations = given === nothing || isempty(given) ? nothing :
        Dict(id => get(given, string(id), "P" * uppercase(string(id; base=36))) for id in wanted)
    text = isempty(rows) ? blank_mpc_report(res["observatory"]; header...) :
           mpc_report(Table(rows), res["observatory"]; trksub_prefix="P", designations, header...)
    emit(; id, event="report", text)
end

function check_report(id, msg)
    kw = (;)
    result_path = get(msg, "result", nothing)
    if result_path !== nothing && isfile(result_path)
        res = JSON.parse(read(result_path, String))
        wcs = [WCS.from_header(h)[1] for h in res["wcs_headers"]]
        nx, ny = res["image_size"]
        known = [(name=k["name"], packed=k["packed"], epochs=[p["epoch"] for p in k["points"]],
                  ra=[p["ra"] for p in k["points"]], dec=[p["dec"] for p in k["points"]]) for k in res["known"]]
        kw = (epochs=[f["epoch"] for f in res["frames"]], field=frame_footprint(wcs, (ny, nx); margin=20), known=known)
    end
    c = check_mpc_report(msg["text"]; observatory=msg["observatory"], kw...)
    emit(; id, event="report_check", ok=c.ok, cleaned=c.cleaned,
         issues=[Dict("code" => string(i.code), "line" => i.line, "detail" => i.detail) for i in c.issues])
end

# What the frames say about where they were taken, and the MPC codes that
# fits — before analyzing, so the interface can pick the telescope's profile
# (or ask, the first time). A .zip is unpacked here, as analyze would.
function inspect(id, msg)
    source = msg["source"]
    paths = if source isa AbstractVector
        collect(String, source)
    elseif isdir(source)
        AsteroidPipeline._fits_in(source)
    else
        extract_set(source, dirname(abspath(source))).paths
    end
    isempty(paths) && error("no FITS frames in $source")
    site = FITS(f -> frame_site(read_header(f[1])), first(paths), "r")
    candidates = try
        first(observatory_candidates(site, read_obscodes(download_obscodes())), 6)
    catch e
        @warn "no MPC observatory list" exception = e
        []
    end
    emit(; id, event="inspection", n_frames=length(paths), telescope=site.telescope, instrument=site.instrument,
         origin=site.origin, observatory=site.observatory, lon=site.lon, lat=site.lat, elev=site.elev,
         candidates=[Dict("code" => c.code, "name" => c.name, "score" => c.score,
                          "distance_km" => isnan(c.distance_km) ? nothing : c.distance_km) for c in candidates])
end

const COMMANDS = Dict("inspect" => inspect, "analyze" => analyze, "build_report" => build_report,
                      "check_report" => check_report)

emit(event="ready")
for line in eachline(stdin)
    isempty(strip(line)) && continue
    id = nothing
    try
        msg = JSON.parse(line)
        id = get(msg, "id", nothing)
        cmd = get(msg, "cmd", "")
        cmd == "quit" && break
        haskey(COMMANDS, cmd) || error("unknown command: $cmd")
        COMMANDS[cmd](id, msg)
    catch e
        @error "command failed" exception = (e, catch_backtrace())
        emit(; id, event="error", message=sprint(showerror, e))
    end
end
