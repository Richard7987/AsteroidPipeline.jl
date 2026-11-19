#=
Whether raw frames reduced here (plain_frame, calibrate_frames) give what
the observatory's own reduction gives:

    julia --project=. examples/calibration_test.jl

Declared before running:

Data. Las Cumbres Observatory's public archive: every set of the nights of
2024-10-01, -04 and -09 (the archive's listings of 2024-10-02, -05 and -10
UTC, surveyed for this) whose target is a minor planet and that has at
least 5 exposures of 15 s or more in one filter on one camera (the filter
with the most, ties alphabetically): (19)
Fortuna, (109) Felicitas and (187) Lamberta on 0.4 m telescopes (QHY600
cameras, one image per file) and (203) Pompeja on 1 m telescopes (Sinistro
cameras, four amplifiers per file, each with its overscan), 14 sets. Of
each, the first 5 such exposures, raw (tile-compressed); the same 5 as the
archive's pipeline (BANZAI) reduced them; and raw calibration frames of the
same camera from the same night, or the nearest of the 3 before it: the
first 5 bias frames, 3 darks, 5 sky flats in the filter. A set without
flats in those nights is left out.

Three reductions of each set, each searched by analyze_set with its
defaults, the MPC code nearest the headers' site, the MPC's orbits:
"raw" (plain_frame only: the amplifiers assembled, overscan subtracted),
"here" (calibrate_frames with the calibration frames), "BANZAI".

Measured, per reduction: whether analyze_set marks the set reliable; the
astrometric fit's scatter; the limiting magnitude; the catalogued asteroids
identified, the target among them; their positions against JPL Horizons'
from the headers' site (median error, mean RA and Dec error, errors over
pos_err); the target's magnitudes, here against BANZAI's, over their
mag_err combined; the unidentified tracklets.

Expected if the reduction here is right: as many sets reliable as with
BANZAI, the same asteroids identified, positions as close to JPL's, limiting
magnitudes within 0.1 of BANZAI's, magnitudes differing by no more than
their uncertainties say (deviation of the ratios about 1 or less).

Frames are kept in data/real/lco, results cached in data/real/calibration_test.

Changed after the first 4 sets: calibration frames were chosen by camera,
kind, filter and night only, and the 0.4 m cameras read some in another
mode (the sky flats of Felicitas sq31 2024-10-04 and Lamberta sq33
2024-10-01 full frame, 9600x6422, the exposures a central 2400x2400), so
calibrate_frames could not combine them and both sets failed. Calibration
frames are now only those read in the exposures' mode (the archive's
headers: CONFMODE and CCDSUM equal). The two sets done before the change
(Fortuna sq34, Felicitas sq34) had their calibrations in that mode already.

Run 1 (14 sets; 3 left out: Pompeja fa01 and fa06 without Y flats, fa15 by
a Gaia query refused, 401): all 11 reliable in each reduction, but the
search tested little. The 0.4 m sets' first 5 exposures span 1-2 min, so
any main-belt asteroid is subtracted with the stars; on 1 m frames Pompeja
moves near the slowest rate 14 min can show. Any tracklet at all came
from 4 of the 11 sets with BANZAI's frames, 4 raw and 3 here; Pompeja was
found in 3 sets with BANZAI's frames, 1 with each of the others; the frames
here gave 61 unidentified tracklets on one set (fa14) whose astrometric fit
had failed. The fit failed, or only corrected an offset, on every frame of
6 sets raw and 3 here, against 1 frame of 1 set with BANZAI's: raw
headers hold the telescope's pointing, 9-106" off and, on the 4 sets
measured, turned 0.37-0.67 deg and 0.15-0.6% out in scale, and BANZAI's
frames carry its own solution. Taken from the raw headers, BANZAI's
pixels failed as the frames here did, and with BANZAI's solution the
frames here matched as many stars as BANZAI's own: the failures were the
astrometry, not the reduction. refine_wcs_set now also solves a frame from
its stars when the header's offset alone fails (solve_wcs near the
header's centre and scale).

Run 2, declared before it was run (julia --project=. examples/calibration_test.jl 2):
the same archive and rules on other nights, 2024-10-15, -18, -22, -25, -28
and 11-01, which gives 6 sets of minor planets: (19) Fortuna sq34,
(450649) fa08, (201) Penelope sq30, (617) Patroclus sq35, 2024 TX13 fa08
and (3460) Shawna sq40; the comet C/2023 A3 is not one. Of each, 5
exposures spread over the sequence's first hour (the first, the last
within the hour and 3 evenly between by number), so that moving objects
move. Calibration frames in the exposures' readout mode, as above.
Measured besides the search, frame by frame and for each reduction, with
the frames' raw headers' WCS for all three:
- the astrometric fit: frames refined, stars matched, scatter, and the
  solution's largest difference from BANZAI's over the frame's corners and
  centre;
- the photometric zero point against Gaia (photometric_zeropoint, through
  BANZAI's WCS) and its scatter over the frame's stars, the error a flat
  field leaves;
- the sky's signal to noise in the central eighth, over BANZAI's.
A set whose catalogue or ephemeris services fail is tried again, twice,
5 min apart. Expected: refined frames wherever BANZAI's pixels are, with
solutions within 0.3" of BANZAI's; zero-point scatter here no more than
0.01 mag above BANZAI's, and the raw frames' above both (or the measure
cannot tell a flat's work); sky signal to noise at least 0.9 of BANZAI's;
for the search as in run 1. Results cached in data/real/calibration_test2.

Changed in run 2 before any of its results: the first launch failed on 4 of
the 6 sets with errors, not results (the script lacked `using WCS`; two of
Patroclus's exposures have no BANZAI product, so only exposures BANZAI
reduced are chosen; 2024 TX13 was read out 2x2 binned from a subarray,
whose DETSEC plain_frame did not scale, now fixed). The other two were left
out by the rules above (no w flats for (450649); Shawna fewer than 5
exposures in an hour).

Run 2 (4 sets): all reliable in each reduction. Frames refined 19 of 20 in
each, the same one (Penelope's first) failing in all three. Solutions a
median 0.74" and at most 1.25" from BANZAI's, and BANZAI's own pixels,
from the same headers, 0.75" and 1.24": the 0.3" expected was wrongly
framed, BANZAI's solution is not the truth (positions against JPL are).
Zero-point scatter a median 0.041 mag here, 0.033 BANZAI's, 0.059 raw;
here less BANZAI's a median -0.000, at most +0.082 (one frame). Sky signal
to noise here 0.98-0.99 of BANZAI's. Limiting magnitudes here less
BANZAI's +0.03, -0.02, -0.01 (Fortuna: no tracklets in any reduction).
Catalogued asteroids: 9 with BANZAI's frames, 7 here, 7 raw; the target
identified in 1 of 4 sets in each (Patroclus). Against JPL: here 35
points, median error 0.070", errors over pos_err deviation 1.05 (robust
0.89); BANZAI's 45, 0.102", 1.98 (0.93); raw 29, 0.096". Patroclus's
magnitudes here less BANZAI's a median -0.001, deviation over the combined
mag_err 0.30. Of the 2 asteroids found with BANZAI's frames and not here:
2010 KW118 (Patroclus) has the same signal to noise at its positions in
both (4-8 in a 4 px aperture), is detected frame by frame in 1 frame here
and 2 in BANZAI's, and was found only by stacking: at the threshold. 2001
DQ96 (Penelope) went with the first frame, whose fit failed from the raw
header; that frame gave 0 detections here, 51 with BANZAI's frames, which
carry BANZAI's solution.

Changed after run 2, in refine_wcs_set: a frame left with few stars,
started again from the best frame's solution carried over by the headers,
was searched only within rematch_arcsec (2"), and raw headers' pointing
differs between exposures by more (Penelope's first frame 12" off, the
others 18-23"); an offset is now searched (refine_wcs's default) when that
fails. Penelope's first frame is then refined in all three reductions
(121-123 stars, 1.25" from BANZAI's, as the others). And the joint fit
gave a frame without stars of its own the others' through the headers,
unchecked: on a synthetic pair with headers 15" apart that frame came out
refined and 15" off. The joint fit now needs 3 stars of each frame's own.
Run 2 is kept as it was; both changes are checked on new data in run 3.

Run 3, declared before it was run (julia --project=. examples/calibration_test.jl 3):
the same rules on the nights of 2024-11-05, -08, -12, -15 and -19 (11-22's
listing is longer than the archive pages through), a target a minor planet
if JPL's Small-Body Database calls it an asteroid: 13 sets, 1769
Carlostorres ep06 and ep07 (MuSCAT, 2 m), Shawna sq37 and Xanthippe sq30 on
11-05; 2024 VC ep02-ep05, (447755) fa12 and fa15 and Xanthippe sq30 on
11-08; Xanthippe sq31 on 11-12 and on 11-19. The comets C/2022 L2,
C/2023 A3, 29P and P/2024 R3 are not. Measured as in run 2. Expected:
refined frames wherever BANZAI's pixels are, solutions no further from
BANZAI's than BANZAI's pixels' from the same headers plus 0.1"; positions
against JPL no worse than with BANZAI's frames (median); zero point, sky and
limiting magnitude as in run 2; the same catalogued asteroids, any found
with one reduction and not the other examined by the signal to noise at its
positions in both. Results cached in data/real/calibration_test3.
=#
using AsteroidPipeline, FITSIO, HTTP, JSON, Printf, Serialization, Statistics, WCS
const P = AsteroidPipeline

const API = "https://archive-api.lco.global/frames/"
const RUN = isempty(ARGS) ? 2 : parse(Int, ARGS[1])
const SETS1 = [("Fortuna", "sq34", "2024-10-09"), ("Felicitas", "sq34", "2024-10-01"), ("Felicitas", "sq31", "2024-10-04"),
              ("Lamberta", "sq33", "2024-10-01"), ("Lamberta", "sq33", "2024-10-04"),
              ("Pompeja", "fa01", "2024-10-01"), ("Pompeja", "fa04", "2024-10-01"), ("Pompeja", "fa06", "2024-10-01"),
              ("Pompeja", "fa08", "2024-10-01"), ("Pompeja", "fa12", "2024-10-02"), ("Pompeja", "fa14", "2024-10-01"),
              ("Pompeja", "fa15", "2024-10-01"), ("Pompeja", "fa16", "2024-10-01"), ("Pompeja", "fa19", "2024-10-02")]
const SETS2 = [("Fortuna", "sq34", "2024-10-15"), ("450649", "fa08", "2024-10-15"), ("Penelope", "sq30", "2024-10-22"),
               ("Patroclus", "sq35", "2024-10-25"), ("2024 TX13", "fa08", "2024-10-28"), ("Shawna", "sq40", "2024-11-01")]
const SETS3 = [("1769 Carlostorres", "ep06", "2024-11-05"), ("1769 Carlostorres", "ep07", "2024-11-05"), ("Shawna", "sq37", "2024-11-05"),
               ("Xanthippe", "sq30", "2024-11-05"), ("2024 VC", "ep02", "2024-11-08"), ("2024 VC", "ep03", "2024-11-08"),
               ("2024 VC", "ep04", "2024-11-08"), ("2024 VC", "ep05", "2024-11-08"), ("447755", "fa12", "2024-11-08"),
               ("447755", "fa15", "2024-11-08"), ("Xanthippe", "sq30", "2024-11-08"), ("Xanthippe", "sq31", "2024-11-12"),
               ("Xanthippe", "sq31", "2024-11-19")]
const SETS = (SETS1, SETS2, SETS3)[RUN]
const DATA = joinpath(@__DIR__, "..", "data", "real", "lco")
const CACHE = joinpath(@__DIR__, "..", "data", "real", RUN == 1 ? "calibration_test" : "calibration_test$RUN")
const PROGRESS = joinpath(homedir(), ".claude", "progress", "a2-calibration")
const ARMS = ("raw", "here", "BANZAI")

get_json(url) = JSON.parse(String(HTTP.get(url; retry=true, retries=8, read_idle_timeout=120).body))

function query(; kw...)
    out = Any[]
    url = API * "?" * HTTP.URIs.escapeuri(Dict(string(k) => string(v) for (k, v) in kw))
    while url !== nothing
        page = get_json(url)
        append!(out, page["results"])
        url = page["next"]
    end
    return sort(out; by=r -> r["observation_date"])
end

function fetch(frame, dir)
    mkpath(dir)
    path = joinpath(dir, frame["filename"])
    isfile(path) && return path
    url = frame["url"]
    write(path * ".part", HTTP.get(url; retry=true, retries=8, read_idle_timeout=300).body)
    mv(path * ".part", path)
    return path
end

day_after(d) = string(P.Dates.Date(d) + P.Dates.Day(1))
nights_back(d) = [string(P.Dates.Date(d) - P.Dates.Day(k)) for k in 0:3]

# raw science, BANZAI's, and calibration frames of a set
function download_set(target, camera, night)
    exposures = filter(r -> r["exposure_time"] >= 15 && r["observation_day"] == night,
                       query(; public=true, reduction_level=0, configuration_type="EXPOSE", target_name=target,
                             instrument_id=camera, start=night, :end => day_after(day_after(night)), limit=100))
    by_filter = Dict{String,Vector{Any}}()
    for r in exposures
        push!(get!(by_filter, r["primary_optical_element"], Any[]), r)
    end
    band, frames = first(sort(collect(by_filter); by=kv -> (-length(kv[2]), kv[1])))
    length(frames) >= 5 || error("fewer than 5 exposures of 15 s or more in one filter")
    frames = RUN == 1 ? frames[1:5] : let
        t(r) = P.Dates.DateTime(first(r["observation_date"], 19))
        reduced(r) = any(id -> (x = get_json(API * "$(id)/"); x["reduction_level"] == 91 && x["configuration_type"] == "EXPOSE"),
                         r["related_frames"])
        hour = filter(r -> t(r) - t(frames[1]) <= P.Dates.Hour(1) && reduced(r), frames)
        length(hour) >= 5 || error("fewer than 5 exposures in the first hour")
        hour[round.(Int, range(1, length(hour); length=5))]
    end
    dir = joinpath(DATA, replace("$(target)_$(camera)_$(night)", " " => "_"))
    raw = [fetch(r, joinpath(dir, "raw")) for r in frames]
    banzai = map(frames) do r
        reduced = [get_json(API * "$(id)/") for id in r["related_frames"]]
        fetch(first(filter(x -> x["reduction_level"] == 91 && x["configuration_type"] == "EXPOSE", reduced)),
              joinpath(dir, "banzai"))
    end
    readout(r) = (h = get_json(API * "$(r["id"])/headers/")["data"]; (get(h, "CONFMODE", nothing), get(h, "CCDSUM", nothing)))
    mode = readout(frames[1])
    function calibration(kind, n; band=nothing)
        for d in nights_back(night)
            found = filter(r -> r["observation_day"] == d && (band === nothing || r["primary_optical_element"] == band),
                           query(; public=true, reduction_level=0, configuration_type=kind, instrument_id=camera,
                                 start=string(P.Dates.Date(d) - P.Dates.Day(1)), :end => day_after(day_after(d)),
                                 limit=100))
            found = first(Iterators.filter(r -> readout(r) == mode, found), n)
            isempty(found) || return [fetch(r, joinpath(dir, lowercase(kind))) for r in found]
        end
        return String[]
    end
    bias = calibration("BIAS", 5); dark = calibration("DARK", 3); flat = calibration("SKYFLAT", 5; band)
    isempty(flat) && error("no $band sky flats in the 3 nights before")
    return (band=band, raw=raw, banzai=banzai, bias=bias, dark=dark, flat=flat, dir=dir)
end

codes = read_obscodes(download_obscodes())
orbits = nothing

function site_of(path)
    h = FITS(f -> read_header(f[1]), plain_frame(path; outdir=mktempdir()), "r")
    s = frame_site(h)
    return (code=first(observatory_candidates(s, codes)).code, lon=s.lon, lat=s.lat, elev=s.elev)
end

function analyze(name, frames, site)
    work = mktempdir()
    r = analyze_set(frames; workdir=work, extract_to=work, observatory=site.code, mpcorb=orbits, name=name)
    rows = collect(r.candidates)
    used = r.frames_used
    return (reliable=r.reliable, limiting_mag=r.limiting_mag, warnings=r.warnings,
            known=[(object=name, rows=[(epoch=c.epoch, ra=c.ra, dec=c.dec, mag=c.mag, mag_err=c.mag_err, pos_err=c.pos_err)
                                       for c in rows if c.id == id])
                   for (id, name) in r.known],
            unknown=length(r.unknown_ids))
end

img2d(f) = f[findfirst(h -> h isa ImageHDU && ndims(h) == 2 && all(>(0), size(h)), collect(f))]
read_frame(p) = FITS(f -> (image=permutedims(Float64.(read(img2d(f)))), header=read_header(img2d(f), String)), p)
robust_std(x) = (v = filter(isfinite, vec(x)); 1.4826 * median(abs.(v .- median(v))))

# frame by frame, each reduction from the raw headers' WCS; BANZAI's WCS as the reference
function image_measures(s, here)
    scratch = mktempdir()
    raw_plain = [plain_frame(p; outdir=scratch) for p in s.raw]
    header_wcs = [P.load_wcs(read_frame(p).header) for p in raw_plain]
    banzai = [read_frame(p) for p in s.banzai]
    banzai_wcs = [P.load_wcs(b.header) for b in banzai]
    epoch = frame_epoch(FITS(f -> read_header(f[1]), raw_plain[1]))
    stars = P._field_reference_stars(banzai[1].image, banzai_wcs[1], epoch)
    function sky_snr(img)
        n1, n2 = size(img)
        c = img[n1 ÷ 2 - n1 ÷ 16 + 1:n1 ÷ 2 + n1 ÷ 16, n2 ÷ 2 - n2 ÷ 16 + 1:n2 ÷ 2 + n2 ÷ 16]
        median(filter(isfinite, vec(c))) / robust_std(c)
    end
    out = Dict{String,Any}()
    for (arm, images) in zip(ARMS, ([read_frame(p).image for p in raw_plain], [read_frame(p).image for p in here],
                                    [b.image for b in banzai]))
        masks = [fill_value_mask(img) for img in images]
        fits = refine_wcs_set(images, header_wcs, stars; masks)
        off = map(zip(fits, banzai_wcs, images)) do (r, w, img)
            r.refined || r.shift_only || return NaN
            n1, n2 = size(img)
            maximum(((1.0, 1.0), (n2, 1.0), (1.0, n1), (n2, n1), (n2 / 2, n1 / 2))) do (x, y)
                a, b = pix_to_world(r.wcs, Float64[x, y]), pix_to_world(w, Float64[x, y])
                hypot((a[1] - b[1]) * cosd(b[2]), a[2] - b[2]) * 3600
            end
        end
        zp = [photometric_zeropoint(img, w, stars; mask=m) for (img, w, m) in zip(images, banzai_wcs, masks)]
        out[arm] = (refined=[r.refined for r in fits], shift_only=[r.shift_only for r in fits],
                    matches=[r.n_matches for r in fits], rms=[r.rms_arcsec for r in fits], off=off,
                    zeropoint=[z.zeropoint for z in zp], zp_stars=[z.n_stars for z in zp], zp_scatter=[z.scatter for z in zp],
                    sky_snr=sky_snr.(images))
    end
    rm(scratch; recursive=true, force=true)
    return out
end

service_failure(e) = occursin(r"http status|IOError|timed? ?out|Timeout|EOFError|Connect"i, sprint(showerror, e))
# the archive names some targets "1769 Carlostorres", analyze_set "Carlostorres"
is_target(t, object) = occursin(replace(t, r"^\d+ (?=[A-Z][a-z])" => ""), object)

mkpath(CACHE); mkpath(dirname(PROGRESS))
results = map(enumerate(SETS)) do (k, (target, camera, night))
    write(PROGRESS, "$(k - 1) $(length(SETS)) calibracion\n")
    name = replace("$(target)_$(camera)_$(night)", " " => "_")
    cache = joinpath(CACHE, name * ".jls")
    isfile(cache) && return deserialize(cache)
    global orbits
    orbits === nothing && (orbits = read_mpcorb(download_mpcorb()))
    result = nothing
    for attempt in 1:(RUN == 1 ? 1 : 3)
        result = try
            s = download_set(target, camera, night)
            site = site_of(s.raw[1])
            here = calibrate_frames(s.raw; outdir=joinpath(s.dir, "here"), s.bias, s.dark, s.flat)
            images = RUN == 1 ? nothing : image_measures(s, here)
            arms = Dict(arm => analyze(name, frames, site)
                        for (arm, frames) in zip(ARMS, (s.raw, here, s.banzai)))
            @printf("%s: %s, %d bias, %d dark, %d flat; reliable %s\n", name, s.band, length(s.bias), length(s.dark),
                    length(s.flat), join([arms[a].reliable for a in ARMS], "/"))
            flush(stdout)
            (name=name, target=target, camera=camera, band=s.band, site=site, arms=arms, images=images)
        catch e
            println(name, ": left out, ", first(sprint(showerror, e), 200))
            flush(stdout)
            failed = (name=name, failed=first(sprint(showerror, e), 200))
            attempt < 3 && service_failure(e) && RUN != 1 && (sleep(300); continue)
            failed
        end
        break
    end
    serialize(cache, result)
    result
end
write(PROGRESS, "$(length(SETS)) $(length(SETS)) calibracion\n")

function jpl_radec(object, jds, site)
    m = match(r"^\((\d+)\)", object)
    command = m === nothing ? (occursin(r"^\d{4} ", object) ? "'DES=$(strip(object));'" : "'$(strip(object));'") : "'$(m[1]);'"
    q = Dict("COMMAND" => command, "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'coord@399'", "COORD_TYPE" => "GEODETIC",
             "SITE_COORD" => @sprintf("'%.6f,%.6f,%.4f'", site.lon, site.lat, something(site.elev, 0.0) / 1000),
             "QUANTITIES" => "'1'", "ANG_FORMAT" => "DEG", "TIME_TYPE" => "UT", "CSV_FORMAT" => "YES",
             "OBJ_DATA" => "NO", "EXTRA_PREC" => "YES")
    out = Tuple{Float64,Float64}[]
    for result in P._horizons(q, [@sprintf("%.7f", jd) for jd in jds])
        for line in split(strip(P._horizons_table(result)), '\n')
            f = strip.(split(line, ','))
            push!(out, (parse(Float64, f[4]), parse(Float64, f[5])))
        end
    end
    length(out) == length(jds) || error("Horizons gave $(length(out)) positions for $(length(jds)) times")
    return out
end

println()
done = filter(r -> !haskey(r, :failed), results)
for r in results
    haskey(r, :failed) && println(r.name, ": left out, ", r.failed)
end
robust(z) = 1.4826 * median(abs.(z .- median(z)))
for arm in ARMS
    errs = Float64[]; z = Float64[]; dra = Float64[]; ddec = Float64[]
    lims = Float64[]; nknown = 0; targets = 0; unknown = 0; reliable = 0
    for r in done
        a = r.arms[arm]
        reliable += a.reliable
        a.reliable || continue
        isfinite(a.limiting_mag) && push!(lims, a.limiting_mag)
        nknown += length(a.known); unknown += a.unknown
        targets += any(k -> is_target(r.target, k.object), a.known)
        for k in a.known
            truth = try
                jpl_radec(k.object, [p.epoch for p in k.rows], r.site)
            catch e
                println(r.name, " ", arm, " ", k.object, ": no ephemeris, ", first(sprint(showerror, e), 100))
                continue
            end
            for (p, (ra, dec)) in zip(k.rows, truth)
                x = (p.ra - ra) * cosd(dec) * 3600; y = (p.dec - dec) * 3600
                push!(dra, x); push!(ddec, y); push!(errs, x, y)
                isfinite(p.pos_err) && push!(z, x / p.pos_err, y / p.pos_err)
            end
        end
    end
    @printf("%-6s reliable %2d of %d | target identified in %2d | %3d catalogued, %3d unidentified | limiting mag median %.2f\n",
            arm, reliable, length(done), targets, nknown, unknown, isempty(lims) ? NaN : median(lims))
    isempty(errs) && continue
    @printf("       vs JPL: %d points, median |error| %.3f\", mean RA %+.3f\" Dec %+.3f\", error/pos_err std %.2f robust %.2f\n",
            length(dra), median(abs.(errs)), mean(dra), mean(ddec), isempty(z) ? NaN : std(z), isempty(z) ? NaN : robust(z))
end

println("\nper set: reliable raw/here/BANZAI, limiting mag here - BANZAI, target magnitude here - BANZAI")
zmag = Float64[]
for r in done
    h, b = r.arms["here"], r.arms["BANZAI"]
    dl = h.limiting_mag - b.limiting_mag
    th = [k for k in h.known if is_target(r.target, k.object)]
    tb = [k for k in b.known if is_target(r.target, k.object)]
    dm = Float64[]
    if !isempty(th) && !isempty(tb)
        for p in th[1].rows, q in tb[1].rows
            abs(p.epoch - q.epoch) < 1e-6 && isfinite(p.mag) && isfinite(q.mag) || continue
            push!(dm, p.mag - q.mag)
            push!(zmag, (p.mag - q.mag) / hypot(p.mag_err, q.mag_err))
        end
    end
    @printf("%-28s %s  %+.2f  %s\n", r.name, join([r.arms[a].reliable ? "y" : "n" for a in ARMS], "/"), dl,
            isempty(dm) ? "-" : @sprintf("%+.3f (%d points)", median(dm), length(dm)))
end
isempty(zmag) || @printf("target magnitude, here - BANZAI over the combined mag_err: %d points, std %.2f robust %.2f, median %+.2f\n",
                         length(zmag), std(zmag), robust(zmag), median(zmag))

finite_median(x) = (v = filter(isfinite, x); isempty(v) ? NaN : median(v))
imaged = filter(r -> haskey(r, :images) && r.images !== nothing, done)
isempty(imaged) || println("\nframe by frame, from the raw headers' WCS: frames refined (offset only), median stars, median scatter,",
                           "\nlargest difference from BANZAI's WCS (median, max), zero-point scatter (median), sky S/N over BANZAI's")
for r in imaged, arm in ARMS
    m = r.images[arm]
    off = filter(isfinite, m.off)
    @printf("%-28s %-6s %d (%d) of %d, %4.0f stars, %.2f\", off %s, zp scatter %.3f (%d stars), sky S/N %.2f\n",
            r.name, arm, count(m.refined), count(m.shift_only), length(m.refined), median(m.matches),
            finite_median(m.rms), isempty(off) ? "-" : @sprintf("%.2f\" / %.2f\"", median(off), maximum(off)),
            finite_median(m.zp_scatter), round(Int, median(m.zp_stars)),
            median(m.sky_snr ./ r.images["BANZAI"].sky_snr))
end
for arm in ARMS
    isempty(imaged) && break
    off = filter(isfinite, reduce(vcat, [r.images[arm].off for r in imaged]))
    sc = filter(isfinite, reduce(vcat, [r.images[arm].zp_scatter for r in imaged]))
    dsc = filter(isfinite, reduce(vcat, [r.images[arm].zp_scatter .- r.images["BANZAI"].zp_scatter for r in imaged]))
    doff = filter(isfinite, reduce(vcat, [r.images[arm].off .- r.images["BANZAI"].off for r in imaged]))
    @printf("%-6s refined %d of %d frames, off median %.2f\" max %.2f\"; zp scatter median %.3f, minus BANZAI's median %+.3f max %+.3f\n",
            arm, sum(count(r.images[arm].refined) for r in imaged), sum(length(r.images[arm].refined) for r in imaged),
            isempty(off) ? NaN : median(off), isempty(off) ? NaN : maximum(off), isempty(sc) ? NaN : median(sc),
            isempty(dsc) ? NaN : median(dsc), isempty(dsc) ? NaN : maximum(dsc))
    isempty(doff) || @printf("       off minus BANZAI's pixels': median %+.2f\" max %+.2f\"\n", median(doff), maximum(doff))
end
