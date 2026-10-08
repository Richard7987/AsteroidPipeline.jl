"""
    inject_movers(fits_paths, out_dir; n=50, mag_range=(19.0, 22.5), speed_range=(5.0, 100.0),
                  rng=default_rng(), timestamp_key=nothing, exptime_key="EXPTIME",
                  reference_stars=nothing, margin=50) -> (paths, truth)

Write copies of a same-field sequence `fits_paths` into `out_dir` with `n`
synthetic asteroids added — known truth to measure the pipeline against,
at any brightness and speed, instead of the handful of real objects a few
sets happen to contain (6 across all the local IASC sets, which can't say
how completeness falls off with magnitude).

Each synthetic object moves in a straight line at constant speed, uniform
in `speed_range` (arcsec/hour) in a uniformly random direction, from a
uniformly random start inside the first frame (`margin` pixels from its
edges); its magnitude is uniform in `mag_range`. Each frame is calibrated
on its own — WCS refined ([`refine_wcs`](@ref)) and zero point fitted
([`photometric_zeropoint`](@ref)) against Gaia DR3 (`reference_stars`, or
fetched once — a live request) — so the object lands at its true sky
position and with its true brightness *in that frame*, transparency
changes included. It is drawn as a Gaussian of the frame's own measured
FWHM, smeared along its motion during the exposure (`exptime_key`; fast
objects really do trail, and whether the detector finds trails is part of
what this measures), with Poisson noise from its own counts (header `GAIN`,
default 1). Detector gaps are not avoided: an object crossing one is lost
there, as a real one would be.

Returns the written paths (same order) and a `truth` table, one row per
object per frame: `object`, `frame`, `ra`, `dec` (degrees, mid-exposure),
`epoch` (Julian Date, mid-exposure — the same convention
[`run_pipeline`](@ref) now uses), `mag`, `speed` (arcsec/hour), `x`, `y`
(pixels in that frame). Score a run with [`injection_recovery`](@ref).
"""
function inject_movers(fits_paths::AbstractVector{<:AbstractString}, out_dir::AbstractString;
                       n::Integer=50, mag_range=(19.0, 22.5), speed_range=(5.0, 100.0),
                       rng::AbstractRNG=default_rng(), timestamp_key::Union{Nothing,AbstractString}=nothing,
                       exptime_key::AbstractString="EXPTIME", reference_stars=nothing, margin::Integer=50)
    frames = map(fits_paths) do path
        FITS(path, "r") do f
            h = read_header(f[1])
            (image=permutedims(Float64.(read(f[1]))), header=h, header_str=read_header(f[1], String),
             epoch=frame_epoch(h; timestamp_key, exptime_key),
             exptime=haskey(h, exptime_key) ? Float64(h[exptime_key]) : 0.0,
             gain=haskey(h, "GAIN") ? Float64(h["GAIN"]) : 1.0)
        end
    end
    epochs = [fr.epoch for fr in frames]

    calib = map(enumerate(frames)) do (k, fr)
        mask = fill_value_mask(fr.image)
        wcs = load_wcs(fr.header_str)
        reference_stars === nothing &&
            (reference_stars = _field_reference_stars(fr.image, wcs, epochs[1]))
        refined = refine_wcs(fr.image, wcs, reference_stars; mask)
        refined.refined && (wcs = refined.wcs)
        zp = photometric_zeropoint(fr.image, wcs, reference_stars; mask)
        isnan(zp.zeropoint) && error("no photometric zero point for $(fits_paths[k])")
        (wcs=wcs, zeropoint=zp.zeropoint, fwhm=_measure_fwhm(fr.image, mask))
    end

    n1, n2 = size(frames[1].image)
    object = Int[]; frame = Int[]; ra = Float64[]; dec = Float64[]; epoch = Float64[]
    mag = Float64[]; speed = Float64[]; xs = Float64[]; ys = Float64[]
    images = [copy(fr.image) for fr in frames]
    for obj in 1:n
        x0, y0 = margin + rand(rng) * (n2 - 2margin), margin + rand(rng) * (n1 - 2margin)
        ra0, dec0 = pix_to_world(calib[1].wcs, [x0, y0])
        θ = 2π * rand(rng)
        v = speed_range[1] + rand(rng) * (speed_range[2] - speed_range[1])   # arcsec/h
        m = mag_range[1] + rand(rng) * (mag_range[2] - mag_range[1])
        # sky position at time t (JD), moving v arcsec/h at position angle θ
        at(t) = let dt_h = (t - epochs[1]) * 24
            (ra0 + v * dt_h * sin(θ) / 3600 / cosd(dec0), dec0 + v * dt_h * cos(θ) / 3600)
        end
        for (k, fr) in enumerate(frames)
            c = calib[k]
            half = fr.exptime / 2 / 86400
            sky_mid = at(epochs[k])
            px = world_to_pix(c.wcs, collect(sky_mid))
            p_start = world_to_pix(c.wcs, collect(at(epochs[k] - half)))
            p_end = world_to_pix(c.wcs, collect(at(epochs[k] + half)))
            counts = 10^(-0.4 * (m - c.zeropoint))
            _add_trail!(images[k], p_start, p_end, counts, c.fwhm, fr.gain, rng)
            push!(object, obj); push!(frame, k); push!(ra, sky_mid[1]); push!(dec, sky_mid[2])
            push!(epoch, epochs[k]); push!(mag, m); push!(speed, v); push!(xs, px[1]); push!(ys, px[2])
        end
    end

    mkpath(out_dir)
    paths = map(enumerate(fits_paths)) do (k, path)
        out = joinpath(out_dir, basename(path))
        FITS(out, "w") do f
            write(f, permutedims(images[k]); header=frames[k].header)
        end
        out
    end
    return paths, Table(; object, frame, ra, dec, epoch, mag, speed, x=xs, y=ys)
end

# Add `counts` spread as a Gaussian (FWHM `fwhm` px) smeared uniformly
# along the segment from `p_start` to `p_end` (x, y pixels) — the trail an
# object moving during the exposure leaves — plus Poisson noise on what was
# added (Gaussian approximation, variance = signal / gain).
function _add_trail!(image, p_start, p_end, counts::Real, fwhm::Real, gain::Real, rng::AbstractRNG)
    σ = fwhm / (2 * sqrt(2 * log(2)))
    len = hypot(p_end[1] - p_start[1], p_end[2] - p_start[2])
    nsub = max(1, ceil(Int, len / 0.25))              # sub-steps every quarter pixel
    r = ceil(Int, 5σ)
    n1, n2 = size(image)
    xlo = max(1, floor(Int, min(p_start[1], p_end[1])) - r); xhi = min(n2, ceil(Int, max(p_start[1], p_end[1])) + r)
    ylo = max(1, floor(Int, min(p_start[2], p_end[2])) - r); yhi = min(n1, ceil(Int, max(p_start[2], p_end[2])) + r)
    (xlo <= xhi && ylo <= yhi) || return image
    model = zeros(yhi - ylo + 1, xhi - xlo + 1)
    norm = counts / nsub / (2π * σ^2)
    for s in 1:nsub
        f = nsub == 1 ? 0.5 : (s - 1) / (nsub - 1)
        cx = p_start[1] + f * (p_end[1] - p_start[1]); cy = p_start[2] + f * (p_end[2] - p_start[2])
        for (jj, j) in enumerate(xlo:xhi), (ii, i) in enumerate(ylo:yhi)
            model[ii, jj] += norm * exp(-((j - cx)^2 + (i - cy)^2) / (2σ^2))
        end
    end
    for (jj, j) in enumerate(xlo:xhi), (ii, i) in enumerate(ylo:yhi)
        mv = model[ii, jj]
        image[i, j] += mv + sqrt(max(mv, 0.0) / gain) * randn(rng)
    end
    return image
end

# Median FWHM (pixels) of bright, unsaturated, isolated stars in `image`,
# from second moments inside a 7x7 box around each detection.
function _measure_fwhm(image::AbstractMatrix, mask)
    dets = detect_sources(image; threshold=20.0, mask, min_sharpness=0.3)
    n1, n2 = size(image)
    bg = median(image[.!mask])
    σs = Float64[]
    for d in dets
        x, y = round(Int, d.xcen), round(Int, d.ycen)
        (8 < x <= n2 - 8 && 8 < y <= n1 - 8) || continue
        any(@view mask[y-4:y+4, x-4:x+4]) && continue           # touches a gap or a saturated core
        count(e -> hypot(e.xcen - d.xcen, e.ycen - d.ycen) < 12, dets) == 1 || continue
        w = max.(image[y-3:y+3, x-3:x+3] .- bg, 0.0)
        total = sum(w)
        total > 0 || continue
        dx = (-3:3)'; dy = -3:3
        mx = sum(w .* dx) / total; my = sum(w .* dy) / total
        varx = sum(w .* (dx .- mx) .^ 2) / total; vary = sum(w .* (dy .- my) .^ 2) / total
        push!(σs, sqrt((varx + vary) / 2))
    end
    isempty(σs) && error("no star bright and isolated enough to measure the PSF")
    return 2 * sqrt(2 * log(2)) * median(σs)
end

"""
    injection_recovery(candidates, truth; radius=1.5, min_frames=3) -> Table

Score a [`run_pipeline`](@ref) run on [`inject_movers`](@ref) output: a
synthetic object is recovered if one tracklet in `candidates` lies within
`radius` arcsec of its true position in at least `min_frames` frames —
the same frame-by-frame standard the demos apply to real known objects.
Returns one row per object: `object`, `mag`, `speed`, `recovered`,
`tracklet` (the matching tracklet's `id`, or 0).
"""
function injection_recovery(candidates, truth; radius::Real=1.5, min_frames::Integer=3)
    rows = collect(candidates)
    objects = sort(unique(truth.object))
    out_mag = Float64[]; out_speed = Float64[]; recovered = Bool[]; tracklet = Int[]
    for obj in objects
        t = filter(r -> r.object == obj, collect(truth))
        pos = Dict(r.frame => (r.ra, r.dec) for r in t)
        hits = Dict{Int,Int}()
        for r in rows
            haskey(pos, r.frame) || continue
            p = pos[r.frame]
            _angular_distance_arcsec(r.ra, r.dec, p[1], p[2]) <= radius &&
                (hits[r.id] = get(hits, r.id, 0) + 1)
        end
        best = isempty(hits) ? (0, 0) : findmax(hits)
        push!(out_mag, t[1].mag); push!(out_speed, t[1].speed)
        push!(recovered, best[1] >= min_frames); push!(tracklet, best[1] >= min_frames ? best[2] : 0)
    end
    return Table(object=objects, mag=out_mag, speed=out_speed, recovered=recovered, tracklet=tracklet)
end
