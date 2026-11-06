"""
    aperture_flux(image, x, y; aperture_radius=5.0, annulus=(10.0, 15.0), mask=nothing) -> Float64

Background-subtracted flux inside a circle of `aperture_radius` pixels at
`(x, y)` (the same `x` = column, `y` = row sense as [`detect_sources`](@ref)'s
`xcen`/`ycen`), with the background taken as the median of the annulus
between `annulus[1]` and `annulus[2]` pixels — local, so a frame's large-scale
background gradients don't leak in. Pixels flagged in `mask` are left out of
both. `NaN` if the aperture leaves the image or its annulus has no valid pixel.
"""
function aperture_flux(image::AbstractMatrix{<:Real}, x::Real, y::Real;
                       aperture_radius::Real=5.0, annulus=(10.0, 15.0),
                       mask::Union{Nothing,AbstractMatrix{Bool}}=nothing)
    px = _aperture_pixels(image, x, y, aperture_radius, annulus, mask)
    px === nothing && return NaN
    return sum(px.inside) - median(px.ring) * length(px.inside)
end

# The valid pixel values inside the aperture and in the background annulus
# around (x, y); `nothing` if the annulus leaves the image or either is empty.
function _aperture_pixels(image, x, y, aperture_radius, annulus, mask)
    n1, n2 = size(image)
    r_out = annulus[2]
    i0, i1 = floor(Int, y - r_out), ceil(Int, y + r_out)
    j0, j1 = floor(Int, x - r_out), ceil(Int, x + r_out)
    (1 <= i0 && i1 <= n1 && 1 <= j0 && j1 <= n2) || return nothing
    inside = Float64[]; ring = Float64[]
    for j in j0:j1, i in i0:i1
        mask === nothing || !mask[i, j] || continue
        r = hypot(j - x, i - y)
        if r <= aperture_radius
            push!(inside, image[i, j])
        elseif annulus[1] <= r <= annulus[2]
            push!(ring, image[i, j])
        end
    end
    (isempty(inside) || isempty(ring)) && return nothing
    return (inside=inside, ring=ring)
end

"""
    photometric_zeropoint(image, wcs, stars; mask=nothing, aperture_radius=5.0,
                          mag_range=(15.0, 20.0), match_radius=1.5) -> (zeropoint, n_stars, scatter)

The magnitude zero point of a frame against reference stars measured in it:
each catalog star in `mag_range` (`stars` with `ra`, `dec`, `gmag`, e.g. from
[`gaia_reference_stars`](@ref)) is matched to the nearest detection within
`match_radius` pixels of its position through `wcs` — which therefore has to
be accurate, i.e. [`refine_wcs`](@ref)'s output for real IASC/PS1 headers —
and `zeropoint = median(gmag + 2.5 log10(aperture_flux))`, after rejecting
3σ outliers. `scatter` is the robust (MAD) spread of the individual stars'
values, a direct measure of how well the frame calibrates.

`mag_range` starts at 15 because brighter stars saturate a PS1 45 s exposure
(their cores are clipped flat — see [`fill_value_mask`](@ref)), and ends at 20
to keep each star's own photon noise small. The catalog band (Gaia G) is not
the frame's (PS1 w for IASC sets), so magnitudes on this zero point are
G-equivalent with a colour term left in — the same convention Astrometrica's
reports use ("21.2 G" in IASC reports), and why reports built on it carry
band `'G'`.

`NaN` zero point (and `n_stars = 0`) if no star matches.
"""
function photometric_zeropoint(image::AbstractMatrix{<:Real}, wcs::WCSTransform, stars;
                               mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                               aperture_radius::Real=5.0, mag_range=(15.0, 20.0),
                               match_radius::Real=1.5)
    detections = detect_sources(image; threshold=10.0, mask)
    n1, n2 = size(image)
    zps = Float64[]
    for s in stars
        (mag_range[1] <= s.gmag <= mag_range[2]) || continue
        x, y = world_to_pix(wcs, Float64[s.ra, s.dec])
        (1 <= x <= n2 && 1 <= y <= n1) || continue
        dist, i = isempty(detections) ? (Inf, 0) :
                  findmin(hypot(d.xcen - x, d.ycen - y) for d in detections)
        dist <= match_radius || continue
        flux = aperture_flux(image, detections[i].xcen, detections[i].ycen; aperture_radius, mask)
        flux > 0 && push!(zps, s.gmag + 2.5 * log10(flux))
    end
    isempty(zps) && return (zeropoint=NaN, n_stars=0, scatter=NaN)
    keep = trues(length(zps))
    while true
        m = median(zps[keep])
        σ = 1.4826 * median(abs.(zps[keep] .- m))
        new_keep = keep .& (abs.(zps .- m) .<= max(3σ, 1e-6))   # monotone: always terminates
        (new_keep == keep || count(new_keep) < 3) && break
        keep = new_keep
    end
    kept = zps[keep]
    m = median(kept)
    return (zeropoint=m, n_stars=length(kept), scatter=1.4826 * median(abs.(kept .- m)))
end

"""
    candidate_magnitudes(fits_paths, candidates; timestamp_key=nothing,
                         reference_stars=nothing, aperture_radius=5.0, wcs=nothing) -> Vector{Float64}

Calibrated (Gaia G-equivalent) magnitude of every row of `candidates` — an
[`astrometric_calibrate`](@ref) table from a [`run_pipeline`](@ref) call on
the same `fits_paths` — measured on each row's own *raw* frame, not on a
difference image: each frame's WCS is refined ([`refine_wcs`](@ref)) and its
zero point fitted ([`photometric_zeropoint`](@ref)) against Gaia DR3 stars
(`reference_stars`, or fetched once for the first frame's field — a live
network request), then each row's sky position is photometered at its pixel
in that frame. Working from `ra`/`dec` rather than the table's `x`/`y` keeps
this independent of which pixel grid the pipeline linked on. `wcs`, one
solution per frame, is used instead of refining each frame on its own, which
fails where a header is further off than [`refine_wcs`](@ref) searches:
[`analyze_set`](@ref) passes the frames' joint fit.

`NaN` where a row's aperture falls off the frame or measures no positive flux.
Each frame gets its own zero point, so a frame taken through thin cloud
doesn't make an object look fainter. `examples/iasc_demo.jl` compares the
result with Astrometrica's.
"""
candidate_magnitudes(fits_paths::AbstractVector{<:AbstractString}, candidates; kw...) =
    _measure_candidates(fits_paths, candidates; kw...).mags

# The candidates' magnitudes and the limit of the frames they are measured on
function _measure_candidates(fits_paths, candidates; timestamp_key::Union{Nothing,AbstractString}=nothing,
                             reference_stars=nothing, aperture_radius::Real=5.0, wcs=nothing)
    rows = collect(candidates)
    mags = fill(NaN, length(rows))
    limits = Float64[]
    for (k, path) in enumerate(fits_paths)
        idx = [i for (i, r) in enumerate(rows) if r.frame == k]
        isempty(idx) && continue
        image, header, epoch = FITS(path, "r") do f
            permutedims(Float64.(read(f[1]))), read_header(f[1], String),
            frame_epoch(read_header(f[1]); timestamp_key)
        end
        mask = fill_value_mask(image)
        frame_wcs = wcs === nothing ? load_wcs(header) : wcs[k]
        reference_stars === nothing &&
            (reference_stars = _field_reference_stars(image, frame_wcs, epoch))
        if wcs === nothing
            refined = refine_wcs(image, frame_wcs, reference_stars; mask)
            refined.refined && (frame_wcs = refined.wcs)
        end
        zp = photometric_zeropoint(image, frame_wcs, reference_stars; mask, aperture_radius)
        isnan(zp.zeropoint) && continue
        push!(limits, _limiting_magnitude(image, mask, zp.zeropoint, _source_radius(image, mask, aperture_radius)))
        for i in idx
            x, y = world_to_pix(frame_wcs, Float64[rows[i].ra, rows[i].dec])
            flux = aperture_flux(image, x, y; aperture_radius, mask)
            flux > 0 && (mags[i] = zp.zeropoint - 2.5 * log10(flux))
        end
    end
    finite = filter(isfinite, limits)
    return (mags=mags, limit=isempty(finite) ? NaN : median(finite))
end

# An aperture of one FWHM holds most of a point source's light and little sky,
# whatever the size of the pixels; `fallback` if no star shows the FWHM.
function _source_radius(image, mask, fallback)
    try
        return _measure_fwhm(image, mask)
    catch
        return fallback
    end
end

# The magnitude of a source 5σ above the sky's noise in an aperture of
# `aperture_radius`: what this frame can show, from its own noise and zero point.
function _limiting_magnitude(image, mask, zeropoint, aperture_radius)
    sky = image[.!mask]
    isempty(sky) && return NaN
    sigma = 1.4826 * median(abs.(sky .- median(sky)))
    sigma > 0 || return NaN
    return zeropoint - 2.5 * log10(5 * sigma * sqrt(pi) * aperture_radius)
end

"""
    below_limit(candidates, limiting_mag) -> Set{Int}

The tracklet ids in `candidates` whose median magnitude is fainter than
`limiting_mag`, what [`analyze_set`](@ref) returns: nothing a single exposure
can show, so most likely noise lined up by chance. They are worth looking at
last, not to be discarded: a stack of frames can hold what one doesn't.
"""
function below_limit(candidates, limiting_mag::Real)
    isfinite(limiting_mag) || return Set{Int}()
    rows = collect(candidates)
    faint = Set{Int}()
    for id in unique(r.id for r in rows)
        mags = filter(isfinite, [r.mag for r in rows if r.id == id])
        !isempty(mags) && median(mags) > limiting_mag && push!(faint, id)
    end
    return faint
end
