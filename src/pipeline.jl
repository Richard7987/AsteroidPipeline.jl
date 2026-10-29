"""
    run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                 timestamp_key::Union{Nothing,AbstractString}=nothing,
                 threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                 aperture_radius::Real=3.0, max_speed::Real=Inf,
                 match_radius::Real=2.0, min_frames::Union{Nothing,Integer}=nothing,
                 reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                 quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing,
                 photometric_outlier_threshold::Real=0.2,
                 mask_fill_values::Bool=false, difference_stack::Bool=false,
                 refine_astrometry::Bool=false, reference_stars=nothing,
                 min_speed::Real=0.0, max_flux_ratio::Real=Inf, max_residual::Real=Inf,
                 min_sharpness::Real=-Inf, exptime_key::Union{Nothing,AbstractString}="EXPTIME",
                 subpixel_alignment::Bool=false, max_background_excess::Real=Inf,
                 detector_artifact_radius::Real=0.0)

Run the local stages of the pipeline — detection, linking, and
astrometric calibration — on a time-ordered sequence of FITS frames from
the same field.

For each frame, the primary header supplies the WCS (via [`load_wcs`](@ref))
and the epoch, a UTC Julian Date at mid-exposure, from [`frame_epoch`](@ref)
with `timestamp_key` and `exptime_key`. Detections are linked across frames
with [`link_candidates`](@ref) and calibrated to sky coordinates with
[`astrometric_calibrate`](@ref).

If `reference` is `nothing` (the default), [`detect_sources`](@ref) runs
directly on each frame. Passing a deep reference stack — an
`(image, sigma, mask, psf, wcs)` named tuple, with `image`/`sigma`/`mask`
from [`build_reference`](@ref), `wcs` the `target_wcs` it was built on and
`psf` from `estimate_psf(image)` — switches to ZOGY difference imaging
([`zogy_subtract`](@ref)): each frame is reprojected onto `reference`'s
grid, subtracted, and detected on the resulting `S_corr`, which reaches
below the single-frame noise floor. `threshold` is then a significance in
sigma, and every frame is calibrated with `reference.wcs`. Each frame's PSF
is estimated after reprojection, on the pixels `zogy_subtract` differences,
and bright stars in the frame and the reference are passed as
`n_sources`/`r_sources` for its astrometric-noise term.

On the `reference` path, a frame whose `S_corr` standard deviation (over
its valid region) exceeds `quality_max_std` is treated as having found
nothing, with a warning. The plain `std` is deliberate: a MAD-based spread
misses a frame whose excess is a few hundred bright residuals. The default
`1.5` was set on ZTF data (good frames ~1.1-1.2); tune it per survey.
With `min_frames=nothing` (the default), linking requires every frame that
was not gated; pass `min_frames` to set it explicitly.

On the raw path, a frame whose [`photometric_scale`](@ref) factor deviates
from the median by more than `photometric_outlier_threshold` (20%) gets a
warning, never an exclusion. It catches a uniform flux-scale change
(transparency, focus); a thin passing cloud can leave it unchanged.

A frame whose header has no WCS is solved from its pointing
([`solve_wcs`](@ref), against Gaia DR3 stars fetched once per field); if
that fails, with `plate_solve_api_key` (a nova.astrometry.net API key) it is
solved by [`plate_solve`](@ref) instead.

These options, all off by default and all on the raw path only, are what
IASC/Pan-STARRS1 sets need:

- `mask_fill_values`: mask constant-valued plateaus (the gaps between
  detector cells; see [`fill_value_mask`](@ref)) before detection.
- `refine_astrometry`: replace each frame's header WCS with one fitted to
  reference stars in all frames together ([`refine_wcs_set`](@ref)), against
  `reference_stars` if given, else Gaia DR3 stars covering the first frame
  ([`gaia_reference_stars`](@ref), a network request). A frame whose
  detections far outnumber the reference stars gets a warning
  ([`_check_frame_quality`](@ref)).
- `difference_stack`: detect on each frame minus the median of the whole
  sequence ([`stack_difference`](@ref)): static sources cancel, so
  `threshold` can go much lower. Needs at least 3 frames. A fast object's
  trail is merged into one detection, satellite trails and saturated-star
  residuals dropped (see [`_merge_clusters`](@ref)).
- `subpixel_alignment` (difference path): align frames by fractional
  offsets, with interpolation (see [`stack_difference`](@ref)'s `subpixel`).
- `min_speed` (pixels per day), `max_flux_ratio` and `max_residual`
  (pixels) filter tracklets after linking, per IASC's "true signature"
  tests (see [`_filter_tracklets`](@ref)); `max_flux_ratio=2.5` is IASC's
  1-magnitude limit. With `refine_astrometry`, `max_residual` is raised to
  three times the frames' typical astrometric residual (at most 2 px): a
  straight line can't be checked more finely than the positions are measured.
- `detector_artifact_radius` (pixels; difference path): drop a detection
  that recurs at the same *detector* pixel in two or more other frames — a
  hot pixel or bad column (see [`_drop_detector_fixed`](@ref)).
- `max_background_excess` (σ; difference path): drop a detection whose
  surroundings in the raw frame sit more than this far above the sky — a
  bright star's halo, ghost or spike (see [`_drop_on_bright_background`](@ref)).
- `min_sharpness`: drop single-pixel spikes — cosmic rays, hot pixels —
  at detection (see [`detect_sources`](@ref)); `0.3` on IASC data.

`fits_paths` must be in time order. `psf_threshold` and `psf_min_separation`
go to `estimate_psf` for each frame's PSF (`reference` path only). All other
keyword arguments are forwarded to `detect_sources` and `link_candidates`.

Returns the candidate table from `astrometric_calibrate` (columns `id`,
`frame`, `x`, `y`, `ra`, `dec`, `epoch`), ready for
[`crossmatch_catalog`](@ref), which is left to the caller since it needs
the network and a choice of catalog.
"""
function run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                       timestamp_key::Union{Nothing,AbstractString}=nothing,
                       threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                       aperture_radius::Real=3.0, max_speed::Real=Inf,
                       match_radius::Real=2.0, min_frames::Union{Nothing,Integer}=nothing,
                       reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                       quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing,
                       photometric_outlier_threshold::Real=0.2,
                       mask_fill_values::Bool=false, difference_stack::Bool=false,
                       refine_astrometry::Bool=false, reference_stars=nothing,
                       min_speed::Real=0.0, max_flux_ratio::Real=Inf, max_residual::Real=Inf,
                       min_sharpness::Real=-Inf, exptime_key::Union{Nothing,AbstractString}="EXPTIME",
                       subpixel_alignment::Bool=false, max_background_excess::Real=Inf,
                       detector_artifact_radius::Real=0.0)
    detections_per_frame, wcs_per_frame, timestamps, n_gated, astrometric_rms_px = _detect_all_frames(
        fits_paths; timestamp_key, threshold, box_size, aperture_radius,
        reference, psf_threshold, psf_min_separation, quality_max_std, plate_solve_api_key,
        photometric_outlier_threshold, mask_fill_values, difference_stack,
        refine_astrometry, reference_stars, min_sharpness, exptime_key, subpixel_alignment,
        max_background_excess, detector_artifact_radius)

    effective_min_frames = min_frames === nothing ? length(fits_paths) - n_gated : min_frames
    tracklets = link_candidates(detections_per_frame, timestamps;
                                 max_speed, match_radius, min_frames=effective_min_frames)
    # at most 2 px: chance alignments sit a median 2.6 px off a line
    isfinite(astrometric_rms_px) && (max_residual = max(max_residual, min(3 * astrometric_rms_px, 2.0)))
    tracklets = _filter_tracklets(tracklets, detections_per_frame, timestamps; min_speed, max_flux_ratio, max_residual,
                                  min_frames=effective_min_frames)
    return astrometric_calibrate(tracklets, wcs_per_frame, timestamps)
end

"""
    _detect_all_frames(fits_paths; timestamp_key, threshold, box_size, aperture_radius,
                        reference, psf_threshold, psf_min_separation,
                        quality_max_std, plate_solve_api_key)
        -> (detections_per_frame, wcs_per_frame, timestamps, n_gated, astrometric_rms_px)

The per-frame detection stage shared by [`run_pipeline`](@ref) and
[`search_field`](@ref), so both use one pass over `fits_paths`. See
`run_pipeline` for the keywords. On the raw path each frame's `GAIN` header
keyword (default `1.0`) goes to [`detect_sources`](@ref), so `flux_err`
includes shot noise, which [`find_variable_sources`](@ref)'s chi-squared
test needs.

`n_gated` is how many frames the `quality_max_std` gate excluded (always
`0` on the raw path): both callers lower their `min_frames` default by it.
`astrometric_rms_px` is the median refined-astrometry residual in pixels
(`NaN` without `refine_astrometry`).
"""
function _detect_all_frames(fits_paths::AbstractVector{<:AbstractString};
                             timestamp_key::Union{Nothing,AbstractString}=nothing,
                             threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                             aperture_radius::Real=3.0,
                             reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                             quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing,
                             photometric_outlier_threshold::Real=0.2,
                             mask_fill_values::Bool=false, difference_stack::Bool=false,
                             refine_astrometry::Bool=false, reference_stars=nothing,
                             min_sharpness::Real=-Inf,
                             exptime_key::Union{Nothing,AbstractString}="EXPTIME",
                             subpixel_alignment::Bool=false, max_background_excess::Real=Inf,
                             detector_artifact_radius::Real=0.0)
    reference === nothing || !(mask_fill_values || difference_stack || refine_astrometry) ||
        throw(ArgumentError("mask_fill_values, difference_stack and refine_astrometry apply only " *
                            "without `reference` (the ZOGY path has its own registration and masks)"))
    detections_per_frame = typeof(_empty_detections())[]
    wcs_per_frame = WCSTransform[]
    timestamps = Float64[]
    n_gated = 0
    rms_px = Float64[]       # refined astrometry's residual per frame, in pixels
    # difference_stack needs every frame before it can detect in any of them
    stacked_images = Matrix{Float64}[]
    stacked_masks = Union{Nothing,BitMatrix}[]
    refine_images = Matrix{Float64}[]
    refine_masks = Union{Nothing,BitMatrix}[]

    # The reference's stars, for zogy_subtract's astrometric-noise term.
    r_sources = reference === nothing ? nothing :
                detect_sources(permutedims(reference.image); threshold=psf_threshold)

    solve_cache = Ref{Any}(nothing)
    for path in fits_paths
        FITS(path, "r") do f
            hdu = f[1]

            # FITSIO reads (x, y); detect_sources wants (y, x). The reference and
            # ZOGY code work in (x, y) and are permuted only for detection.
            raw = Float64.(read(hdu))
            epoch = frame_epoch(read_header(hdu); timestamp_key, exptime_key)
            frame_wcs = try
                first(_header_or_solved_wcs(permutedims(raw), read_header(hdu, String), read_header(hdu),
                                            epoch, solve_cache))
            catch e
                plate_solve_api_key === nothing && rethrow(e)
                plate_solve(path; api_key=plate_solve_api_key)
            end
            gain = haskey(read_header(hdu), "GAIN") ? read_key(hdu, "GAIN")[1] : 1.0

            push!(timestamps, epoch)

            if reference === nothing
                image = permutedims(raw)
                mask = mask_fill_values ? fill_value_mask(image) : nothing
                if refine_astrometry
                    reference_stars === nothing &&
                        (reference_stars = _field_reference_stars(image, frame_wcs, timestamps[end]))
                    # refitted below, all frames together
                    push!(refine_images, image)
                    push!(refine_masks, mask)
                end
                push!(wcs_per_frame, frame_wcs)
                if difference_stack
                    push!(stacked_images, image)
                    push!(stacked_masks, mask)
                else
                    push!(detections_per_frame,
                          detect_sources(image; threshold, box_size, aperture_radius, gain, mask, min_sharpness))
                end
            else
                _, sigma_n = estimate_background(permutedims(raw);
                                                  location=SourceExtractorBackground(), rms=MADStdRMS())

                resampled, frame_mask = Reproject.reproject((raw, frame_wcs), reference.wcs;
                                                              shape_out=size(reference.image))
                # One NaN from the reprojection's border turns all of S_corr NaN;
                # filled with the valid level, as 0 would bias the background.
                resampled[.!frame_mask] .= median(resampled[frame_mask])

                # Measured on the resampled pixels, which are what zogy_subtract differences.
                psf_n = estimate_psf(resampled; threshold=psf_threshold, min_separation=psf_min_separation)
                n_sources = detect_sources(permutedims(resampled); threshold=psf_threshold)

                _, s_corr = zogy_subtract(resampled, reference.image;
                                           psf_n=psf_n, psf_r=reference.psf,
                                           sigma_n=sigma_n, sigma_r=reference.sigma, gain_n=gain,
                                           n_sources=n_sources, r_sources=r_sources)
                valid = frame_mask .& reference.mask
                # Plain std, not MAD — see the docstring for why.
                frame_std = std(s_corr[valid])
                if frame_std > quality_max_std
                    @warn "skipping frame: S_corr std exceeds quality_max_std" path frame_std quality_max_std
                    image = permutedims(zeros(size(s_corr)))
                    n_gated += 1
                else
                    image = permutedims(s_corr .* valid)
                end
                push!(wcs_per_frame, reference.wcs)
                # No gain: S_corr is a significance map, not counts.
                push!(detections_per_frame, detect_sources(image; threshold, box_size, aperture_radius))
            end
        end
    end

    if !isempty(refine_images)
        refined = refine_wcs_set(refine_images, wcs_per_frame, reference_stars; masks=refine_masks)
        for (path, image, r) in zip(fits_paths, refine_images, refined)
            if r.refined
                push!(rms_px, r.rms_arcsec / _pixel_scale_arcsec(r.wcs, size(image)))
            elseif r.shift_only
                @warn "astrometric refinement: too few stars; only the header WCS's offset corrected" path n_matches=r.n_matches
            else
                @warn "astrometric refinement failed; keeping the header WCS" path n_matches=r.n_matches
            end
            _check_frame_quality(path, r.n_detections, length(reference_stars))
        end
        wcs_per_frame = [r.wcs for r in refined]
    end

    if difference_stack
        masks = all(isnothing, stacked_masks) ? nothing : stacked_masks
        differences, diff_masks = stack_difference(stacked_images, wcs_per_frame; masks, subpixel=subpixel_alignment)
        # No gain: a mover's shot noise is negligible at the faint end this mode is for.
        per_frame = [_merge_clusters(detect_sources(diff; threshold, box_size, aperture_radius, mask=m, min_sharpness),
                                     _CLUSTER_RADIUS, _MAX_TRAIL_EXTENT)
                     for (diff, m) in zip(differences, diff_masks)]
        detector_artifact_radius > 0 && (per_frame = _drop_detector_fixed(per_frame, detector_artifact_radius, timestamps))
        if isfinite(max_background_excess)
            per_frame = [_drop_on_bright_background(d, img, m, max_background_excess)
                         for (d, img, m) in zip(per_frame, stacked_images, stacked_masks)]
        end
        append!(detections_per_frame, per_frame)
    end

    if reference === nothing && length(fits_paths) >= 2
        detections_per_frame = _to_common_grid(detections_per_frame, wcs_per_frame)
        wcs_per_frame = fill(wcs_per_frame[1], length(wcs_per_frame))
    end

    # Skipped for difference_stack: its detections are what *changed*
    # between frames, so there's no ensemble of stars to compare.
    if reference === nothing && !difference_stack && length(fits_paths) >= 2
        scales = photometric_scale(detections_per_frame)
        reference_scale = median(scales)
        for (k, path) in enumerate(fits_paths)
            deviation = abs(scales[k] - reference_scale) / reference_scale
            if deviation > photometric_outlier_threshold
                @warn "frame's photometric scale deviates from the field's median" path scale=scales[k] reference_scale deviation photometric_outlier_threshold
            end
        end
    end

    return detections_per_frame, wcs_per_frame, timestamps, n_gated, isempty(rms_px) ? NaN : median(rms_px)
end

"""
    search_field(fits_paths; <run_pipeline's keywords through photometric_outlier_threshold>,
                 variability_position_tolerance::Real=2.0,
                 variability_min_frames::Union{Nothing,Integer}=nothing,
                 variability_chi2_threshold::Real=10.0,
                 variability_normalize::Bool=(reference === nothing),
                 variability_max_relative_error::Real=0.10,
                 variability_systematic_error_fraction::Real=0.03)
        -> (movers=<table>, variables=<table>)

Run [`run_pipeline`](@ref)'s mover search and [`find_variable_sources`](@ref)'s
variable search on one shared detection pass over `fits_paths`
([`_detect_all_frames`](@ref)).

Keywords through `photometric_outlier_threshold` are `run_pipeline`'s; the
`variability_*` ones go to `find_variable_sources` without the prefix.
`variability_min_frames=nothing` means every non-gated frame, as for
`min_frames`. `variability_normalize` is off on the ZOGY path, where
`S_corr` isn't on a flux scale.

Returns `(movers=..., variables=...)`, each an `astrometric_calibrate`
table (columns `id`, `frame`, `x`, `y`, `ra`, `dec`, `epoch`); `variables`
is ready for `crossmatch_catalog(...; :vsx)`.
"""
function search_field(fits_paths::AbstractVector{<:AbstractString};
                       timestamp_key::Union{Nothing,AbstractString}=nothing,
                       threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                       aperture_radius::Real=3.0, max_speed::Real=Inf,
                       match_radius::Real=2.0, min_frames::Union{Nothing,Integer}=nothing,
                       reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                       quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing,
                       photometric_outlier_threshold::Real=0.2,
                       variability_position_tolerance::Real=2.0,
                       variability_min_frames::Union{Nothing,Integer}=nothing,
                       variability_chi2_threshold::Real=10.0,
                       variability_normalize::Bool=(reference === nothing),
                       variability_max_relative_error::Real=0.10,
                       variability_systematic_error_fraction::Real=0.03)
    detections_per_frame, wcs_per_frame, timestamps, n_gated, _ = _detect_all_frames(
        fits_paths; timestamp_key, threshold, box_size, aperture_radius,
        reference, psf_threshold, psf_min_separation, quality_max_std, plate_solve_api_key,
        photometric_outlier_threshold)

    effective_min_frames = min_frames === nothing ? length(fits_paths) - n_gated : min_frames
    tracklets = link_candidates(detections_per_frame, timestamps;
                                 max_speed, match_radius, min_frames=effective_min_frames)
    movers = astrometric_calibrate(tracklets, wcs_per_frame, timestamps)

    effective_variability_min_frames = variability_min_frames === nothing ?
                                        length(fits_paths) - n_gated : variability_min_frames
    variable_groups = find_variable_sources(
        detections_per_frame, timestamps;
        position_tolerance=variability_position_tolerance,
        min_frames=effective_variability_min_frames,
        chi2_threshold=variability_chi2_threshold,
        normalize=variability_normalize,
        max_relative_error=variability_max_relative_error,
        systematic_error_fraction=variability_systematic_error_fraction)
    variables = astrometric_calibrate(variable_groups, wcs_per_frame, timestamps)

    return (movers=movers, variables=variables)
end

"""
    _to_common_grid(detections_per_frame, wcs_per_frame) -> Vector

Every frame's detections re-expressed in the *first* frame's pixel grid
(through each frame's own WCS, then the first frame's inverse), so the
whole raw-path sequence shares one coordinate system — afterwards the
first frame's WCS applies to every frame's detections.

[`link_candidates`](@ref) and [`find_variable_sources`](@ref) compare
positions in pixels on one shared grid; without it, each frame's pointing
offset bends a straight track. Only the frames' *relative* WCS accuracy
matters: an error common to all of them cancels.

`x`/`y` become the transformed position rounded to whole pixels; `xcen`/`ycen`
keep it at full precision.
"""
function _to_common_grid(detections_per_frame, wcs_per_frame)
    w1 = wcs_per_frame[1]
    return map(zip(detections_per_frame, wcs_per_frame)) do (dets, w)
        (w === w1 || isempty(dets)) && return dets
        pix = Matrix{Float64}(undef, 2, length(dets))
        for (i, d) in enumerate(dets)
            pix[1, i], pix[2, i] = _position(d)
        end
        common = world_to_pix(w1, pix_to_world(w, pix))
        Table(dets; x=round.(Int, common[1, :]), y=round.(Int, common[2, :]),
              xcen=common[1, :], ycen=common[2, :])
    end
end

# Grouping radius and widest trail (pixels) for `_merge_clusters` on the
# difference_stack path. 60 px is ~15" at PS1's scale: an object moving
# ~1200"/h trails that far in a 45 s IASC exposure.
const _CLUSTER_RADIUS = 15.0
const _MAX_TRAIL_EXTENT = 60.0

"""
    _drop_detector_fixed(per_frame, radius, timestamps) -> Vector

Each frame's detections (own pixel grid) minus those with a detection at
the same pixel, within `radius`, in two or more other frames. A moving
object can't sit on one detector pixel for three exposures; a hot pixel or
bad column does, and only seems to move on the sky because the telescope
dithers between exposures — so the difference images don't cancel it.
The frames it recurs in must also cover half the sequence's time span: in
exposures seconds apart, a slow mover stays within `radius` of itself in
the frames next to it.
"""
function _drop_detector_fixed(per_frame, radius::Real, timestamps)
    xy = [[_position(d) for d in dets] for dets in per_frame]
    span = maximum(timestamps) - minimum(timestamps)
    return map(enumerate(per_frame)) do (k, dets)
        keep = map(xy[k]) do p
            others = [j for j in eachindex(xy) if j != k && any(q -> hypot(q[1] - p[1], q[2] - p[2]) <= radius, xy[j])]
            t = timestamps[[k; others]]
            length(others) < 2 || maximum(t) - minimum(t) < span / 2
        end
        dets[keep]
    end
end

"""
    _drop_on_bright_background(detections, image, mask, max_excess) -> table

`detections` minus those whose surroundings in the raw `image` — the median
of an 8-20 px annulus, in units of the frame's own sky noise — sit more
than `max_excess` above the sky. Real moving objects are found on clean
sky; what a difference image leaves near a very bright star is its halo,
an optical ghost, or a diffraction spike that turns between exposures, all
of which can move just plausibly enough to link. Real objects on IASC sets
sit within ±0.12σ; ghosts and spikes reached 1.7-2.4σ.
"""
function _drop_on_bright_background(detections, image, mask, max_excess::Real)
    valid = mask === nothing ? vec(image) : image[.!mask]
    sky = median(valid)
    σ = 1.4826 * median(abs.(valid .- sky))
    σ > 0 || return detections
    n1, n2 = size(image)
    keep = map(detections) do d
        x, y = round(Int, d.xcen), round(Int, d.ycen)
        ring = Float64[]
        for i in max(1, y - 20):min(n1, y + 20), j in max(1, x - 20):min(n2, x + 20)
            (mask === nothing || !mask[i, j]) || continue
            8 <= hypot(i - y, j - x) <= 20 && push!(ring, image[i, j])
        end
        isempty(ring) || (median(ring) - sky) / σ <= max_excess
    end
    return detections[keep]
end

"""
    _check_frame_quality(path, n_detections, n_reference_stars) -> Float64

Warn when a frame holds far more sources at [`refine_wcs`](@ref)'s
threshold (10σ) than reference stars in its field — the signature of
defective detector cells (speckled noise, gradients) rather than a crowded
field: good sets show a few hundred to ~900, a defective one 20,000 or
more. Returns the ratio.
"""
function _check_frame_quality(path, n_detections::Integer, n_reference_stars::Integer)
    ratio = n_detections / max(n_reference_stars, 1)
    ratio > _MAX_DETECTIONS_PER_REFERENCE_STAR &&
        @warn "frame has far more sources than reference stars: likely detector defects; treat this set's candidates as unreliable" path n_detections n_reference_stars
    return ratio
end

const _MAX_DETECTIONS_PER_REFERENCE_STAR = 30

"""
    set_quality(fits_paths; reference_stars=nothing) -> (reliable, ratios)

Whether a same-field set is clean enough to trust its candidates: each
frame's count of 10σ sources (fill-value gaps masked) over the number of
reference stars in its field — Gaia DR3 stars from
[`gaia_reference_stars`](@ref) if `reference_stars` isn't given (a live
request). `reliable` is `false` if any frame exceeds
`_MAX_DETECTIONS_PER_REFERENCE_STAR` (see [`_check_frame_quality`](@ref),
whose warning [`run_pipeline`](@ref) emits on the same test). A cheap
check to run before reporting a set's candidates.
"""
function set_quality(fits_paths::AbstractVector{<:AbstractString}; reference_stars=nothing)
    ratios = map(fits_paths) do path
        image, header, hdr = FITS(path, "r") do f
            permutedims(Float64.(read(f[1]))), read_header(f[1], String), read_header(f[1])
        end
        mask = fill_value_mask(image)
        reference_stars === nothing &&
            (reference_stars = _field_reference_stars(image, load_wcs(header), frame_epoch(hdr)))
        length(detect_sources(image; threshold=10.0, mask)) / max(length(reference_stars), 1)
    end
    return (reliable=all(<=(_MAX_DETECTIONS_PER_REFERENCE_STAR), ratios), ratios=ratios)
end
