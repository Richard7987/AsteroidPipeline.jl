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
                 min_speed::Real=0.0, max_flux_ratio::Real=Inf, min_sharpness::Real=-Inf,
                 exptime_key::Union{Nothing,AbstractString}="EXPTIME")

Run the local stages of the pipeline — detection, linking, and
astrometric calibration — on a time-ordered sequence of FITS frames from
the same field.

For each frame: the primary image is read; the primary header supplies the
WCS solution (via [`load_wcs`](@ref)) and the observation epoch, taken
as a UTC Julian Date at mid-exposure, by [`frame_epoch`](@ref) with this
call's `timestamp_key` and `exptime_key` — see it for why the default
reads PS1's `DATE-OBS` rather than its `MJD-OBS` (TAI, 37 s off), and why
half the exposure time is added (both surveys stamp the start; on a real
IASC set that left every epoch 22.5 s early against Astrometrica's report
of the same object). Detections are then linked across frames with
[`link_candidates`](@ref) and calibrated to sky coordinates with
[`astrometric_calibrate`](@ref).

If `reference` is `nothing` (the default), [`detect_sources`](@ref) runs
directly on each science frame — limited to sources well above the
frame's own background noise. Passing a deep reference stack (a
`(image, sigma, mask, psf, wcs)` named tuple — `image`/`sigma`/`mask`/`wcs`
from [`build_reference`](@ref) (`wcs` is the `target_wcs` it was built
on), `psf` from `estimate_psf(image)` on that same reference image)
switches to ZOGY difference imaging ([`zogy_subtract`](@ref)): each frame
is reprojected onto `reference`'s pixel grid (frames dither by several
pixels even within one night — see `src/reference.jl` — so this is not
optional), subtracted against it, and `detect_sources` runs on the
resulting `S_corr` instead of the science image, extending detection to
sources below the single-frame noise floor. `threshold` is then a
detection significance in sigma directly, since `S_corr` has zero mean and
unit variance by construction; pixels outside either frame's footprint are
excluded from detection. Because reprojection puts every frame's
detections on `reference.wcs` rather than each frame's own WCS,
`astrometric_calibrate` is given `reference.wcs` for every frame in this
path. The empirical PSF is estimated from each frame *after*
reprojection, from the same resampled pixels `zogy_subtract` actually
differences — interpolation subtly reshapes a PSF, and measuring it
beforehand understated that mismatch, leaving systematic subtraction
residuals at every bright star (found via real ZTF data: ~85-97% of
excess `S_corr` detections in a real test field fell within 15 px of a
bright reference star). Bright stars in both the frame and `reference`
are also detected and passed to [`zogy_subtract`](@ref) as
`n_sources`/`r_sources`, so its astrometric-noise term can account for
residual sub-pixel misregistration at exactly those positions instead of
`S_corr` reporting an inflated, uncorrected significance there.

A per-frame quality gate applies in the `reference` path: if a frame's
`S_corr` standard deviation (plain `Statistics.std`, over its valid,
in-footprint region) exceeds `quality_max_std`, that frame is treated as
having found nothing (a warning is emitted) rather than contributing
possibly-spurious detections. `S_corr` is only *nominally* unit-variance —
real registration and PSF-matching imperfections inflate it — and on real
ZTF data, good frames measured std ~1.1-1.2 while one frame with an
apparent transient atmospheric issue (see git log `3e31b95`) measured
~2.0 and produced roughly twice as many bright detections as every other
frame. Plain `std`, not a robust (MAD-based) one, deliberately: that real
anomalous frame's excess showed up as ~232 point-like bright residuals
rather than a bulk noise increase, and a MAD-based spread — robust to
exactly that kind of minority-of-pixels outlier — measured ~0.997 for it,
indistinguishable from the good frames' ~0.98-1.10; plain `std` is what
actually separates them here. The default of `1.5` sits with margin on
both sides of that real gap; it is calibrated from that one dataset, not
a universal constant — tune per survey/conditions.

A gated frame contributes zero detections, and `link_candidates`
(via `min_frames`) requires every frame to match by default — so a
single real gated frame used to make no tracklet reachable at all unless
the caller manually lowered `min_frames`. `min_frames`'s default
(`nothing`) now accounts for this automatically: the value actually used
is `length(fits_paths)` minus however many frames got gated this run,
computed *after* detection, not the eager keyword default — pass
`min_frames` explicitly to override this and get the old, literal
behavior.

The raw (no-`reference`) path has no equivalent gate built into
`S_corr`'s own statistics, but a related signal is available there:
[`photometric_scale`](@ref)'s per-frame flux-scale factor (computed after
detection, from every frame's own detections). If any frame's factor
relative to the field's median deviates by more than
`photometric_outlier_threshold` (default `0.2`, i.e. 20%), a warning is
emitted — this is only a warning, never an automatic exclusion, unlike
`quality_max_std`, precisely to avoid repeating that gate's own
combinatorial side effect on `link_candidates`'s `min_frames` (see the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#The-quality-gate's-combinatorial-side-effect-on-tracklet-count)).

Tested directly (not left unvalidated) against the one real, confirmed
anomaly available: the same field-451 frame `quality_max_std` catches (a
likely passing cloud — `S_corr` std ~2.0 vs. ~1.1-1.2 for the other four,
and ~232 excess bright residuals). That frame's *raw-path*
`photometric_scale`, measured directly, differs from the field median by
only 0.64% — nowhere near the default 20% threshold, at any reasonable
setting of it. This is a real, informative negative result, not just "it
didn't flag, so it's unvalidated": it shows `photometric_outlier_threshold`
and `quality_max_std` are sensitive to genuinely different failure modes,
not two redundant checks for the same thing. A passing cloud during ZOGY
differencing shows up as elevated per-pixel residual noise in `S_corr` —
exactly what `quality_max_std` measures — but a *raw*, undifferenced
frame's stars can still read at normal relative brightness to each other
even under a cloud thin enough not to have caused uniform extinction
across the whole exposure, which is what `photometric_scale`'s
ensemble-ratio approach would need to see. `photometric_outlier_threshold`
is left in place for the failure mode it *would* catch (a real, uniform
per-frame flux-scale shift — heavier cloud, real transparency loss,
guiding/focus problems) — that specific scenario is still unvalidated,
since no real example of it has turned up in this project's data yet —
but it is no longer accurate to say this check has never been tested
against a real anomaly at all.

If a frame's header has no WCS, [`load_wcs`](@ref) raises an error;
passing `plate_solve_api_key` (a nova.astrometry.net API key) makes that
frame fall back to [`plate_solve`](@ref) instead of failing outright —
see its docstring for what that involves (a live network round trip,
polling until the frame solves or times out).

These options, all off by default and all only on the raw (no-`reference`)
path, make the pipeline usable on real IASC/Pan-STARRS1 image sets — each
was needed on a real 2025 IASC practice set (ps1-NewPractice_3, field
`XY54_p10`), where the defaults recovered neither of its two real moving
objects and reported positions ~7" off:

- `mask_fill_values`: mask constant-valued plateaus (the gaps between
  detector cells; see [`fill_value_mask`](@ref)) before detection.
- `refine_astrometry`: replace each frame's header WCS with one fitted to
  reference stars detected in the frame ([`refine_wcs`](@ref)) — against
  `reference_stars` if given, else Gaia DR3 stars covering the first
  frame, fetched once ([`gaia_reference_stars`](@ref), a live network
  request). On that set the header WCS was off by 6.8"; refined, an
  unknown object's positions agreed with Astrometrica's to 0.08-0.18".
- `difference_stack`: detect on each frame minus the median of the whole
  sequence ([`stack_difference`](@ref)) instead of on the raw frame —
  static sources cancel, so `threshold` can go much lower without
  drowning in stars. Needs at least 3 frames. Clustered detections
  (satellite trails, saturated-star residuals) are then dropped — see
  [`_drop_clustered`](@ref).
- `min_speed` (pixels per day) and `max_flux_ratio` filter tracklets
  after linking, per IASC's own "true signature" tests (see
  [`_filter_tracklets`](@ref)): use `max_flux_ratio=2.5` for IASC's
  1-magnitude limit.
- `min_sharpness`: drop single-pixel spikes — cosmic rays, hot pixels —
  at detection (see [`detect_sources`](@ref)); `0.3` on IASC data.

`fits_paths` must already be given in time order. `psf_threshold` and
`psf_min_separation` are forwarded to `estimate_psf` for each frame's own
PSF (only used when `reference` is given — `reference.psf` is estimated
separately, by the caller, when building it). All other keyword arguments
are forwarded to `detect_sources` and `link_candidates`.

Returns the candidate table from `astrometric_calibrate` (columns `id`,
`frame`, `x`, `y`, `ra`, `dec`, `epoch`), ready to pass to
[`crossmatch_catalog`](@ref). Cross-matching is not performed here, since
it requires network access and a choice of catalog and radius.
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
                       min_speed::Real=0.0, max_flux_ratio::Real=Inf, min_sharpness::Real=-Inf,
                       exptime_key::Union{Nothing,AbstractString}="EXPTIME")
    detections_per_frame, wcs_per_frame, timestamps, n_gated = _detect_all_frames(
        fits_paths; timestamp_key, threshold, box_size, aperture_radius,
        reference, psf_threshold, psf_min_separation, quality_max_std, plate_solve_api_key,
        photometric_outlier_threshold, mask_fill_values, difference_stack,
        refine_astrometry, reference_stars, min_sharpness, exptime_key)

    effective_min_frames = min_frames === nothing ? length(fits_paths) - n_gated : min_frames
    tracklets = link_candidates(detections_per_frame, timestamps;
                                 max_speed, match_radius, min_frames=effective_min_frames)
    tracklets = _filter_tracklets(tracklets, detections_per_frame, timestamps; min_speed, max_flux_ratio)
    return astrometric_calibrate(tracklets, wcs_per_frame, timestamps)
end

"""
    _detect_all_frames(fits_paths; timestamp_key, threshold, box_size, aperture_radius,
                        reference, psf_threshold, psf_min_separation,
                        quality_max_std, plate_solve_api_key)
        -> (detections_per_frame, wcs_per_frame, timestamps, n_gated)

The per-frame detection stage shared by [`run_pipeline`](@ref) (which
links these into movers) and [`search_field`](@ref) (which additionally
looks for stationary variable sources) — factored out so both share one
detection pass over `fits_paths` rather than repeating the expensive
reprojection/PSF/ZOGY work. See `run_pipeline`'s docstring for the
meaning of every keyword. On the raw (no-`reference`) path, each frame's
own `GAIN` header keyword (default `1.0` if absent) is passed to
[`detect_sources`](@ref) so `flux_err` includes source shot noise, not
just background noise — needed for [`find_variable_sources`](@ref)'s
chi-squared test to have a realistic error bar to test against. Not
applied on the `reference` (ZOGY) path, since `S_corr` there is already a
normalized detection-significance map, not physical counts. On that same
raw path, `photometric_outlier_threshold` is forwarded to a post-loop
[`photometric_scale`](@ref) check — see `run_pipeline`'s docstring.

`n_gated` is how many frames the `quality_max_std` gate excluded (always
`0` on the raw path, which has no such gate) — both callers use it to
auto-adjust their own `min_frames` default, since a gated frame otherwise
silently makes no tracklet/variable-candidate reachable at all.

`detections_per_frame` is concretely typed from `detect_sources`'s own
return type (confirmed identical on both its empty and non-empty return
paths, not assumed) rather than left as `Vector{Any}` — a real
type-stability fix found while profiling for other speedups, though it
didn't measurably change this function's own real-data timing; see
[Design refinements](https://richard7987.github.io/AsteroidPipeline.jl/dev/design-refinements).
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
                             exptime_key::Union{Nothing,AbstractString}="EXPTIME")
    reference === nothing || !(mask_fill_values || difference_stack || refine_astrometry) ||
        throw(ArgumentError("mask_fill_values, difference_stack and refine_astrometry apply only " *
                            "without `reference` (the ZOGY path has its own registration and masks)"))
    detections_per_frame = typeof(_empty_detections())[]
    wcs_per_frame = WCSTransform[]
    timestamps = Float64[]
    n_gated = 0
    # difference_stack needs every frame before it can detect in any of them
    stacked_images = Matrix{Float64}[]
    stacked_masks = Union{Nothing,BitMatrix}[]

    # Computed once, not per frame: reference.image is the same every
    # iteration, and this is what zogy_subtract's astrometric-noise term
    # matches each frame's own bright stars against.
    r_sources = reference === nothing ? nothing :
                detect_sources(permutedims(reference.image); threshold=psf_threshold)

    for path in fits_paths
        FITS(path, "r") do f
            hdu = f[1]

            # FITSIO.jl preserves FITS's column-major storage, so `read`
            # returns an (x, y)-ordered array; detect_sources (via
            # Photometry.PeakMesh) expects the (y, x) ordering used
            # throughout the rest of the pipeline. reference/zogy code
            # (src/reference.jl, src/zogy.jl) works in raw (x, y) order
            # throughout and only permutes at this same final step.
            raw = Float64.(read(hdu))
            frame_wcs = try
                load_wcs(read_header(hdu, String))
            catch e
                plate_solve_api_key === nothing && rethrow(e)
                plate_solve(path; api_key=plate_solve_api_key)
            end
            gain = haskey(read_header(hdu), "GAIN") ? read_key(hdu, "GAIN")[1] : 1.0

            push!(timestamps, frame_epoch(read_header(hdu); timestamp_key, exptime_key))

            if reference === nothing
                image = permutedims(raw)
                mask = mask_fill_values ? fill_value_mask(image) : nothing
                if refine_astrometry
                    reference_stars === nothing &&
                        (reference_stars = _field_reference_stars(image, frame_wcs, timestamps[end]))
                    refined = refine_wcs(image, frame_wcs, reference_stars; mask)
                    if refined.refined
                        frame_wcs = refined.wcs
                    else
                        @warn "astrometric refinement failed; keeping the header WCS" path n_matches=refined.n_matches
                    end
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
                # Reprojection legitimately leaves NaN at the border where
                # the footprint doesn't fully cover the target grid; left
                # in, a single NaN poisons zogy_subtract's whole-array
                # median/fft, turning all of S_corr to NaN. Filled with the
                # valid region's own level, not 0, since 0 would instead
                # bias zogy_subtract's background-median computation
                # (negligible at real ZTF scale, but not in general).
                resampled[.!frame_mask] .= median(resampled[frame_mask])

                # PSF and source positions measured on the resampled
                # (post-reprojection) pixels, matching what zogy_subtract
                # actually differences (see the docstring above).
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
                # No gain here: S_corr is already a normalized
                # detection-significance map, not physical counts, so
                # Poisson noise doesn't apply to it.
                push!(detections_per_frame, detect_sources(image; threshold, box_size, aperture_radius))
            end
        end
    end

    if difference_stack
        masks = all(isnothing, stacked_masks) ? nothing : stacked_masks
        differences, diff_masks = stack_difference(stacked_images, wcs_per_frame; masks)
        # No gain either: a difference image's noise is the frame's and the
        # median's together, and a mover's own shot noise is negligible at
        # the faint end this mode exists for.
        for (diff, m) in zip(differences, diff_masks)
            dets = detect_sources(diff; threshold, box_size, aperture_radius, mask=m, min_sharpness)
            push!(detections_per_frame, _drop_clustered(dets, _CLUSTER_RADIUS))
        end
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

    return detections_per_frame, wcs_per_frame, timestamps, n_gated
end

"""
    search_field(fits_paths; <all run_pipeline keywords>,
                 variability_position_tolerance::Real=2.0,
                 variability_min_frames::Union{Nothing,Integer}=nothing,
                 variability_chi2_threshold::Real=10.0,
                 variability_systematic_error_fraction::Real=0.03)
        -> (movers=<table>, variables=<table>)

Run [`run_pipeline`](@ref)'s asteroid-candidate search and
[`find_variable_sources`](@ref)'s variable/transient search together,
sharing a single detection pass over `fits_paths` (via
[`_detect_all_frames`](@ref)) rather than repeating the expensive
per-frame work (reprojection, PSF estimation, ZOGY differencing) twice —
the efficient choice for a real observing run that wants both outputs,
since `run_pipeline` alone would need a second full pass over the same
frames to also find variables.

All keywords through `photometric_outlier_threshold` are exactly
`run_pipeline`'s (see its docstring, including how `min_frames`'s default
now accounts for quality-gated frames); `variability_position_tolerance`,
`variability_min_frames`, `variability_chi2_threshold`,
`variability_normalize`, `variability_max_relative_error`, and
`variability_systematic_error_fraction` are forwarded to
[`find_variable_sources`](@ref) as its `position_tolerance`, `min_frames`,
`chi2_threshold`, `normalize`, `max_relative_error`, and
`systematic_error_fraction`. `variability_min_frames` defaults the same
way `min_frames` does — `nothing` means "every non-gated frame must
match". `variability_normalize` defaults to `reference === nothing` — off
on the ZOGY path, per `find_variable_sources`'s own docstring (`S_corr`
isn't on a physical flux scale, so ensemble photometric normalization
doesn't apply there).

Returns a named tuple `(movers=..., variables=...)`, each an
`astrometric_calibrate` candidate table (columns `id`, `frame`, `x`, `y`,
`ra`, `dec`, `epoch`) — `movers` is identical to what `run_pipeline` would
return for the same arguments; `variables` is ready for
`crossmatch_catalog(...; :vsx)` the same way. `run_pipeline` itself
remains the right entry point when only movers are needed.
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
    detections_per_frame, wcs_per_frame, timestamps, n_gated = _detect_all_frames(
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

[`link_candidates`](@ref) and [`find_variable_sources`](@ref) both
compare positions across frames in pixels and assume one shared grid;
the ZOGY path gets that from reprojecting onto `reference`, but the raw
path used to link each frame's *own* pixel positions directly. Frames
dither between exposures, so a real object's straight-line track was
bent by each frame's pointing offset. Found on a real IASC set (field
`XY54_p10`, dithers up to (6, 9) px): a G≈21 object detected in all four
frames was never linked, because those offsets exceeded the 2"
(7.8 px) `match_radius` tuned to PS1's astrometric precision. Only
*relative* WCS accuracy matters here — the header's 6.8" absolute error
on that same set is common to all four frames and cancels.

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

# Neighbourhood radius (pixels) for `_drop_clustered` on the difference_stack path.
const _CLUSTER_RADIUS = 15.0
