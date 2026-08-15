"""
    run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                 timestamp_key::AbstractString="MJD-OBS",
                 threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                 aperture_radius::Real=3.0, max_speed::Real=Inf,
                 match_radius::Real=2.0, min_frames::Integer=length(fits_paths),
                 reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                 quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing)

Run the local stages of the pipeline — detection, linking, and
astrometric calibration — on a time-ordered sequence of FITS frames from
the same field.

For each frame: the primary image is read; the primary header supplies the
WCS solution (via [`load_wcs`](@ref)) and the observation epoch, taken
from the `timestamp_key` keyword as a Modified Julian Date and converted
to Julian Date. Detections are then linked across frames with
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

If a frame's header has no WCS, [`load_wcs`](@ref) raises an error;
passing `plate_solve_api_key` (a nova.astrometry.net API key) makes that
frame fall back to [`plate_solve`](@ref) instead of failing outright —
see its docstring for what that involves (a live network round trip,
polling until the frame solves or times out).

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
                       timestamp_key::AbstractString="MJD-OBS",
                       threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                       aperture_radius::Real=3.0, max_speed::Real=Inf,
                       match_radius::Real=2.0, min_frames::Integer=length(fits_paths),
                       reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0,
                       quality_max_std::Real=1.5, plate_solve_api_key::Union{Nothing,AbstractString}=nothing)
    detections_per_frame = []
    wcs_per_frame = WCSTransform[]
    timestamps = Float64[]

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

            if reference === nothing
                image = permutedims(raw)
                push!(wcs_per_frame, frame_wcs)
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
                else
                    image = permutedims(s_corr .* valid)
                end
                push!(wcs_per_frame, reference.wcs)
            end
            push!(detections_per_frame, detect_sources(image; threshold, box_size, aperture_radius))

            mjd = read_key(hdu, timestamp_key)[1]
            push!(timestamps, mjd + 2400000.5)
        end
    end

    tracklets = link_candidates(detections_per_frame, timestamps; max_speed, match_radius, min_frames)
    return astrometric_calibrate(tracklets, wcs_per_frame, timestamps)
end
