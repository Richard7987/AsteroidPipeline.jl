"""
    run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                 timestamp_key::AbstractString="MJD-OBS",
                 threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                 aperture_radius::Real=3.0, max_speed::Real=Inf,
                 match_radius::Real=2.0, min_frames::Integer=length(fits_paths),
                 reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0)

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
path. The empirical PSF is estimated from each frame *before*
reprojection (interpolation would otherwise bias its shape), which means
`psf_n` describes the pre-reprojection image while `zogy_subtract` is fed
the post-reprojection one — a known approximation, adequate because
`reference.psf` was built the same way (`build_reference` also reprojects
before this pipeline ever measures a PSF from it).

`fits_paths` must already be given in time order. `psf_threshold` and
`psf_min_separation` are forwarded to `estimate_psf` for each frame's own
PSF (only used when `reference` is given — `reference.psf` is estimated
separately, by the caller, when building it). All other keyword arguments
are forwarded to `detect_sources` and `link_candidates`.

TODO: `link_candidates` derives its trial velocity from frames 1 and 2
only (see its docstring); closely spaced first frames amplify velocity
error when extrapolated across a long time baseline. A robust fit across
all frames would remove this sensitivity to frame spacing and ordering.

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
                       reference=nothing, psf_threshold::Real=20.0, psf_min_separation::Real=40.0)
    detections_per_frame = []
    wcs_per_frame = WCSTransform[]
    timestamps = Float64[]

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
            frame_wcs = load_wcs(read_header(hdu, String))
            gain = haskey(read_header(hdu), "GAIN") ? read_key(hdu, "GAIN")[1] : 1.0

            if reference === nothing
                image = permutedims(raw)
                push!(wcs_per_frame, frame_wcs)
            else
                _, sigma_n = estimate_background(permutedims(raw);
                                                  location=SourceExtractorBackground(), rms=MADStdRMS())
                psf_n = estimate_psf(raw; threshold=psf_threshold, min_separation=psf_min_separation)

                resampled, frame_mask = Reproject.reproject((raw, frame_wcs), reference.wcs;
                                                              shape_out=size(reference.image))
                # Reprojection legitimately leaves NaN at the thin border
                # where the rotated/dithered footprint doesn't fully cover
                # the target grid (frame_mask marks exactly this). Left
                # in, a single NaN poisons zogy_subtract's whole-array
                # median/fft computations, turning the *entire* S_corr
                # into NaN — the fill value itself doesn't matter for
                # detection correctness (these pixels are excluded via
                # frame_mask below regardless), only that it isn't NaN.
                resampled[.!frame_mask] .= 0.0
                _, s_corr = zogy_subtract(resampled, reference.image;
                                           psf_n=psf_n, psf_r=reference.psf,
                                           sigma_n=sigma_n, sigma_r=reference.sigma, gain_n=gain)
                image = permutedims(s_corr .* frame_mask .* reference.mask)
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
