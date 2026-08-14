"""
    run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                 timestamp_key::AbstractString="MJD-OBS",
                 threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                 aperture_radius::Real=3.0, max_speed::Real=Inf,
                 match_radius::Real=2.0, min_frames::Integer=length(fits_paths))

Run the local stages of the pipeline — detection, linking, and
astrometric calibration — on a time-ordered sequence of FITS frames from
the same field.

For each frame: the primary image is read and detected on with
[`detect_sources`](@ref); the primary header supplies the WCS solution
(via [`load_wcs`](@ref)) and the observation epoch, taken from the
`timestamp_key` keyword as a Modified Julian Date and converted to
Julian Date. Detections are then linked across frames with
[`link_candidates`](@ref) and calibrated to sky coordinates with
[`astrometric_calibrate`](@ref).

`fits_paths` must already be given in time order. All other keyword
arguments are forwarded to `detect_sources` and `link_candidates`.

Returns the candidate table from `astrometric_calibrate` (columns `id`,
`frame`, `x`, `y`, `ra`, `dec`, `epoch`), ready to pass to
[`crossmatch_catalog`](@ref). Cross-matching is not performed here, since
it requires network access and a choice of catalog and radius.
"""
function run_pipeline(fits_paths::AbstractVector{<:AbstractString};
                       timestamp_key::AbstractString="MJD-OBS",
                       threshold::Real=5.0, box_size::NTuple{2,<:Integer}=(5, 5),
                       aperture_radius::Real=3.0, max_speed::Real=Inf,
                       match_radius::Real=2.0, min_frames::Integer=length(fits_paths))
    detections_per_frame = []
    wcs_per_frame = WCSTransform[]
    timestamps = Float64[]

    for path in fits_paths
        FITS(path, "r") do f
            hdu = f[1]

            # FITSIO.jl preserves FITS's column-major storage, so `read`
            # returns an (x, y)-ordered array; detect_sources (via
            # Photometry.PeakMesh) expects the (y, x) ordering used
            # throughout the rest of the pipeline.
            image = permutedims(read(hdu))
            push!(detections_per_frame, detect_sources(image; threshold, box_size, aperture_radius))

            push!(wcs_per_frame, load_wcs(read_header(hdu, String)))

            mjd = read_key(hdu, timestamp_key)[1]
            push!(timestamps, mjd + 2400000.5)
        end
    end

    tracklets = link_candidates(detections_per_frame, timestamps; max_speed, match_radius, min_frames)
    return astrometric_calibrate(tracklets, wcs_per_frame, timestamps)
end
