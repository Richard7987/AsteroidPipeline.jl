"""
    detect_sources(image::AbstractMatrix{<:Real}; threshold::Real,
                    box_size::NTuple{2,Integer}=(5, 5), aperture_radius::Real=3.0)

Detect point sources in a FITS frame.

The background level and RMS noise are estimated with
`BackgroundMeshes.estimate_background` (`SourceExtractorBackground` location,
`MADStdRMS` scale — the same estimators used by SourceExtractor/`photutils`).
Local maxima at least `threshold` sigma above the background, within grid
boxes of `box_size` pixels, are extracted with `Photometry.PeakMesh`. Each
detection's flux is then measured with circular aperture photometry of
radius `aperture_radius` pixels on the background-subtracted image.

Returns a table with columns `x`, `y` (pixel position), `peak` (background-
subtracted peak pixel value), `flux` (aperture sum), and `flux_err`
(`Photometry.photometry`'s propagated aperture error from the uniform
per-pixel `noise` used for detection — not a full per-pixel variance map,
the same caveat [`light_curve`](@ref) documents), to be consumed by
[`link_candidates`](@ref) for inter-frame motion matching or
[`find_variable_sources`](@ref) for flux-variability matching.
"""
function detect_sources(image::AbstractMatrix{<:Real}; threshold::Real,
                         box_size::NTuple{2,<:Integer}=(5, 5), aperture_radius::Real=3.0)
    background, noise = estimate_background(image; location=SourceExtractorBackground(), rms=MADStdRMS())
    subtracted = image .- background

    if noise <= 0
        return Table(x=Int[], y=Int[], peak=Float64[], flux=Float64[], flux_err=Float64[])
    end

    finder = PeakMesh(box_size, threshold)
    error_map = fill(noise, size(image))
    peaks = extract_sources(finder, subtracted, error_map)

    apertures = [CircularAperture(row.x, row.y, aperture_radius) for row in peaks]
    photom = isempty(apertures) ? nothing : photometry(apertures, subtracted, error_map)
    flux = isempty(apertures) ? Float64[] : [row.aperture_sum for row in photom]
    flux_err = isempty(apertures) ? Float64[] : [row.aperture_sum_err for row in photom]

    return Table(x=peaks.x, y=peaks.y, peak=peaks.value, flux=flux, flux_err=flux_err)
end
