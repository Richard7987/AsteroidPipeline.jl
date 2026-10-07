"""
    detect_sources(image::AbstractMatrix{<:Real}; threshold::Real,
                    box_size::NTuple{2,Integer}=(5, 5), aperture_radius::Real=3.0,
                    gain::Union{Nothing,Real}=nothing,
                    mask::Union{Nothing,AbstractMatrix{Bool}}=nothing, min_sharpness::Real=-Inf)

Detect point sources in a FITS frame.

The background level and RMS noise are estimated with
`BackgroundMeshes.estimate_background` (`SourceExtractorBackground` location,
`MADStdRMS` scale — the same estimators used by SourceExtractor/`photutils`).
Local maxima at least `threshold` sigma above the background, within grid
boxes of `box_size` pixels, are extracted with `Photometry.PeakMesh`. Each
detection's flux is then measured with circular aperture photometry of
radius `aperture_radius` pixels on the background-subtracted image.

Each aperture is centered on a sub-pixel-refined position (a flux-weighted
centroid within `aperture_radius` of `PeakMesh`'s own integer-pixel peak
— falls back to that raw peak whenever the refinement isn't trustworthy;
see [`_refine_centroid`](@ref)), not the raw integer peak itself: a real
star's true position doesn't move frame to frame, but which *integer*
pixel reads highest does, under noise — and since `aperture_radius` is
comparable to a typical PSF scale, that alone changes how much flux a
fixed aperture encloses, frame to frame, for no real reason. This was the
actual cause of a real, measured false-positive floor in
[`find_variable_sources`](@ref)'s variability test (see the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#The-centroid-fix-barely-moved-the-false-positive-floor-—-the-real-cause-was-a-systematic-error-floor)) — the refinement only affects where the aperture
is centered, never the returned `x`/`y` (still `PeakMesh`'s own integer
position, unchanged). That refined position is returned separately, as
`xcen`/`ycen`: for astrometry the integer peak alone is a real error
source — up to half a pixel, i.e. ~0.13" at Pan-STARRS1's 0.257"/px,
comparable to PS1's own whole astrometric solution (`CERROR` ~0.06" in
real IASC headers). [`link_candidates`](@ref) uses `xcen`/`ycen` when a
table has them.

`mask` (same size as `image`, `true` = invalid) excludes pixels that
carry no real sky signal: they are left out of the background/noise
estimate and set to the background level before peak finding, so they
can neither be detected nor bias the noise. Real IASC/Pan-STARRS1 frames
need this — the gaps between the detector's cells are filled with one
constant value (not `BLANK`, not `NaN`), ~14% of every frame in a real
2025 practice set (see [`fill_value_mask`](@ref)).

`sharpness` is the mean of a peak's four direct neighbours over the peak
itself (background-subtracted): a star or asteroid spreads its light
over several pixels (real IASC objects measured 0.46-0.82 at PS1's
seeing), a cosmic ray or hot pixel doesn't (≈0). Peaks below
`min_sharpness` are dropped — off by default. On a real 2019 IASC set
(`XY14_p10`), single-pixel spikes were 94-98% of every frame's 6σ
detections on [`stack_difference`](@ref) images, about 3,800 in its worst
frame; `min_sharpness=0.3` removes them and keeps every real object.

Returns a table with columns `x`, `y` (pixel position), `peak` (background-
subtracted peak pixel value), `flux` (aperture sum), `flux_err`, and
`xcen`, `ycen` (the sub-pixel centroid, in the same sense as `x`/`y`), and
`sharpness` (see below). By
default (`gain=nothing`) `flux_err` is `Photometry.photometry`'s propagated
aperture error from the uniform per-pixel `noise` used for detection alone
— not a full per-pixel variance map, the same caveat [`light_curve`](@ref)
documents. Passing `gain` (electrons/ADU, from the frame's own `GAIN`
header keyword) adds each pixel's own Poisson (shot) noise,
`sqrt(noise^2 + max(pixel, 0) / gain)`: background noise alone
underestimates a bright star's real flux uncertainty — on real ZTF data
this barely moves the *median* `flux_err/flux` (2.97% vs 3.16%, most
detections being near-threshold and background-noise-dominated either
way), but for the single brightest star in that same frame, shot noise
was ~10x the background-only estimate (0.039% vs 0.004%) — exactly where
underestimating the error bar would most distort
[`find_variable_sources`](@ref)'s chi-squared test. Left `nothing`
(background noise only) for [`link_candidates`](@ref), which only uses
positions and has no use for a flux uncertainty at all.
"""
function detect_sources(image::AbstractMatrix{<:Real}; threshold::Real,
                         box_size::NTuple{2,<:Integer}=(5, 5), aperture_radius::Real=3.0,
                         gain::Union{Nothing,Real}=nothing,
                         mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                         min_sharpness::Real=-Inf)
    if mask === nothing
        background, noise = estimate_background(image; location=SourceExtractorBackground(), rms=MADStdRMS())
    else
        size(mask) == size(image) || throw(DimensionMismatch("mask must have the same size as image"))
        valid = image[.!mask]
        isempty(valid) && return _empty_detections()
        background, noise = estimate_background(valid; location=SourceExtractorBackground(), rms=MADStdRMS())
    end
    subtracted = image .- background
    # Masked pixels sit exactly at the background level after subtraction,
    # so they can neither produce a peak nor add flux to a nearby aperture.
    mask === nothing || (subtracted[mask] .= 0.0)

    noise <= 0 && return _empty_detections()

    finder = PeakMesh(box_size, threshold)
    detection_error_map = fill(noise, size(image))
    peaks = extract_sources(finder, subtracted, detection_error_map)
    sharpness = [_sharpness(subtracted, row.y, row.x) for row in peaks]
    if isfinite(min_sharpness)
        keep = sharpness .>= min_sharpness   # NaN (array edge) never passes
        peaks, sharpness = peaks[keep], sharpness[keep]
    end

    # `PeakMesh` reports x/y in the standard Cartesian sense (x=column,
    # i.e. the array's 2nd dimension; y=row, the 1st — see
    # Photometry.jl's own `extract_sources`, `to_nt(ci) = (x=ci[2],
    # y=ci[1], ...)`). `CircularAperture`, in the very same package,
    # does the opposite internally (its `x` field indexes the array's
    # *1st* dimension, `y` the 2nd — see `bounds`/`overlap` in
    # Photometry.jl's circular.jl). Passing `(row.x, row.y)` straight
    # through silently centers the aperture at the transposed pixel
    # whenever the true position isn't on the row==column diagonal —
    # confirmed by comparing a measured flux against the analytic
    # enclosed-energy integral for an isolated Gaussian, which came back
    # ~140x too small before this swap and matched to ~1% after it.
    # `(row.y, row.x)` here is the fix (dim1, dim2 order); `_refine_centroid`
    # takes and returns positions in that same (dim1, dim2) order.
    centroids = [_refine_centroid(subtracted, row.y, row.x, aperture_radius) for row in peaks]
    apertures = [CircularAperture(d1, d2, aperture_radius) for (d1, d2) in centroids]
    photom_error_map = gain === nothing ? detection_error_map :
                        sqrt.(noise^2 .+ max.(subtracted, 0.0) ./ gain)
    photom = isempty(apertures) ? nothing : photometry(apertures, subtracted, photom_error_map)
    flux = isempty(apertures) ? Float64[] : [row.aperture_sum for row in photom]
    flux_err = isempty(apertures) ? Float64[] : [row.aperture_sum_err for row in photom]
    xcen = Float64[d2 for (_, d2) in centroids]
    ycen = Float64[d1 for (d1, _) in centroids]

    return Table(x=peaks.x, y=peaks.y, peak=peaks.value, flux=flux, flux_err=flux_err, xcen=xcen, ycen=ycen,
                 sharpness=sharpness)
end

_empty_detections() = Table(x=Int[], y=Int[], peak=Float64[], flux=Float64[], flux_err=Float64[],
                            xcen=Float64[], ycen=Float64[], sharpness=Float64[])

"""
    _refine_centroid(subtracted, d1, d2, radius) -> (Float64, Float64)

Refine an integer pixel position `(d1, d2)` (indices into `subtracted`'s
1st and 2nd dimensions respectively) to a sub-pixel flux-weighted
centroid within a `radius`-pixel window around it — the same
first-moment technique [`estimate_psf`](@ref) already uses per star
stamp, applied here per detection instead.

Falls back to the unmodified `(d1, d2)` (as `Float64`s) whenever the
refinement can't be trusted: the window would run off the array edge,
the window's total flux isn't positive, or the computed shift exceeds
`radius` itself (a sign the "centroid" is being pulled toward a neighbor
or noise excursion, not the true peak) — never makes the position worse
than the raw integer peak it started from.
"""
function _refine_centroid(subtracted::AbstractMatrix{<:Real}, d1::Integer, d2::Integer, radius::Real)
    r = max(1, round(Int, radius))
    n1, n2 = size(subtracted)
    (r < d1 <= n1 - r && r < d2 <= n2 - r) || return Float64(d1), Float64(d2)

    stamp = subtracted[d1-r:d1+r, d2-r:d2+r]
    total = sum(stamp)
    total > 0 || return Float64(d1), Float64(d2)

    offs = Float64.(-r:r)
    delta1 = sum(offs .* sum(stamp, dims=2)[:]) / total
    delta2 = sum(offs .* sum(stamp, dims=1)[:]) / total
    (abs(delta1) <= radius && abs(delta2) <= radius) || return Float64(d1), Float64(d2)

    return d1 + delta1, d2 + delta2
end

"""
    fill_value_mask(image::AbstractMatrix{<:Real}; dilate::Integer=1) -> BitMatrix

Mask (`true` = invalid) of every constant-valued plateau in `image`: any
pixel whose whole 3x3 neighbourhood shares exactly its value, plus that
neighbourhood, grown by `dilate` more pixels. Real sky never does this —
even a smooth background carries per-pixel noise — so a plateau is
always something written into the frame rather than measured: a fill
value, or a saturated star's clipped core.

Real IASC/Pan-STARRS1 frames need it: the gaps between the detector's
cells are filled with one constant value (160 ADU in a real 2025
practice set, ~14% of every frame), marked neither with `BLANK` nor
`NaN`, so nothing else flags them. Unmasked, each gap edge is a sharp
step that background estimation and peak finding treat as real signal,
and an object crossing a gap simply vanishes for those frames —
confirmed on the same set, where the known asteroid 2018 LT's
SkyBoT-predicted track ran through a gap. Detecting plateaus instead of
assuming a particular fill value keeps this survey-agnostic.

Pass the result to [`detect_sources`](@ref) as `mask`.
"""
function fill_value_mask(image::AbstractMatrix{<:Real}; dilate::Integer=1)
    n1, n2 = size(image)
    mask = falses(n1, n2)
    r = 1 + dilate
    for j in 2:n2-1, i in 2:n1-1
        v = image[i, j]
        all(image[i+a, j+b] == v for a in -1:1, b in -1:1) || continue
        mask[max(1, i-r):min(n1, i+r), max(1, j-r):min(n2, j+r)] .= true
    end
    return mask
end

# Mean of the four direct neighbours of `(d1, d2)` over its own value;
# NaN at the array edge or for a non-positive peak.
function _sharpness(img::AbstractMatrix, d1::Integer, d2::Integer)
    n1, n2 = size(img)
    (1 < d1 < n1 && 1 < d2 < n2) || return NaN
    v = img[d1, d2]
    v > 0 || return NaN
    return (img[d1 - 1, d2] + img[d1 + 1, d2] + img[d1, d2 - 1] + img[d1, d2 + 1]) / 4v
end
