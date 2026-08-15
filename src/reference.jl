"""
    build_reference(frames, target_wcs, shape) -> (image, sigma, mask)

Build a deep, static-sky reference image on the pixel grid of `target_wcs`
from `frames`, each a `(image, wcs, magzp, sigma)` tuple (raw FITS pixel
order — see [`load_frame`](@ref)).

Every frame is reprojected onto the target grid with `Reproject.reproject`
(bilinear interpolation; each frame's own WCS legitimately differs from
`target_wcs` by a few pixels of dither even within one night), rescaled to
`target_magzp` by `10^(-0.4*(magzp - target_magzp))`, and combined
per pixel with the **median** — chosen specifically because it rejects
whatever moved between epochs (asteroids, satellite trails, cosmic rays),
which a mean would instead bake into the reference as ghost artifacts.

`sigma` is the reference's per-pixel background RMS, propagated from each
frame's own (pre-reprojection) noise estimate and combined as
`1.2533 * median(sigma_i) / sqrt(n)`, the asymptotic standard error of the
median under Gaussian noise (efficiency `sqrt(pi/2)` relative to the mean).
This is an approximation — reprojection correlates adjacent-pixel noise via
interpolation, which this does not model — adequate for ZOGY's per-pixel
normalization but not for e.g. certifying it as a strict statistical
guarantee.

`mask` is `true` where at least one input frame had a valid (in-footprint,
non-NaN) sample at that pixel; elsewhere `image` is filled with 0.
Reprojection legitimately leaves a thin NaN border where a frame's rotated
footprint doesn't fully cover the target grid — real, not a bug — and
callers must exclude `!mask` pixels from detection.
"""
function build_reference(frames, target_wcs::WCSTransform, shape::NTuple{2,Integer})
    target_magzp = frames[1].magzp

    stack = Array{Float64}(undef, shape..., length(frames))
    valid = falses(shape..., length(frames))
    sigmas = Float64[]

    for (k, frame) in enumerate(frames)
        resampled, frame_mask = Reproject.reproject((frame.image, frame.wcs), target_wcs; shape_out=shape)
        scale = 10.0^(-0.4 * (frame.magzp - target_magzp))
        stack[:, :, k] .= resampled .* scale
        valid[:, :, k] .= frame_mask
        push!(sigmas, frame.sigma * scale)
    end

    image = zeros(Float64, shape)
    mask = falses(shape)
    for ci in CartesianIndices(image)
        samples = [stack[ci, k] for k in 1:length(frames) if valid[ci, k]]
        if !isempty(samples)
            image[ci] = median(samples)
            mask[ci] = true
        end
    end

    sigma = 1.2533 * median(sigmas) / sqrt(length(frames))
    return image, sigma, mask
end

"""
    load_frame(path) -> (image, wcs, magzp, sigma)

Read a FITS science frame for use with [`build_reference`](@ref) or ZOGY
differencing. `image` is in raw FITS pixel order (`(NAXIS1, NAXIS2)`,
**not** [`detect_sources`](@ref)'s `(y, x)` convention — see the note in
`src/pipeline.jl`), matching what `Reproject.jl` and `WCS.jl` both expect.
`sigma` is the background RMS estimated on the native (pre-reprojection)
pixel grid, since resampling correlates neighbouring pixels and would bias
a naive RMS computed after the fact.
"""
function load_frame(path::AbstractString)
    return FITS(path, "r") do f
        hdu = f[1]
        image = Float64.(read(hdu))
        wcs = load_wcs(read_header(hdu, String))
        magzp = read_key(hdu, "MAGZP")[1]
        _, sigma = estimate_background(permutedims(image); location=SourceExtractorBackground(), rms=MADStdRMS())
        (image=image, wcs=wcs, magzp=magzp, sigma=sigma)
    end
end
