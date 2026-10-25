"""
    build_reference(frames, target_wcs, shape; workers=nothing) -> (image, sigma, mask)

Build a deep, static-sky reference image on the pixel grid of `target_wcs`
from `frames`, each a `(image, wcs, magzp, sigma)` tuple (raw FITS pixel
order — see [`load_frame`](@ref)).

Every frame is reprojected onto the target grid with `Reproject.reproject`
(bilinear), rescaled to `target_magzp` by `10^(-0.4*(magzp - target_magzp))`,
and combined per pixel with the **median**, which rejects whatever moved
between epochs (asteroids, satellite trails, cosmic rays). Reprojection
dominates the runtime.

By default (`workers=nothing`) frames are reprojected one at a time: not
with `Threads.@threads`, since `wcslib`, which `WCS.jl` wraps, isn't
thread-safe. With `workers` (process ids from `Distributed.addprocs`, each
with `@everywhere using AsteroidPipeline` already done) frames go through
`pmap`; a `WCSTransform` holds C pointers that can't cross processes, so
each WCS is sent as a FITS header string and rebuilt on the worker. On a
30-frame ZTF reference, 8 workers ran 3.2× faster with identical output.
`build_reference` never starts or stops workers itself.

`sigma` is the reference's per-pixel background RMS,
`1.2533 * median(sigma_i) / sqrt(n)`: the standard error of the median
under Gaussian noise. It ignores the noise correlation interpolation
introduces, which is adequate for ZOGY's normalization.

`mask` is `true` where at least one frame had a valid sample; elsewhere
`image` is 0. A frame's rotated footprint leaves a thin border uncovered,
so callers must exclude `!mask` pixels from detection.
"""
function build_reference(frames, target_wcs::WCSTransform, shape::NTuple{2,Integer};
                          workers::Union{Nothing,AbstractVector{Int}}=nothing)
    target_magzp = frames[1].magzp
    nframes = length(frames)

    stack = Array{Float64}(undef, shape..., nframes)
    valid = falses(shape..., nframes)
    sigmas = Vector{Float64}(undef, nframes)

    if workers === nothing
        # not Threads.@threads: wcslib isn't thread-safe
        for k in 1:nframes
            frame = frames[k]
            resampled, frame_mask = Reproject.reproject((frame.image, frame.wcs), target_wcs; shape_out=shape)
            scale = 10.0^(-0.4 * (frame.magzp - target_magzp))
            stack[:, :, k] .= resampled .* scale
            valid[:, :, k] .= frame_mask
            sigmas[k] = frame.sigma * scale
        end
    else
        # a WCSTransform can't cross processes; its header string can
        target_header = WCS.to_header(target_wcs)
        jobs = [(frames[k].image, WCS.to_header(frames[k].wcs), frames[k].magzp) for k in 1:nframes]
        pool = WorkerPool(collect(workers))
        results = pmap(pool, jobs) do job
            image, wcs_header, magzp = job
            _reproject_and_scale(image, wcs_header, magzp, target_header, shape, target_magzp)
        end
        for k in 1:nframes
            resampled, frame_mask = results[k]
            stack[:, :, k] .= resampled
            valid[:, :, k] .= frame_mask
            sigmas[k] = frames[k].sigma * 10.0^(-0.4 * (frames[k].magzp - target_magzp))
        end
    end

    image = zeros(Float64, shape)
    mask = falses(shape)
    buffer = Vector{Float64}(undef, length(frames))
    for ci in CartesianIndices(image)
        n = 0
        for k in 1:length(frames)
            if valid[ci, k]
                n += 1
                buffer[n] = stack[ci, k]
            end
        end
        if n > 0
            image[ci] = median(@view buffer[1:n])
            mask[ci] = true
        end
    end

    sigma = 1.2533 * median(sigmas) / sqrt(length(frames))
    return image, sigma, mask
end

"""
    _reproject_and_scale(image, wcs_header, magzp, target_header, shape, target_magzp)

Worker-side helper for [`build_reference`](@ref)'s `workers` path: rebuilds
both WCS solutions locally from plain header strings (never receives a
`WCSTransform` itself — see `build_reference`'s docstring for why) before
calling `Reproject.reproject`.
"""
function _reproject_and_scale(image, wcs_header::AbstractString, magzp::Real,
                               target_header::AbstractString, shape, target_magzp::Real)
    wcs = load_wcs(wcs_header)
    target_wcs = load_wcs(target_header)
    resampled, frame_mask = Reproject.reproject((image, wcs), target_wcs; shape_out=shape)
    scale = 10.0^(-0.4 * (magzp - target_magzp))
    return resampled .* scale, frame_mask
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
