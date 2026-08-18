"""
    build_reference(frames, target_wcs, shape; workers=nothing) -> (image, sigma, mask)

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
Reprojection dominates this function's real runtime by roughly two
orders of magnitude over the per-pixel combine step below (~24s/frame
vs. ~1s total, on real ZTF data; see the Investigation Log).

By default (`workers=nothing`) reprojection runs sequentially, one frame
at a time. Parallelizing this loop with `Threads.@threads` was tried and
reverted after it crashed (a real segfault, confirmed via a real
multi-threaded run on real data, not a hypothetical): `Reproject.reproject`
itself has no shared mutable state, but it calls into `WCS.jl`'s
`pix_to_world!`, which wraps `wcslib` (a C library) via `ccall` — and
concurrent calls into that library from multiple Julia **threads**
(sharing one process) are not safe. Checking a Julia package's own source
for global state, as was done here, is not sufficient to establish
thread-safety when it wraps a C library; the transitive dependency needs
the same scrutiny, which this hadn't had until the crash forced it.

Multiple **processes** (`Distributed.jl`) sidestep that specific hazard —
each process has its own independent `wcslib` state — but naively passing
a `WCSTransform` itself to a worker (e.g. inside a `pmap` closure) still
segfaults: `WCSTransform` holds pointers into `wcslib`-allocated C memory
that's only valid in the process that created it, and Julia's generic
serialization doesn't know to reconstruct that state on the receiving
end (confirmed via a real crash inside `WCS.jl`'s `getproperty`/
`convert_string` on a worker, deserializing a `WCSTransform` sent from
the main process). Passing `workers` here avoids that: each frame's WCS
is converted to a plain FITS header string (`WCS.to_header`, no pointers,
serializes safely) before being sent, and each worker reconstructs its
own local `WCSTransform` from that string (`load_wcs`) before calling
`Reproject.reproject` — confirmed on real data (real 30-frame ZTF
field-451 reference set, 8 worker processes) to run without crashing,
in 331.8s vs. an ~720s sequential baseline (the ~24s/frame figure above)
— about 2.2x, not full 8x core-count scaling, since each worker still
does its own `load_wcs` parsing per call and `pmap`'s own scheduling and
serialization overhead isn't free. `workers` must be pre-existing worker
process ids (e.g. from `Distributed.addprocs`), each of which the caller
must have already loaded this package on (`@everywhere using
AsteroidPipeline`) — `build_reference` itself never spawns or manages
worker processes, since process-pool lifecycle is an environment concern
the caller controls, not something a data-processing function should own
as a side effect.

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

The per-pixel median-combine below uses one reusable buffer rather than
allocating a fresh array per pixel — a real, measured 2x speedup on real
data (1.39s → 0.68s, 9.46M → 4.20M allocations; see the Investigation
Log for the full comparison), though small next to reprojection's own
cost above.
"""
function build_reference(frames, target_wcs::WCSTransform, shape::NTuple{2,Integer};
                          workers::Union{Nothing,AbstractVector{Int}}=nothing)
    target_magzp = frames[1].magzp
    nframes = length(frames)

    stack = Array{Float64}(undef, shape..., nframes)
    valid = falses(shape..., nframes)
    sigmas = Vector{Float64}(undef, nframes)

    if workers === nothing
        # Sequential, not Threads.@threads — see the docstring above.
        for k in 1:nframes
            frame = frames[k]
            resampled, frame_mask = Reproject.reproject((frame.image, frame.wcs), target_wcs; shape_out=shape)
            scale = 10.0^(-0.4 * (frame.magzp - target_magzp))
            stack[:, :, k] .= resampled .* scale
            valid[:, :, k] .= frame_mask
            sigmas[k] = frame.sigma * scale
        end
    else
        # WCSTransform never crosses the wire — see the docstring above.
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
    # Reusable buffer, not a fresh array per pixel — see the docstring above.
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
