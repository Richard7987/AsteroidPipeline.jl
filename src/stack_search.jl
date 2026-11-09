"""
    stack_search(images, wcss, epochs; masks=nothing, max_rate=100.0,
                 peak_threshold=3.0, min_score=7.0) -> (tracklets, threshold)

Find moving objects too faint to be detected in each frame on its own, by
adding their light along trial motions. `images` are a same-field sequence
(`(y, x)` matrices, as [`detect_sources`](@ref) takes) with their solutions
`wcss` and `epochs` (Julian Dates), in time order; `masks` (`true` = invalid)
as for [`stack_difference`](@ref), whose normalized differences this searches.
`max_rate` is in arcsec/hour; the slowest motions searched cross twice the
FWHM over the sequence, as slower objects subtract themselves in its median.

Each difference is filtered with a Gaussian of the frames' measured FWHM and
renormalized to unit noise. Peaks above `peak_threshold` in the first and the
last frame are paired into motions, and each frame is read at the predicted
position. The score is the faintest frame's value times `√N`: every frame
must show the object, which keeps one frame's artifact from making a
detection. On Pan-STARRS1 sets with injected asteroids it raised those
recovered from 82 to 101, all added at V 21.4-23.0.

Several tracks through one spot are dropped, all of them: a point source
makes one, a patch of extended light many.

What noise can reach is measured on the same frames: the search is repeated
with the frames out of time order, where nothing real lines up, and the
threshold is 1.1 times the best score that reaches, at least `min_score`. On
a field whose noise is not Gaussian, the threshold rises and little is found,
instead of false candidates.

The artifacts [`run_pipeline`](@ref) drops are dropped here too: peaks that
recur on one detector pixel ([`_drop_detector_fixed`](@ref)) are masked
before the search, and a track passing within 3 FWHM of a masked pixel, or
over a bright star's light ([`_drop_on_bright_background`](@ref)), is left out.
So is one whose light, within 4 FWHM, is drawn out into a line (axis ratio
above 1.7 in most frames): a bright star's diffraction spike.

Returns the tracklets as [`astrometric_calibrate`](@ref) does (`id`, `frame`,
`x`, `y`, `ra`, `dec`, `epoch`, positions refined to a fraction of a pixel on
the filtered images) and the threshold used.
"""
function stack_search(images::AbstractVector{<:AbstractMatrix{<:Real}}, wcss, epochs::AbstractVector{<:Real};
                      masks=nothing, max_rate::Real=100.0,
                      peak_threshold::Real=3.0, min_score::Real=7.0)
    n = length(images)
    masks === nothing && (masks = [falses(size(img)) for img in images])
    differences, dmasks = stack_difference(images, wcss; masks)
    fwhm = median(_source_radius(img, m, 5.0) for (img, m) in zip(images, masks))
    maps = [_matched_filter(d, m, fwhm) for (d, m) in zip(differences, dmasks)]
    dmasks = _mask_detector_fixed(maps, dmasks, fwhm, peak_threshold, epochs)

    n1, n2 = size(images[1])
    centre = [n2 / 2, n1 / 2]
    sky = pix_to_world(wcss[1], centre)
    offsets = [Tuple(round.(Int, world_to_pix(w, sky) .- centre)) for w in wcss]
    span = maximum(epochs) - minimum(epochs)
    span > 0 || return (tracklets=astrometric_calibrate([], wcss, epochs), threshold=NaN)
    τ = (epochs .- mean(epochs)) ./ span
    scale = _pixel_scale_arcsec(wcss[1], size(images[1]))
    # end-to-end displacement over the sequence, pixels
    dmin = 2fwhm
    dmax = max_rate * span * 24 / scale
    clear(c, order) = _clear_of_masks(c, order, dmasks, offsets, τ, ceil(Int, 3fwhm))

    kept_from = min_score / 1.1
    best = 0.0
    for order in _scrambled_orders(n)
        found = _pair_search(maps[order], dmasks[order], offsets[order], τ, dmin, dmax, peak_threshold, kept_from)
        for c in found
            c.score > best && clear(c, order) && (best = c.score)
        end
    end
    threshold = max(min_score, 1.1best)

    order = collect(1:n)
    r = round(Int, 4fwhm)
    tracks = filter(_pair_search(maps, dmasks, offsets, τ, dmin, dmax, peak_threshold, kept_from)) do c
        c.score >= threshold && clear(c, order) &&
            !any(k -> _on_bright_light(images[k], masks[k], _track_pixel(c, k, offsets, τ)), 1:n) &&
            median(_elongation(maps[k], dmasks[k], _track_pixel(c, k, offsets, τ)..., r) for k in 1:n) <= 1.7
    end
    shared(a, b) = any(k -> hypot((_track_pixel(a, k, offsets, τ) .- _track_pixel(b, k, offsets, τ))...) <= fwhm, 1:n)
    tracks = [a for a in tracks if !any(b -> b !== a && shared(a, b), tracks)]
    tracklets = map(tracks) do c
        map(1:n) do k
            x, y = _refined_peak(maps[k], dmasks[k], _track_pixel(c, k, offsets, τ)...)
            (frame=k, x=x, y=y)
        end
    end
    return (tracklets=astrometric_calibrate(tracklets, wcss, epochs), threshold=threshold)
end

# `S` convolved with a Gaussian of `fwhm`, in units of its own robust noise
# (masked pixels 0): a matched filter for a point source.
function _matched_filter(S::AbstractMatrix, mask::AbstractMatrix{Bool}, fwhm::Real)
    σ = fwhm / (2sqrt(2log(2)))
    r = ceil(Int, 3σ)
    kernel = [exp(-i^2 / (2σ^2)) for i in -r:r]
    A = Float32.(ifelse.(mask .| .!isfinite.(S), 0.0, S))
    n1, n2 = size(A)
    B = zeros(Float32, n1, n2)
    C = zeros(Float32, n1, n2)
    Threads.@threads for x in 1:n2
        for i in -r:r, y in max(1, 1 - i):min(n1, n1 - i)
            @inbounds B[y, x] += kernel[i+r+1] * A[y+i, x]
        end
    end
    Threads.@threads for x in 1:n2
        for i in -r:r
            1 <= x + i <= n2 || continue
            w = kernel[i+r+1]
            @inbounds for y in 1:n1
                C[y, x] += w * B[y, x+i]
            end
        end
    end
    valid = C[.!mask]
    isempty(valid) && return C
    noise = 1.4826f0 * median(abs.(valid .- median(valid)))
    noise > 0 && (C ./= noise)
    C[mask] .= 0
    return C
end

# Local maxima at or above `threshold`, (x, y), over a 5x5 neighbourhood.
function _local_peaks(M::AbstractMatrix, mask::AbstractMatrix{Bool}, threshold::Real)
    n1, n2 = size(M)
    peaks = Tuple{Int,Int}[]
    for x in 3:n2-2, y in 3:n1-2
        v = M[y, x]
        (v >= threshold && !mask[y, x]) || continue
        all(M[y+i, x+j] <= v for j in -2:2, i in -2:2) && push!(peaks, (x, y))
    end
    return peaks
end

# A hot pixel or bad column stays on one detector pixel while the telescope
# dithers: masked, as run_pipeline drops it.
function _mask_detector_fixed(maps, masks, fwhm, threshold, epochs)
    peaks = [_local_peaks(M, m, threshold) for (M, m) in zip(maps, masks)]
    tables = [Table(xcen=Float64[p[1] for p in ps], ycen=Float64[p[2] for p in ps]) for ps in peaks]
    kept = _drop_detector_fixed(tables, 2.0, epochs)
    r = ceil(Int, fwhm)
    out = copy.(masks)
    for k in eachindex(peaks)
        keep = Set(zip(kept[k].xcen, kept[k].ycen))
        n1, n2 = size(out[k])
        for (x, y) in peaks[k]
            (Float64(x), Float64(y)) in keep && continue
            out[k][max(1, y - r):min(n1, y + r), max(1, x - r):min(n2, x + r)] .= true
        end
    end
    return out
end

# The highest unmasked value within a pixel of (x, y), -Inf if none.
function _best_near(M, mask, x::Integer, y::Integer)
    n1, n2 = size(M)
    best = -Inf32
    for j in -1:1, i in -1:1
        yy, xx = y + i, x + j
        (1 <= yy <= n1 && 1 <= xx <= n2 && !mask[yy, xx]) || continue
        best = max(best, M[yy, xx])
    end
    return best
end

# Frame k's pixel on track `c` (positions at τ = 0 in frame 1's grid).
_track_pixel(c, k, offsets, τ) =
    (round(Int, c.x + c.dx * τ[k]) + offsets[k][1], round(Int, c.y + c.dy * τ[k]) + offsets[k][2])

# Tracks from each pair of peaks in the first and last frame, `dmin` to
# `dmax` pixels apart, scored in every frame; those scoring at least `kept_from`,
# one per position, the best.
function _pair_search(maps, masks, offsets, τ, dmin, dmax, threshold, kept_from)
    n = length(maps)
    first_peaks = [(x - offsets[1][1], y - offsets[1][2]) for (x, y) in _local_peaks(maps[1], masks[1], threshold)]
    last_peaks = [(x - offsets[n][1], y - offsets[n][2]) for (x, y) in _local_peaks(maps[n], masks[n], threshold)]
    T = @NamedTuple{x::Float64, y::Float64, dx::Float64, dy::Float64, score::Float64}
    chunks = Iterators.partition(first_peaks, max(1, cld(length(first_peaks), 4Threads.nthreads())))
    tasks = map(chunks) do chunk
        Threads.@spawn begin
            out = T[]
            values = zeros(Float32, n)
            for p in chunk, q in last_peaks
                dx, dy = (q[1] - p[1]) / (τ[n] - τ[1]), (q[2] - p[2]) / (τ[n] - τ[1])
                dmin <= hypot(dx, dy) <= dmax || continue
                c = (x=p[1] - dx * τ[1], y=p[2] - dy * τ[1], dx=dx, dy=dy, score=0.0)
                for k in 1:n
                    values[k] = _best_near(maps[k], masks[k], _track_pixel(c, k, offsets, τ)...)
                end
                score = minimum(values) * sqrt(n)
                score >= kept_from && push!(out, merge(c, (score=Float64(score),)))
            end
            out
        end
    end
    found = sort(reduce(vcat, fetch.(tasks); init=T[]); by=c -> -c.score)
    kept = T[]
    cells = Dict{Tuple{Int,Int},Vector{Int}}()
    for c in found
        cell = (floor(Int, c.x / 4), floor(Int, c.y / 4))
        near = any(hypot(kept[i].x - c.x, kept[i].y - c.y) <= 3
                   for di in -1:1, dj in -1:1 for i in get(cells, (cell[1] + di, cell[2] + dj), Int[]))
        near && continue
        push!(kept, c)
        push!(get!(cells, cell, Int[]), length(kept))
    end
    return kept
end

# No masked pixel within `r` of the track in any frame, nor the image's edge.
function _clear_of_masks(c, order, masks, offsets, τ, r)
    for (k, f) in enumerate(order)
        x = round(Int, c.x + c.dx * τ[k]) + offsets[f][1]
        y = round(Int, c.y + c.dy * τ[k]) + offsets[f][2]
        n1, n2 = size(masks[f])
        (1 + r <= x <= n2 - r && 1 + r <= y <= n1 - r) || return false
        any(@view masks[f][y-r:y+r, x-r:x+r]) && return false
    end
    return true
end

# Axis ratio of the positive light within `r` of (x, y), from its second
# moments: about 1.0-1.45 for injected point sources, 1.85-2.2 for the
# diffraction spikes stacking found on Pan-STARRS1 sets.
function _elongation(M, mask, x::Integer, y::Integer, r::Integer)
    n1, n2 = size(M)
    sxx = syy = sxy = total = 0.0
    for j in -r:r, i in -r:r
        hypot(i, j) <= r || continue
        yy, xx = y + i, x + j
        (1 <= yy <= n1 && 1 <= xx <= n2 && !mask[yy, xx]) || continue
        v = Float64(M[yy, xx])
        v > 0 || continue
        sxx += v * j^2; syy += v * i^2; sxy += v * i * j; total += v
    end
    total > 0 || return 1.0
    sxx /= total; syy /= total; sxy /= total
    root = sqrt(((sxx - syy) / 2)^2 + sxy^2)
    major, minor = (sxx + syy) / 2 + root, (sxx + syy) / 2 - root
    return sqrt(major / max(minor, 1e-9))
end

_on_bright_light(image, mask, (x, y)) =
    isempty(_drop_on_bright_background(Table(xcen=[Float64(x)], ycen=[Float64(y)]), image, mask, 1.0))

# Frame orders in which no constant motion lines up, for N >= 3 frames.
function _scrambled_orders(n::Integer)
    swap(start) = (o = collect(1:n); for i in start:2:n-1; o[i], o[i+1] = o[i+1], o[i]; end; o)
    return [swap(1), swap(2)]
end

# The peak nearest (x, y), to a fraction of a pixel by a parabola through it
# and its neighbours on each axis.
function _refined_peak(M, mask, x::Integer, y::Integer)
    n1, n2 = size(M)
    bx, by, best = x, y, -Inf32
    for j in -1:1, i in -1:1
        yy, xx = y + i, x + j
        (1 <= yy <= n1 && 1 <= xx <= n2 && !mask[yy, xx]) || continue
        M[yy, xx] > best && ((bx, by, best) = (xx, yy, M[yy, xx]))
    end
    function vertex(a, b, c)
        d = a - 2b + c
        return d < 0 ? clamp(0.5 * (a - c) / d, -0.5, 0.5) : 0.0
    end
    fx = 1 < bx < n2 ? vertex(M[by, bx-1], M[by, bx], M[by, bx+1]) : 0.0
    fy = 1 < by < n1 ? vertex(M[by-1, bx], M[by, bx], M[by+1, bx]) : 0.0
    return (bx + fx, by + fy)
end
