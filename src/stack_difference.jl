"""
    stack_difference(images, wcss; masks=nothing, min_valid::Integer=3,
                     bright_star_sigma::Real=50.0, bright_star_radius::Integer=4,
                     halo_factor::Real=3.0, noise_box::Integer=128)
        -> (differences, masks)

Difference each frame of a short same-field sequence against the
per-pixel median of the whole sequence — a numerical version of the
"blink" a human does in Astrometrica, for sequences with no separate
deep reference stack (an IASC image set is 4 exposures of one field,
nothing else). Anything static (stars, galaxies) is in the median and
subtracts away; a moving object is at a different position in each frame,
so the median at any one position barely contains it, and it survives in
that frame's difference only.

`images` are `(y, x)`-ordered matrices (as [`detect_sources`](@ref)
takes), `wcss` their solutions, used only to align the frames: each frame
is shifted onto the first frame's pixel grid by the whole-pixel offset
its WCS implies at the image centre (frames of one sequence dither by a
few pixels; only *relative* WCS accuracy matters here, not absolute).
`masks` (`true` = invalid, e.g. from [`fill_value_mask`](@ref)) excludes
pixels from the median; a pixel needs at least `min_valid` valid frames
(3, so the median of a 4-frame set still out-votes one frame's mover).

Before taking the median each frame is background-subtracted and scaled
to a common flux level — fitted from bright pixels, not assumed — since
transparency changes between exposures (a real IASC set had one frame
~29% fainter than the others): unscaled, every star would leave a
residual in every frame. Whole-pixel alignment and changing seeing
still leave residuals around bright stars, so each connected region
where the median exceeds `bright_star_sigma` times the difference noise
is masked out to `bright_star_radius` plus `halo_factor` times its own
equivalent radius (`sqrt(area / π)`) — a brighter star has a bigger
core above that level *and* a wider halo, so the mask has to scale with
it. Measured on the real IASC set: one saturated star's halo stayed 3σ
above the sky out to ~40 px from a ~12 px core, and a fixed 4 px margin
left its residuals linked into 14 of 16 surviving tracklets.

Each difference is divided by its own *local* noise — a robust (MAD)
estimate in `noise_box`-pixel boxes — so it is a significance map, like
ZOGY's `S_corr`: [`detect_sources`](@ref)'s `threshold` then means the
same thing everywhere. A single frame-wide noise figure doesn't: PS1's
cells differ in noise, and on a real 2019 IASC set (`XY15_p01`) a couple
of noisier cells produced ~3,000 of a frame's ~4,000 detections at 6σ,
most barely over threshold — enough to stall linking entirely.

Returns each frame's normalized difference image, in that frame's
**own** pixel grid (so its own WCS still applies to detections on it),
and each frame's mask of pixels with no valid difference.

Validated on the real IASC practice set it was built for (ps1-NewPractice_3,
field `XY54_p10`), where it recovered both real moving objects that
detecting on the raw frames missed at the same threshold: an unknown
object at G≈21 (later measured in Astrometrica and reported as NHU0001)
and the known Mars-crosser 2018 LT.
"""
function stack_difference(images::AbstractVector{<:AbstractMatrix{<:Real}}, wcss;
                          masks=nothing, min_valid::Integer=3,
                          bright_star_sigma::Real=50.0, bright_star_radius::Integer=4,
                          halo_factor::Real=3.0, noise_box::Integer=128)
    n = length(images)
    n >= min_valid || throw(ArgumentError("stack_difference needs at least min_valid=$min_valid frames, got $n"))
    length(wcss) == n || throw(ArgumentError("images and wcss must have the same length"))
    masks === nothing && (masks = [falses(size(img)) for img in images])
    n1, n2 = size(images[1])

    # whole-pixel shift of each frame relative to frame 1, at frame 1's centre
    centre = [n2 / 2, n1 / 2]
    sky = pix_to_world(wcss[1], centre)
    shifts = [round.(Int, world_to_pix(w, sky) .- centre) for w in wcss]   # (dx, dy)

    backgrounds = [median(img[.!m]) for (img, m) in zip(images, masks)]
    # frame k, background-subtracted, sampled on frame 1's grid (NaN = invalid)
    aligned = [fill(NaN, n1, n2) for _ in 1:n]
    for k in 1:n
        dx, dy = shifts[k]; img = images[k]; m = masks[k]
        nk1, nk2 = size(img)
        for j in 1:n2, i in 1:n1
            ii, jj = i + dy, j + dx
            (1 <= ii <= nk1 && 1 <= jj <= nk2 && !m[ii, jj]) || continue
            aligned[k][i, j] = img[ii, jj] - backgrounds[k]
        end
    end

    template = _nanmedian_stack(aligned, min_valid)
    # flux scale per frame, from pixels clearly on stars in the template
    tvalid = filter(!isnan, template)
    tnoise = 1.4826 * median(abs.(tvalid .- median(tvalid)))
    bright = .!isnan.(template) .& (template .> 20 * tnoise)
    scales = map(aligned) do a
        sel = bright .& .!isnan.(a)
        count(sel) < 20 ? 1.0 : median(a[sel] ./ template[sel])
    end
    template = _nanmedian_stack([a ./ s for (a, s) in zip(aligned, scales)], min_valid)

    differences = Matrix{Float64}[]
    out_masks = BitMatrix[]
    for k in 1:n
        dx, dy = shifts[k]; img = images[k]
        nk1, nk2 = size(img)
        diff = zeros(nk1, nk2)
        invalid = trues(nk1, nk2)
        for jj in 1:nk2, ii in 1:nk1
            i, j = ii - dy, jj - dx
            (1 <= i <= n1 && 1 <= j <= n2) || continue
            t = template[i, j]
            (isnan(t) || masks[k][ii, jj]) && continue
            diff[ii, jj] = img[ii, jj] - backgrounds[k] - scales[k] * t
            invalid[ii, jj] = false
        end
        dvalid = diff[.!invalid]
        dnoise = isempty(dvalid) ? 0.0 : 1.4826 * median(abs.(dvalid .- median(dvalid)))
        _mask_bright!(invalid, template, shifts[k], scales[k] * bright_star_sigma * dnoise,
                      bright_star_radius, halo_factor)
        diff[invalid] .= 0.0
        push!(differences, diff ./ _local_noise(diff, invalid, noise_box))
        push!(out_masks, invalid)
    end
    return differences, out_masks
end

function _nanmedian_stack(stack, min_valid::Integer)
    n1, n2 = size(stack[1])
    out = fill(NaN, n1, n2)
    buf = Float64[]
    for j in 1:n2, i in 1:n1
        empty!(buf)
        for a in stack
            v = a[i, j]
            isnan(v) || push!(buf, v)
        end
        length(buf) >= min_valid && (out[i, j] = median!(buf))
    end
    return out
end

# Mask (in frame k's own grid) a disk around every connected region of
# `template` brighter than `level`: radius `radius` plus `halo_factor`
# times the region's equivalent radius (see `stack_difference`).
function _mask_bright!(invalid, template, shift, level::Real, radius::Integer, halo_factor::Real)
    level > 0 || return invalid
    dx, dy = shift
    nk1, nk2 = size(invalid)
    for (ci, cj, area) in _bright_regions(template, level)
        r = radius + halo_factor * sqrt(area / π)
        ii, jj = ci + dy, cj + dx
        for j in max(1, floor(Int, jj - r)):min(nk2, ceil(Int, jj + r)),
            i in max(1, floor(Int, ii - r)):min(nk1, ceil(Int, ii + r))
            hypot(i - ii, j - jj) <= r && (invalid[i, j] = true)
        end
    end
    return invalid
end

# Connected (4-neighbour) regions of `template` above `level`, as
# (centroid row, centroid column, area in pixels). NaN never counts as
# above, so a saturated core masked out of the template leaves a ring —
# its area still grows with the star's brightness, which is what matters.
function _bright_regions(template::AbstractMatrix, level::Real)
    n1, n2 = size(template)
    seen = falses(n1, n2)
    regions = Tuple{Float64,Float64,Int}[]
    stack = Tuple{Int,Int}[]
    for j in 1:n2, i in 1:n1
        (seen[i, j] || !(template[i, j] > level)) && continue
        seen[i, j] = true
        push!(stack, (i, j))
        si = sj = 0.0
        area = 0
        while !isempty(stack)
            a, b = pop!(stack)
            si += a; sj += b; area += 1
            for (p, q) in ((a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1))
                (1 <= p <= n1 && 1 <= q <= n2 && !seen[p, q] && template[p, q] > level) || continue
                seen[p, q] = true
                push!(stack, (p, q))
            end
        end
        push!(regions, (si / area, sj / area, area))
    end
    return regions
end

# Robust per-box noise of `diff` over its valid pixels, as a full-size map
# (each pixel takes its box's value; boxes with too few valid pixels take
# the whole frame's).
function _local_noise(diff::AbstractMatrix, invalid::AbstractMatrix{Bool}, box::Integer)
    n1, n2 = size(diff)
    mad_sigma(v) = 1.4826 * median(abs.(v .- median(v)))
    valid_all = diff[.!invalid]
    fallback = isempty(valid_all) ? 1.0 : mad_sigma(valid_all)
    fallback > 0 || (fallback = 1.0)
    noise = fill(fallback, n1, n2)
    for j0 in 1:box:n2, i0 in 1:box:n1
        is, js = i0:min(n1, i0 + box - 1), j0:min(n2, j0 + box - 1)
        v = [diff[i, j] for i in is, j in js if !invalid[i, j]]
        length(v) >= box^2 ÷ 4 || continue
        # Floored at half the frame's own noise: a box whose valid pixels
        # are nearly constant (seen for real: one 2019 IASC frame, `XY25_p10`)
        # has a near-zero MAD, and dividing by it blew that frame's whole
        # normalized difference up to ~1e15, burying a real asteroid at 30σ.
        noise[is, js] .= max(mad_sigma(v), fallback / 2)
    end
    return noise
end

"""
    _drop_clustered(detections, radius) -> table

`detections` minus every detection with two or more others within
`radius` pixels of it in the same frame. On a difference image a real
moving object stands alone; what clusters is extended junk — a
satellite trail fragmenting into a row of peaks, the residual ring of
a saturated star — whose pieces [`link_candidates`](@ref) can otherwise
string together with noise from other frames. Measured on a real 2019
IASC set (`XY25_p10`): two broad satellite trails, in two of the four
frames, turned into ~240 spurious 3-frame tracklets; on a typical
difference image (~400 detections over ~6M pixels) a lone object
expects ~0.05 chance neighbours within 15 px, so this costs real
objects essentially nothing.
"""
function _drop_clustered(detections, radius::Real)
    xy = [_position(d) for d in detections]
    keep = map(eachindex(xy)) do i
        count(j -> j != i && hypot(xy[j][1] - xy[i][1], xy[j][2] - xy[i][2]) <= radius, eachindex(xy)) < 2
    end
    return detections[keep]
end
