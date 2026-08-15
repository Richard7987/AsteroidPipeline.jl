"""
    estimate_psf(image; stamp_size=25, threshold=20.0, min_separation=40.0,
                 saturate=Inf) -> Matrix{Float64}

Estimate an empirical point-spread function from `image` (raw FITS pixel
order, matching [`load_frame`](@ref)) by stacking bright, isolated stars.

Stars are found with [`detect_sources`](@ref) at `threshold` sigma (well
above the detection threshold used for the science search itself, since
these must be high-S/N enough to define the PSF cleanly), then filtered to
ones with no neighbour within `min_separation` pixels and no pixel above
`saturate` in their stamp — a blended or saturated star would bias the
combined PSF's shape. Each surviving stamp is background-subtracted, its
sub-pixel centroid found from the flux-weighted first moment, resampled
onto that centre with bilinear interpolation so misaligned integer-pixel
cutouts don't smear the stack, and normalized to unit sum. The final PSF
is the pixel-wise **median** across stamps (again to resist any residual
blend or cosmic ray a single star's stamp might carry), renormalized to
unit sum.

An empirical PSF is used rather than an analytic model (e.g. Gaussian or
Moffat) because ZOGY needs the real PSF shape, wings included, and this
keeps `estimate_psf` usable on any survey's data without fitting a model
family per instrument.
"""
function estimate_psf(image::AbstractMatrix{<:Real}; stamp_size::Integer=25,
                       threshold::Real=20.0, min_separation::Real=40.0,
                       saturate::Real=Inf)
    isodd(stamp_size) || throw(ArgumentError("stamp_size must be odd"))
    half = stamp_size ÷ 2

    sources = detect_sources(permutedims(image); threshold=threshold)
    nx, ny = size(image)

    stamps = Matrix{Float64}[]
    for s in sources
        x, y = round(Int, s.x), round(Int, s.y)
        (half < x <= nx - half && half < y <= ny - half) || continue

        isolated = all(hypot(s.x - o.x, s.y - o.y) >= min_separation
                        for o in sources if !(o.x == s.x && o.y == s.y))
        isolated || continue

        stamp = image[x-half:x+half, y-half:y+half]
        maximum(stamp) < saturate || continue

        background = median(vcat(stamp[1, :], stamp[end, :], stamp[:, 1], stamp[:, end]))
        stamp = stamp .- background

        total = sum(stamp)
        total > 0 || continue
        xs = Float64.(-half:half)
        cx = sum(xs .* sum(stamp, dims=2)[:]) / total
        cy = sum(xs .* sum(stamp, dims=1)[:]) / total

        centered = _shift_bilinear(stamp, -cx, -cy)
        s2 = sum(centered)
        s2 > 0 || continue
        push!(stamps, centered ./ s2)
    end

    isempty(stamps) && error("no isolated, unsaturated stars found for PSF estimation " *
                              "(try lowering threshold or min_separation)")

    combined = dropdims(median(cat(stamps...; dims=3); dims=3); dims=3)
    return combined ./ sum(combined)
end

"""
    _shift_bilinear(stamp, dx, dy) -> Matrix{Float64}

Resample `stamp` shifted by `(dx, dy)` pixels via bilinear interpolation,
extrapolating with 0 beyond the stamp edge. Used to sub-pixel-align star
stamps onto a common centroid before combining.
"""
function _shift_bilinear(stamp::AbstractMatrix{<:Real}, dx::Real, dy::Real)
    itp = extrapolate(interpolate(Float64.(stamp), BSpline(Linear())), 0.0)
    nx, ny = size(stamp)
    return [itp(i - dx, j - dy) for i in 1:nx, j in 1:ny]
end
