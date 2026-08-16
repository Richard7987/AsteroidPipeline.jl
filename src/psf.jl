"""
    estimate_psf(image; stamp_size=25, threshold=20.0, min_separation=40.0,
                 saturate=Inf, fallback::Bool=true) -> Matrix{Float64}

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
family per instrument. But a field sparser or more crowded than expected
can leave zero stars passing the isolation/saturation filter above; when
that happens and `fallback` is true (the default), a warning is emitted
and [`fit_moffat_psf`](@ref) is used instead — an analytic fit trades the
real PSF's exact shape for something usable at all. `fallback=false`
keeps the hard failure instead.
"""
function estimate_psf(image::AbstractMatrix{<:Real}; stamp_size::Integer=25,
                       threshold::Real=20.0, min_separation::Real=40.0,
                       saturate::Real=Inf, fallback::Bool=true)
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

    if isempty(stamps)
        if fallback
            @warn "no isolated, unsaturated stars found for empirical PSF estimation; " *
                  "falling back to an analytic Moffat fit" threshold min_separation
            return fit_moffat_psf(image; stamp_size, threshold, saturate)
        end
        error("no isolated, unsaturated stars found for PSF estimation " *
              "(try lowering threshold or min_separation)")
    end

    combined = dropdims(median(cat(stamps...; dims=3); dims=3); dims=3)
    return combined ./ sum(combined)
end

"""
    fit_moffat_psf(image; stamp_size=25, threshold=20.0, saturate=Inf) -> Matrix{Float64}

Estimate a point-spread function from `image` (same convention as
[`estimate_psf`](@ref)) by fitting an analytic Moffat profile to bright
stars, rather than stacking star cutouts empirically. Used as
[`estimate_psf`](@ref)'s fallback when no stars pass its stricter
isolation/saturation filter (e.g. a crowded or sparse field), so this
function deliberately reuses only the unfiltered front half of that
filter chain — [`detect_sources`](@ref) at `threshold` sigma, plus the
same border margin check — and *drops* the `min_separation` isolation
requirement, which is usually what emptied `estimate_psf`'s star list in
the first place. That's a reasonable trade for a fit but not for a
stack: a median stack of blended-star stamps bakes the neighbour
directly into the combined PSF's shape, whereas a Moffat fit has only 6
free parameters and can't absorb a companion star into them — and a bad
stamp shows up as a poor fit residual rather than silently corrupting the
result.

For each candidate star, a circular Moffat profile
`I(r) = A * (1 + (r/α)^2)^(-β) + B` — the standard atmospheric-seeing PSF
model, with heavier wings than a Gaussian, which ZOGY's matched filter
depends on — is fit to its `stamp_size` cutout via `LsqFit.curve_fit`,
with free amplitude `A`, centre `(x0, y0)`, core width `α`, power-law
index `β`, and background offset `B`. The **median** `α` and `β` across
all stars whose fit converges (resistant to any one star's bad fit, the
same reasoning behind `estimate_psf`'s own use of a median stack)
parameterize the returned PSF: that Moffat shape evaluated on a
`stamp_size` grid centred at zero offset, normalized to unit sum — the
same convention `estimate_psf` returns, so this drops directly into
[`zogy_subtract`](@ref)'s `psf_n`/`psf_r` wherever the empirical estimate
would.

Errors if no candidate star is found at all, or if every fit fails to
converge — no usable star means no PSF, and that stays a hard failure
even here.
"""
function fit_moffat_psf(image::AbstractMatrix{<:Real}; stamp_size::Integer=25,
                         threshold::Real=20.0, saturate::Real=Inf)
    isodd(stamp_size) || throw(ArgumentError("stamp_size must be odd"))
    half = stamp_size ÷ 2

    sources = detect_sources(permutedims(image); threshold=threshold)
    nx, ny = size(image)

    xs = Float64.(-half:half)
    grid_x = repeat(xs, 1, stamp_size)
    grid_y = repeat(xs', stamp_size, 1)
    xy = hcat(vec(grid_x), vec(grid_y))

    moffat(xy, p) = @. p[1] / (1 + ((xy[:, 1] - p[2])^2 + (xy[:, 2] - p[3])^2) / p[4]^2)^p[5] + p[6]

    alphas = Float64[]
    betas = Float64[]

    for s in sources
        x, y = round(Int, s.x), round(Int, s.y)
        (half < x <= nx - half && half < y <= ny - half) || continue

        stamp = image[x-half:x+half, y-half:y+half]
        maximum(stamp) < saturate || continue

        background = median(vcat(stamp[1, :], stamp[end, :], stamp[:, 1], stamp[:, end]))
        stamp = stamp .- background

        peak = maximum(stamp)
        peak > 0 || continue

        p0 = [peak, 0.0, 0.0, 3.0, 2.5, 0.0]
        fit = try
            curve_fit(moffat, xy, vec(stamp), p0)
        catch
            continue
        end
        fit.converged || continue

        alpha, beta = fit.param[4], fit.param[5]
        alpha > 0 && beta > 0 || continue
        push!(alphas, alpha)
        push!(betas, beta)
    end

    isempty(alphas) && error("no star could be fit for an analytic PSF model " *
                              "(try lowering threshold or raising saturate)")

    alpha, beta = median(alphas), median(betas)
    profile = @. 1 / (1 + (grid_x^2 + grid_y^2) / alpha^2)^beta
    return profile ./ sum(profile)
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
