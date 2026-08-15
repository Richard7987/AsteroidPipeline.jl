"""
    zogy_subtract(n_image, r_image; psf_n, psf_r, sigma_n, sigma_r,
                  gain_n=1.0, gain_r=1.0, f_n=1.0, f_r=1.0) -> (D, S_corr)

Proper image subtraction (Zackay, Ofek & Gal-Yam 2016, ApJ 830, 27) between
a science frame `n_image` ("new") and a reference frame `r_image`, both raw
FITS pixel order on the same grid (e.g. `n_image` from [`load_frame`](@ref),
`r_image` from [`build_reference`](@ref)) with matching PSFs `psf_n`,
`psf_r` (see [`estimate_psf`](@ref)) and background RMS `sigma_n`,
`sigma_r`.

Computed entirely in Fourier space (`FFTW.fft`), all as `Δ̂ = σ_n²F_r²|P̂_r|²
+ σ_r²F_n²|P̂_n|²`:

- `D̂ = (F_r P̂_r N̂ − F_n P̂_n R̂) / √Δ̂` — the difference image proper: PSF
  matching and subtraction combined into one deconvolution, avoiding the
  correlated-noise artifacts of "convolve-then-subtract" methods.
- `Ŝ = F_D conj(P̂_D) D̂` with `P̂_D = F_r F_n P̂_r P̂_n / (F_D √Δ̂)` and
  `F_D = F_n F_r / √(σ_n²F_r² + σ_r²F_n²)` — `D` matched-filtered by its own
  PSF, which is what makes a pixel's value in `S` proportional to the flux
  of a point source centred there.
- `S_corr = S / √(V(S_N) + V(S_R) + V_ast(S_N) + V_ast(S_R))` — `S`
  normalized by its own propagated variance (photon noise from both images
  plus, when `n_sources`/`r_sources` are given, sub-pixel astrometric
  registration error), so that under the null hypothesis (no source)
  `S_corr` has zero mean and **unit variance everywhere** and its value at
  any pixel is directly a detection significance in sigma. This is what
  `detect_sources` runs on in the difference-imaging path of
  [`run_pipeline`](@ref).

`n_image` and `r_image` need not be background-subtracted beforehand —
each is subtracted its own median internally, before anything else. This
matters more than it looks: a raw FFT's DC term is the *sum*, not the
mean, of a whole image (often >10^6 pixels here), so even a few ADU of
unmatched sky background between `n_image` and `r_image` — routine, since
sky brightness depends on moon phase and airglow, not the photometric
zeropoint `build_reference` already reconciles — blows up into a
near-constant offset across all of `D` and `S`, swamping `S_corr` almost
everywhere. This was caught by testing against real ZTF data, not by the
synthetic tests, all of which happened to use matched backgrounds.

`f_n`, `f_r` are relative photometric flux scales; left at 1 because
`build_reference` already rescales every reference frame to the science
frame's zeropoint before stacking, so both inputs are on a common scale by
the time they reach here. `gain_n`, `gain_r` (e-/ADU, from each frame's
`GAIN` header) convert background-subtracted counts to Poisson variance for
the `V_N`, `V_R` terms.

Astrometric registration noise (`V_ast`) is included only if `n_sources`
and `r_sources` (tables with `x`, `y` columns, e.g. from
[`detect_sources`](@ref) on `n_image`/`r_image` directly, *not* on `D` or
`S`) are passed as keywords; otherwise it is omitted, which is only valid
if reprojection onto a common grid has already removed essentially all
registration error — reasonable here since `build_reference` already
reprojects reference frames onto the exact science-frame grid.
"""
function zogy_subtract(n_image::AbstractMatrix{<:Real}, r_image::AbstractMatrix{<:Real};
                        psf_n::AbstractMatrix{<:Real}, psf_r::AbstractMatrix{<:Real},
                        sigma_n::Real, sigma_r::Real,
                        gain_n::Real=1.0, gain_r::Real=1.0,
                        f_n::Real=1.0, f_r::Real=1.0,
                        n_sources=nothing, r_sources=nothing)
    size(n_image) == size(r_image) || throw(ArgumentError("n_image and r_image must be the same size"))
    shape = size(n_image)

    n_bg, r_bg = median(n_image), median(r_image)
    n_sub, r_sub = Float64.(n_image) .- n_bg, Float64.(r_image) .- r_bg

    N = fft(n_sub)
    R = fft(r_sub)
    P_n = fft(_pad_kernel(psf_n, shape))
    P_r = fft(_pad_kernel(psf_r, shape))

    Δ = sigma_n^2 * f_r^2 * abs2.(P_r) .+ sigma_r^2 * f_n^2 * abs2.(P_n)
    Δ_safe = max.(Δ, eps(Float64) * maximum(Δ))

    D_hat = (f_r .* P_r .* N .- f_n .* P_n .* R) ./ sqrt.(Δ_safe)
    D = real.(ifft(D_hat))

    F_D = f_n * f_r / sqrt(sigma_n^2 * f_r^2 + sigma_r^2 * f_n^2)
    P_D_hat = f_r * f_n .* P_r .* P_n ./ (F_D .* sqrt.(Δ_safe))
    S_hat = F_D .* conj.(P_D_hat) .* D_hat
    S = real.(ifft(S_hat))

    # S = k_n ⊛ N - k_r ⊛ R by definition; expanding S = F_D conj(P_D) D
    # (D and P_D above) gives these directly — note F_D cancels entirely.
    k_n_hat = f_r^2 * f_n .* conj.(P_n) .* abs2.(P_r) ./ Δ_safe
    k_r_hat = f_r * f_n^2 .* conj.(P_r) .* abs2.(P_n) ./ Δ_safe
    k_n = real.(ifft(k_n_hat))
    k_r = real.(ifft(k_r_hat))

    V_n = _variance_map(n_sub, sigma_n, gain_n)
    V_r = _variance_map(r_sub, sigma_r, gain_r)
    var_S = _convolve_variance(V_n, k_n) .+ _convolve_variance(V_r, k_r)

    if n_sources !== nothing && r_sources !== nothing
        sigma_x, sigma_y = _astrometric_scatter(n_sources, r_sources)
        if sigma_x > 0 || sigma_y > 0
            dSdx, dSdy = _gradient(S)
            var_S = var_S .+ sigma_x^2 .* dSdx .^ 2 .+ sigma_y^2 .* dSdy .^ 2
        end
    end

    S_corr = S ./ sqrt.(max.(var_S, eps(Float64) * maximum(var_S)))
    return D, S_corr
end

"""
    _pad_kernel(kernel, shape) -> Matrix{Float64}

Zero-pad `kernel` (odd-sized, peak at its centre) to `shape` with its peak
wrapped to array index `(1, 1)`, the layout `fft`/`ifft` require for
convolution: place it centred, then `ifftshift`.
"""
function _pad_kernel(kernel::AbstractMatrix{<:Real}, shape::NTuple{2,Integer})
    all(isodd, size(kernel)) || throw(ArgumentError("kernel must have odd dimensions"))
    padded = zeros(Float64, shape)
    hy, hx = size(kernel) .÷ 2
    cy, cx = shape .÷ 2 .+ 1
    padded[cy-hy:cy+hy, cx-hx:cx+hx] .= kernel
    return ifftshift(padded)
end

"""
    _variance_map(image_sub, sigma, gain) -> Matrix{Float64}

Per-pixel variance: background RMS plus Poisson shot noise from
above-background counts, `σ² + max(image_sub, 0) / gain`. `image_sub` must
already be background-subtracted (by the same convention `zogy_subtract`
itself uses, so the two agree on where "above background" starts).
"""
function _variance_map(image_sub::AbstractMatrix{<:Real}, sigma::Real, gain::Real)
    return sigma^2 .+ max.(image_sub, 0) ./ gain
end

"""
    _convolve_variance(V, k) -> Matrix{Float64}

`k² ⊗ V`: how a linear filter `k` propagates per-pixel variance `V`. For
independent per-pixel noise, `Var[sum_y k(x-y) D(y)] = sum_y k(x-y)² V(y)`,
i.e. a convolution of `V` with `k` *squared in real space* — not with `|k̂|²`,
a different quantity — computed here via the convolution theorem applied to
that squared kernel.
"""
function _convolve_variance(V::AbstractMatrix{<:Real}, k::AbstractMatrix{<:Real})
    return real.(ifft(fft(V) .* fft(k .^ 2)))
end

"""
    _astrometric_scatter(n_sources, r_sources; match_radius=3.0) -> (sigma_x, sigma_y)

Standard deviation of the pixel-position offsets between sources in
`n_sources` and their nearest counterpart in `r_sources` within
`match_radius` pixels — an empirical estimate of residual registration
error between the two images. Returns `(0.0, 0.0)` if fewer than 5 pairs
match, too few to trust a scatter estimate.
"""
function _astrometric_scatter(n_sources, r_sources; match_radius::Real=3.0)
    dxs = Float64[]
    dys = Float64[]
    for s in n_sources
        best_d, best = Inf, nothing
        for r in r_sources
            d = hypot(s.x - r.x, s.y - r.y)
            if d < best_d
                best_d, best = d, r
            end
        end
        if best !== nothing && best_d <= match_radius
            push!(dxs, s.x - best.x)
            push!(dys, s.y - best.y)
        end
    end
    length(dxs) >= 5 || return (0.0, 0.0)
    return std(dxs), std(dys)
end

"""
    _gradient(image) -> (dIdx, dIdy)

Central-difference gradient (forward/backward at the edges), for the ZOGY
astrometric noise term, which needs `∂S/∂x` and `∂S/∂y`.
"""
function _gradient(image::AbstractMatrix{<:Real})
    nx, ny = size(image)
    dIdx = similar(image, Float64)
    dIdy = similar(image, Float64)
    for j in 1:ny, i in 1:nx
        il, ir = max(i - 1, 1), min(i + 1, nx)
        dIdx[i, j] = (image[ir, j] - image[il, j]) / (ir - il)
        jl, jr = max(j - 1, 1), min(j + 1, ny)
        dIdy[i, j] = (image[i, jr] - image[i, jl]) / (jr - jl)
    end
    return dIdx, dIdy
end
