#=
Preliminary orbits from astrometry spread over nights: Gauss's method for a
first orbit, then least squares on the heliocentric state, two-body.
=#

# Closer than this (AU, ~1.5 million km) a solution is taken for the observer's own orbit.
const _MIN_DISTANCE = 0.01

# k², the Sun's gravitational parameter in AU³/day² (Gaussian constant k = 0.01720209895)
const _GM_SUN = 0.01720209895^2

# TAI − UTC (s) from each date (JD UTC) on, IERS Bulletin C, 1972 to 2017.
const _LEAP_SECONDS = [(2441317.5, 10), (2441499.5, 11), (2441683.5, 12), (2442048.5, 13), (2442413.5, 14),
                       (2442778.5, 15), (2443144.5, 16), (2443509.5, 17), (2443874.5, 18), (2444239.5, 19),
                       (2444786.5, 20), (2445151.5, 21), (2445516.5, 22), (2446247.5, 23), (2447161.5, 24),
                       (2447892.5, 25), (2448257.5, 26), (2448804.5, 27), (2449169.5, 28), (2449534.5, 29),
                       (2450083.5, 30), (2450630.5, 31), (2451179.5, 32), (2453736.5, 33), (2454832.5, 34),
                       (2456109.5, 35), (2457204.5, 36), (2457754.5, 37)]

# TT − UTC in days at `jd_utc`.
function _tt_minus_utc(jd_utc::Real)
    jd_utc >= _LEAP_SECONDS[1][1] || throw(ArgumentError("epochs before 1972 need ΔT, not leap seconds"))
    i = findlast(l -> l[1] <= jd_utc, _LEAP_SECONDS)
    return (_LEAP_SECONDS[i][2] + 32.184) / 86400
end

# Stumpff functions C(z), S(z), by series near z = 0.
function _stumpff(z::Real)
    if abs(z) < 1e-3
        return (1 / 2 - z / 24 + z^2 / 720 - z^3 / 40320, 1 / 6 - z / 120 + z^2 / 5040 - z^3 / 362880)
    elseif z > 0
        s = sqrt(z)
        return ((1 - cos(s)) / z, (s - sin(s)) / s^3)
    else
        s = sqrt(-z)
        return ((cosh(s) - 1) / -z, (sinh(s) - s) / s^3)
    end
end

"""
    _kepler_propagate(r0, v0, dt) -> (r, v)

Two-body heliocentric position (AU) and velocity (AU/day) `dt` days after
(`r0`, `v0`), by the universal variable (Curtis 2014, *Orbital Mechanics for
Engineering Students*, Algorithms 3.3-3.4): any conic.
"""
function _kepler_propagate(r0::AbstractVector, v0::AbstractVector, dt::Real)
    μ = _GM_SUN
    sμ = sqrt(μ)
    rn = norm(r0)
    vr = dot(r0, v0) / rn
    α = 2 / rn - dot(v0, v0) / μ
    χ = sμ * abs(α) * dt
    χ == 0 && (χ = sμ * dt / rn)
    for _ in 1:60
        z = α * χ^2
        C, S = _stumpff(z)
        F = rn * vr / sμ * χ^2 * C + (1 - α * rn) * χ^3 * S + rn * χ - sμ * dt
        dF = rn * vr / sμ * χ * (1 - z * S) + (1 - α * rn) * χ^2 * C + rn
        step = F / dF
        χ -= step
        abs(step) < 1e-13 * max(1.0, abs(χ)) && break
    end
    z = α * χ^2
    C, S = _stumpff(z)
    f = 1 - χ^2 / rn * C
    g = dt - χ^3 * S / sμ
    r = f .* r0 .+ g .* v0
    r1 = norm(r)
    fdot = sμ / (r1 * rn) * (z * S - 1) * χ
    gdot = 1 - χ^2 / r1 * C
    return (r, fdot .* r0 .+ gdot .* v0)
end

# Unit vector toward (ra, dec), degrees, equatorial.
_line_of_sight(ra, dec) = [cosd(dec) * cosd(ra), cosd(dec) * sind(ra), sind(dec)]

"""
    _gauss_orbits(ρ̂, R, t) -> Vector of (r2, v2)

Gauss's method (Curtis 2014, Algorithm 5.5) for three lines of sight `ρ̂`
from heliocentric observer positions `R` at TT Julian Dates `t`: one state
at the middle epoch per positive root of the eighth-degree equation for the
heliocentric distance, at least `_MIN_DISTANCE` from every observer.
"""
function _gauss_orbits(ρ̂, R, t)
    μ = _GM_SUN
    τ1, τ3 = t[1] - t[2], t[3] - t[2]
    τ = τ3 - τ1
    p1, p2, p3 = cross(ρ̂[2], ρ̂[3]), cross(ρ̂[1], ρ̂[3]), cross(ρ̂[1], ρ̂[2])
    D0 = dot(ρ̂[1], p1)
    abs(D0) < 1e-12 && return Tuple{Vector{Float64},Vector{Float64}}[]
    D = [dot(R[i], p) for i in 1:3, p in (p1, p2, p3)]
    A = (-D[1, 2] * τ3 / τ + D[2, 2] + D[3, 2] * τ1 / τ) / D0
    B = (D[1, 2] * (τ3^2 - τ^2) * τ3 / τ + D[3, 2] * (τ^2 - τ1^2) * τ1 / τ) / (6D0)
    E = dot(R[2], ρ̂[2])
    a = -(A^2 + 2A * E + dot(R[2], R[2]))
    b = -2μ * B * (A + E)
    c = -μ^2 * B^2
    # r⁸ + a r⁶ + b r³ + c = 0: its positive roots, bracketed on a log grid
    poly(r) = r^8 + a * r^6 + b * r^3 + c
    grid = exp10.(range(-3, 3; length=4000))
    roots = Float64[]
    for i in 1:length(grid)-1
        lo, hi = grid[i], grid[i+1]
        sign(poly(lo)) == sign(poly(hi)) && continue
        for _ in 1:100
            mid = (lo + hi) / 2
            sign(poly(mid)) == sign(poly(lo)) ? (lo = mid) : (hi = mid)
        end
        push!(roots, (lo + hi) / 2)
    end
    states = Tuple{Vector{Float64},Vector{Float64}}[]
    for r2 in roots
        ρ1 = ((6 * (D[3, 1] * τ1 / τ3 + D[2, 1] * τ / τ3) * r2^3 + μ * D[3, 1] * (τ^2 - τ1^2) * τ1 / τ3) /
              (6r2^3 + μ * (τ^2 - τ3^2)) - D[1, 1]) / D0
        ρ2 = A + μ * B / r2^3
        ρ3 = ((6 * (D[1, 3] * τ3 / τ1 - D[2, 3] * τ / τ1) * r2^3 + μ * D[1, 3] * (τ^2 - τ3^2) * τ3 / τ1) /
              (6r2^3 + μ * (τ^2 - τ1^2)) - D[3, 3]) / D0
        # the root at ρ ≈ 0 is the observer's own orbit
        min(ρ1, ρ2, ρ3) > _MIN_DISTANCE || continue
        r1v, r2v, r3v = R[1] .+ ρ1 .* ρ̂[1], R[2] .+ ρ2 .* ρ̂[2], R[3] .+ ρ3 .* ρ̂[3]
        f1, f3 = 1 - μ * τ1^2 / (2r2^3), 1 - μ * τ3^2 / (2r2^3)
        g1, g3 = τ1 - μ * τ1^3 / (6r2^3), τ3 - μ * τ3^3 / (6r2^3)
        push!(states, (r2v, (-f3 .* r1v .+ f1 .* r3v) ./ (f1 * g3 - f3 * g1)))
    end
    return states
end

"""
    _lambert(r1, r2, dt) -> v1, or `nothing`

The velocity at `r1` of the two-body arc reaching `r2` `dt` days later,
prograde with respect to the ecliptic, less than one revolution (Curtis 2014,
Algorithm 5.2, universal variable; `z` by bisection, on which the time of
flight increases monotonically).
"""
function _lambert(r1::AbstractVector, r2::AbstractVector, dt::Real)
    μ = _GM_SUN
    n1, n2 = norm(r1), norm(r2)
    ε = deg2rad(_OBLIQUITY_J2000)
    pole = [0.0, -sin(ε), cos(ε)]
    θ = acos(clamp(dot(r1, r2) / (n1 * n2), -1, 1))
    dot(cross(r1, r2), pole) < 0 && (θ = 2π - θ)
    A = sin(θ) * sqrt(n1 * n2 / (1 - cos(θ)))
    y(z) = (cs = _stumpff(z); n1 + n2 + A * (z * cs[2] - 1) / sqrt(cs[1]))
    function flight(z)
        yz = y(z)
        yz < 0 && return -Inf
        C, S = _stumpff(z)
        return ((yz / C)^1.5 * S + A * sqrt(yz)) / sqrt(μ)
    end
    lo, hi = -400.0, 4π^2 - 1e-9
    (flight(lo) < dt < flight(hi)) || return nothing
    for _ in 1:200
        mid = (lo + hi) / 2
        flight(mid) < dt ? (lo = mid) : (hi = mid)
    end
    yz = y((lo + hi) / 2)
    f, g = 1 - yz / n1, A * sqrt(yz / μ)
    return (r2 .- f .* r1) ./ g
end

# First states from a grid of topocentric distances at the first and last
# observation (systematic ranging: Virtanen, Muinonen & Bowell 2001), the
# arc between them by Lambert's problem; the `keep` best by weighted
# residuals, as states at `t0`.
function _ranging_orbits(ra, dec, t, R, σ, t0; keep::Integer=5,
                         distances=exp10.(range(log10(2_MIN_DISTANCE), log10(20.0); length=40)))
    n = length(ra)
    ρ̂1, ρ̂n = _line_of_sight(ra[1], dec[1]), _line_of_sight(ra[n], dec[n])
    w = repeat(1 ./ σ .^ 2; inner=2)
    tried = Tuple{Float64,Vector{Float64}}[]
    for ρ1 in distances, ρn in distances
        r1, rn = R[1] .+ ρ1 .* ρ̂1, R[n] .+ ρn .* ρ̂n
        v1 = _lambert(r1, rn, t[n] - t[1])
        v1 === nothing && continue
        r0, v0 = _kepler_propagate(r1, v1, t0 - t[1])
        x = vcat(r0, v0)
        c = sum(w .* _orbit_residuals(x, t0, ra, dec, R, t) .^ 2)
        isfinite(c) && push!(tried, (c, x))
    end
    sort!(tried; by=first)
    return [x for (_, x) in tried[1:min(keep, length(tried))]]
end

# Astrometric (RA, Dec) in degrees of the state (r0, v0) at t0 seen from R at
# t (TT), the light time solved by iteration.
function _orbit_radec(r0, v0, t0, R, t)
    Δ = 0.0
    p = r0
    for _ in 1:3
        p, _ = _kepler_propagate(r0, v0, t - Δ / _AU_PER_DAY_LIGHT - t0)
        Δ = norm(p .- R)
    end
    v = p .- R
    return (mod(atand(v[2], v[1]), 360), asind(v[3] / norm(v)))
end

# Residuals (arcsec), observed minus computed: (Δα cos δ, Δδ) for each observation.
function _orbit_residuals(x, t0, ra, dec, R, t)
    out = Vector{Float64}(undef, 2length(ra))
    for i in eachindex(ra)
        α, δ = _orbit_radec(x[1:3], x[4:6], t0, R[i], t[i])
        out[2i-1] = 3600 * mod(ra[i] - α + 180, 360) - 3600 * 180
        out[2i-1] *= cosd(dec[i])
        out[2i] = 3600 * (dec[i] - δ)
    end
    return out
end

function _orbit_jacobian(x, t0, ra, dec, R, t)
    J = Matrix{Float64}(undef, 2length(ra), 6)
    for j in 1:6
        h = j <= 3 ? 1e-7 * max(1.0, abs(x[j])) : 1e-9 * max(0.01, abs(x[j]))
        xp, xm = copy(x), copy(x)
        xp[j] += h; xm[j] -= h
        # residuals are observed - computed: the model's derivative is minus theirs
        J[:, j] = -(_orbit_residuals(xp, t0, ra, dec, R, t) .- _orbit_residuals(xm, t0, ra, dec, R, t)) ./ 2h
    end
    return J
end

# Levenberg-Marquardt on the state, observations weighted by 1/σ²; `use`
# selects the observations fitted.
function _orbit_least_squares(x, t0, ra, dec, R, t, σ, use)
    w = repeat(1 ./ σ .^ 2; inner=2) .* repeat(Float64.(use); inner=2)
    cost(x) = sum(w .* _orbit_residuals(x, t0, ra, dec, R, t) .^ 2)
    λ = 1e-3
    c = cost(x)
    for _ in 1:100
        res = _orbit_residuals(x, t0, ra, dec, R, t)
        J = _orbit_jacobian(x, t0, ra, dec, R, t)
        N = J' * (w .* J)
        g = J' * (w .* res)
        # a start far off can make the derivatives overflow
        all(isfinite, N) && all(isfinite, g) || break
        improved = false
        for _ in 1:20
            step = (N + λ * Diagonal(diag(N))) \ g
            xn = x .+ step
            cn = cost(xn)
            if isfinite(cn) && cn < c
                converged = c - cn < 1e-10 * c
                x, c, λ, improved = xn, cn, λ / 10, true
                converged && return (x, c)
                break
            end
            λ *= 10
        end
        improved || break
    end
    return (x, c)
end

# Osculating heliocentric ecliptic J2000 elements of a state (AU, AU/day, equatorial).
function _state_to_elements(r_eq, v_eq)
    μ = _GM_SUN
    ε = deg2rad(_OBLIQUITY_J2000)
    rot(p) = [p[1], cos(ε) * p[2] + sin(ε) * p[3], -sin(ε) * p[2] + cos(ε) * p[3]]
    r, v = rot(r_eq), rot(v_eq)
    rn = norm(r)
    h = cross(r, v)
    i = acosd(clamp(h[3] / norm(h), -1, 1))
    node = mod(atand(h[1], -h[2]), 360)
    evec = cross(v, h) ./ μ .- r ./ rn
    e = norm(evec)
    a = 1 / (2 / rn - dot(v, v) / μ)
    nvec = [-h[2], h[1], 0.0]
    nn = norm(nvec)
    peri = nn == 0 || e == 0 ? 0.0 : mod(acosd(clamp(dot(nvec, evec) / (nn * e), -1, 1)) * (evec[3] < 0 ? -1 : 1), 360)
    ν = e == 0 ? 0.0 : acos(clamp(dot(evec, r) / (e * rn), -1, 1)) * (dot(r, v) < 0 ? -1 : 1)
    M = if e < 1
        E = 2atan(sqrt((1 - e) / (1 + e)) * tan(ν / 2))
        mod(rad2deg(E - e * sin(E)), 360)
    else
        NaN
    end
    return (a=a, e=e, i=i, node=node, peri=peri, M=M, q=a * (1 - e))
end

"""
    fit_orbit(ra, dec, epochs, observatory; sigma=0.5, outlier_sigma=3.0)

A preliminary orbit from astrometry: `ra`, `dec` (degrees) at `epochs`
(Julian Dates UTC, 1972 or later), each from the MPC site `observatory` (one
code, or one per observation). The observers' heliocentric positions come
from JPL Horizons (a live request).

First orbits come from Gauss's method on the first, middle and last
observation (one per root of its eighth-degree equation) and from systematic
ranging (Virtanen, Muinonen & Bowell 2001): a grid of distances at the first
and last observation, joined by Lambert's problem, of which the five fitting
best are kept. Each is refined
by least squares on the heliocentric state at the middle epoch, with
observations weighted by `sigma` (arcsec, one or per observation), and the
one with the smallest weighted residuals kept, among those not within
0.01 AU of the observer (Gauss's spurious root, the observer's own orbit). Observations whose residual
exceeds `outlier_sigma` times their `sigma` are then left out and the fit
repeated, until none is (Carpino, Milani & Chesley 2003 reject on the same
principle, with a χ² threshold).

The motion is two-body around the Sun, with light time: no planet perturbs
it, so the orbit predicts well over weeks to months for a main-belt object,
less for one passing near a planet. Observations from one night alone don't
constrain it; the arc should cover several nights.

Returns a named tuple: `epoch` (JD TT of the state), `state` (heliocentric
position, AU, and velocity, AU/day, ICRF equatorial), `covariance` (6×6,
from the weighted fit), `elements` (`a`, `e`, `i`, `node`, `peri`, `M`, `q`:
osculating, heliocentric ecliptic J2000, degrees), `sigma_elements` (their
1σ, by linear propagation of the covariance), `residuals` (arcsec, `ra`
times cos δ and `dec`, per observation), `used` (`false` for the outliers),
`rms` (arcsec, of the used ones), `arc_days`.
"""
function fit_orbit(ra::AbstractVector{<:Real}, dec::AbstractVector{<:Real}, epochs::AbstractVector{<:Real},
                   observatory; sigma=0.5, outlier_sigma::Real=3.0)
    n = length(ra)
    (length(dec) == n && length(epochs) == n) || throw(ArgumentError("ra, dec and epochs must have the same length"))
    n >= 3 || throw(ArgumentError("an orbit needs at least 3 observations"))
    order = sortperm(epochs)
    ra, dec, epochs = Float64.(ra[order]), Float64.(dec[order]), Float64.(epochs[order])
    codes = observatory isa AbstractString ? fill(String(observatory), n) : String.(collect(observatory)[order])
    σ = sigma isa Real ? fill(Float64(sigma), n) : Float64.(collect(sigma)[order])
    t = [jd + _tt_minus_utc(jd) for jd in epochs]
    R = Vector{Vector{Float64}}(undef, n)
    for code in unique(codes)
        k = findall(==(code), codes)
        R[k] = collect.(_observer_heliocentric_tt(t[k], code))
    end
    return _fit_orbit(ra, dec, t, R, σ, outlier_sigma)
end

# fit_orbit given TT epochs `t` and observer positions `R`, sorted in time.
# Faster starts, for many fits: ranging over `distances` only, or one state
# `from` = (state, epoch) already near.
function _fit_orbit(ra, dec, t, R, σ, outlier_sigma; distances=nothing, from=nothing)
    n = length(ra)
    mid = (1 + n) ÷ 2
    t0 = t[mid]
    starts = if from !== nothing
        r, v = _kepler_propagate(from[1][1:3], from[1][4:6], t0 - from[2])
        [vcat(r, v)]
    elseif distances !== nothing
        _ranging_orbits(ra, dec, t, R, σ, t0; distances=distances)
    else
        picks = (1, mid, n)
        ρ̂ = [_line_of_sight(ra[i], dec[i]) for i in picks]
        vcat([vcat(r, v) for (r, v) in _gauss_orbits(ρ̂, [R[i] for i in picks], [t[i] for i in picks])],
             _ranging_orbits(ra, dec, t, R, σ, t0))
    end
    isempty(starts) && error("no first orbit found: the arc may be too short")

    use = trues(n)
    best = nothing
    for x0 in starts
        x, c = _orbit_least_squares(x0, t0, ra, dec, R, t, σ, use)
        far = all(i -> norm(x[1:3] .- R[i]) > _MIN_DISTANCE, 1:n)
        far && isfinite(c) && (best === nothing || c < best[2]) && (best = (x, c))
    end
    best === nothing && error("no orbit converged away from the observer")
    x = best[1]
    while true
        res = _orbit_residuals(x, t0, ra, dec, R, t)
        norms = [hypot(res[2i-1], res[2i]) / σ[i] for i in 1:n]
        worst = argmax(i -> use[i] ? norms[i] : -Inf, 1:n)
        (use[worst] && norms[worst] > outlier_sigma && count(use) > 3) || break
        use[worst] = false
        x, _ = _orbit_least_squares(x, t0, ra, dec, R, t, σ, use)
    end

    res = _orbit_residuals(x, t0, ra, dec, R, t)
    w = repeat(1 ./ σ .^ 2; inner=2) .* repeat(Float64.(use); inner=2)
    J = _orbit_jacobian(x, t0, ra, dec, R, t)
    covariance = inv(Symmetric(J' * (w .* J)))
    elements = _state_to_elements(x[1:3], x[4:6])
    # d(elements)/d(state), numerically, for their 1σ
    names = keys(elements)
    G = zeros(length(names), 6)
    for j in 1:6
        h = j <= 3 ? 1e-8 : 1e-10
        xp, xm = copy(x), copy(x)
        xp[j] += h; xm[j] -= h
        ep, em = _state_to_elements(xp[1:3], xp[4:6]), _state_to_elements(xm[1:3], xm[4:6])
        for (k, name) in enumerate(names)
            d = ep[name] - em[name]
            name in (:node, :peri, :M) && (d = mod(d + 180, 360) - 180)
            G[k, j] = d / 2h
        end
    end
    σe = sqrt.(max.(diag(G * covariance * G'), 0.0))
    used_res = [hypot(res[2i-1], res[2i]) for i in 1:n if use[i]]
    return (epoch=t0, state=x, covariance=Matrix(covariance), elements=elements,
            sigma_elements=NamedTuple{names}(Tuple(σe)),
            residuals=Table(ra=res[1:2:end], dec=res[2:2:end]), used=use,
            rms=sqrt(sum(abs2, used_res) / (2length(used_res))), arc_days=t[end] - t[1])
end

"""
    orbit_ephemeris(orbit, epochs, observatory) -> Table

Where a [`fit_orbit`](@ref) orbit puts the object as seen from the MPC site
`observatory` at `epochs` (JD UTC): `ra`, `dec` (degrees), `delta` and `r`
(AU) and `sigma` (arcsec), the semi-major axis of its 1σ uncertainty ellipse
on the sky from the fit's covariance, linearly propagated. A live request to
JPL Horizons for the observer's positions.
"""
function orbit_ephemeris(orbit, epochs::AbstractVector{<:Real}, observatory::AbstractString)
    t = [jd + _tt_minus_utc(jd) for jd in epochs]
    R = collect.(_observer_heliocentric_tt(t, observatory))
    x, t0 = orbit.state, orbit.epoch
    out = Table(ra=Float64[], dec=Float64[], delta=Float64[], r=Float64[], sigma=Float64[])
    for (ti, Ri) in zip(t, R)
        α, δ = _orbit_radec(x[1:3], x[4:6], t0, Ri, ti)
        J = zeros(2, 6)
        for j in 1:6
            h = j <= 3 ? 1e-7 : 1e-9
            xp, xm = copy(x), copy(x)
            xp[j] += h; xm[j] -= h
            ap, dp = _orbit_radec(xp[1:3], xp[4:6], t0, Ri, ti)
            am, dm = _orbit_radec(xm[1:3], xm[4:6], t0, Ri, ti)
            J[1, j] = 3600 * (mod(ap - am + 180, 360) - 180) * cosd(δ) / 2h
            J[2, j] = 3600 * (dp - dm) / 2h
        end
        p, _ = _kepler_propagate(x[1:3], x[4:6], ti - t0)
        push!(out, (ra=α, dec=δ, delta=norm(p .- Ri), r=norm(p), sigma=sqrt(maximum(eigvals(Symmetric(J * orbit.covariance * J'))))))
    end
    return out
end
