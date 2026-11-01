"""
    light_curve(fits_paths, ra::Real, dec::Real;
                timestamp_key::AbstractString="MJD-OBS",
                aperture_radius::Real=3.0) -> (times, flux, flux_err)

Build a light curve for a confirmed object at a fixed sky position by
forced aperture photometry: for each frame in `fits_paths`, its WCS
(via [`load_wcs`](@ref)) gives the pixel position of `(ra, dec)`, and flux
is measured there directly rather than at a `detect_sources`-found peak —
appropriate for rotation-period follow-up, where the object's position is
already known precisely (from [`astrometric_calibrate`](@ref) output) and
barely moves between exposures taken close together in time, unlike the
original discovery epochs. For a target that moves appreciably over the
follow-up sequence, `ra`/`dec` need to be recomputed per frame from an
ephemeris — not done here.

Background/noise come from `Photometry.Background.estimate_background`
as in [`detect_sources`](@ref); `flux_err` is `Photometry.photometry`'s
propagated aperture error from a uniform per-pixel `noise`, not a full
per-pixel variance map.

Returns `times` (Julian Date), `flux`, and `flux_err`, ready for
[`recover_rotation_period`](@ref).
"""
function light_curve(fits_paths::AbstractVector{<:AbstractString}, ra::Real, dec::Real;
                      timestamp_key::AbstractString="MJD-OBS", aperture_radius::Real=3.0)
    times = Float64[]
    flux = Float64[]
    flux_err = Float64[]

    for path in fits_paths
        FITS(path, "r") do f
            hdu = f[1]
            wcs = load_wcs(read_header(hdu, String))
            x, y = WCS.world_to_pix(wcs, Float64[ra, dec])

            # raw (x, y) order, as world_to_pix gives positions; not permuted
            image = Float64.(read(hdu))
            background, noise = estimate_background(image; location=SourceExtractorBackground(), rms=MADStdRMS())
            error_map = fill(noise, size(image))

            ap = CircularAperture(x, y, aperture_radius)
            result = photometry(ap, image .- background, error_map)
            push!(flux, result.aperture_sum)
            push!(flux_err, result.aperture_sum_err)

            mjd = read_key(hdu, timestamp_key)[1]
            push!(times, mjd + 2400000.5)
        end
    end

    return (times=times, flux=flux, flux_err=flux_err)
end

"""
    moving_light_curve(fits_paths, positions; aperture_radius=5.0, annulus=(10.0, 15.0),
                       recentre_snr=5.0, timestamp_key=nothing, exptime_key="EXPTIME",
                       gain=nothing, reference_stars=nothing) -> Table

The light curve of a moving object: its calibrated brightness in every frame
of `fits_paths`, measured where it is in each, whether or not it was detected
there.

`positions` are rows with `ra`, `dec` (degrees) and `epoch` (JD UTC) — the
object's tracklet from [`astrometric_calibrate`](@ref), or an ephemeris.
Linear motion in RA and Dec, fitted to them by least squares, gives its
position at each frame's mid-exposure ([`frame_epoch`](@ref)); one position
means a fixed target. Within a night this is better than an orbit, since it
comes from the object itself; it doesn't follow curvature over days.

Each frame's WCS is refined against Gaia DR3 stars (`reference_stars`, or
fetched once for the first frame's field); a frame without one, or whose
header WCS can't be refined, is solved from its pointing
([`solve_wcs`](@ref)); one where none is found is kept, with a warning, as
a row of `NaN`. Frames are measured in parallel on Julia's threads (start
Julia with `-t auto` to use every core); the result doesn't depend on how
many. The zero points come from the
same comparison stars in every frame: catalogue stars inside all of them,
measured with the object's aperture at signal-to-noise 30 or more, away
from its track, and dropping any whose brightness varies far more than the
others' (a saturated star). Each star's colour differs from the frames'
band, unfiltered or not, by its own constant amount, which a fixed set of
stars cancels from frame to frame; a set changing with each frame's
detections would not. With fewer than 3 such stars, each frame's zero point
is fitted on its own ([`photometric_zeropoint`](@ref)), as is a frame in which
fewer than 3 were measured (under cloud), the others keeping theirs. Where
the aperture's signal-to-noise reaches
`recentre_snr`, it is moved to the object's flux-weighted centroid; elsewhere
the flux is measured at the predicted position (forced photometry), which
keeps faint frames in the curve without biasing them upward. The background
is the median of `annulus`; `flux_err` combines its scatter, the error of
that median and, with `gain` (e-/ADU, else the header's `GAIN` if it has
one), the object's shot noise; `mag_err` adds the zero point's error,
`zeropoint_err` (the comparison stars' spread about the frame's offset,
over the square root of their number).

Returns a table, one row per frame: `frame`, `epoch`, `ra`, `dec` (where it
was measured), `x`, `y`, `recentred`, `star_nearby`, `flux`, `flux_err`, `mag`,
`mag_err` (Gaia G-equivalent; `NaN` where the frame has no zero point or the
flux isn't positive), `zeropoint` and `zeropoint_err`. `star_nearby` marks a frame where a
reference star bright enough to add a tenth of the object's flux (at most
2.5 mag fainter) lies within two aperture radii of it: crossing a star
makes it look brighter, so leave those frames out before passing
`epoch` and `mag` (or `flux`) to [`recover_rotation_period`](@ref).
"""
function moving_light_curve(fits_paths::AbstractVector{<:AbstractString}, positions;
                            aperture_radius::Real=5.0, annulus=(10.0, 15.0), recentre_snr::Real=5.0,
                            timestamp_key::Union{Nothing,AbstractString}=nothing,
                            exptime_key::Union{Nothing,AbstractString}="EXPTIME",
                            gain::Union{Nothing,Real}=nothing, reference_stars=nothing)
    pts = collect(positions)
    isempty(pts) && throw(ArgumentError("positions is empty"))
    t = [Float64(p.epoch) for p in pts]
    ra0 = Float64(pts[1].ra)
    ra_at, dec_at = if length(pts) == 1
        (_ -> ra0, _ -> Float64(pts[1].dec))
    else
        a_ra, b_ra = _linfit(t, [ra0 + mod(p.ra - ra0 + 180, 360) - 180 for p in pts])
        a_dec, b_dec = _linfit(t, [Float64(p.dec) for p in pts])
        (e -> mod(a_ra + b_ra * e, 360), e -> a_dec + b_dec * e)
    end

    stars = reference_stars
    solve_cache = _SolveCache()
    fits_lock = ReentrantLock()
    function measure(k, path)
        image, header, hdr = lock(fits_lock) do
            FITS(path, "r") do f
                permutedims(Float64.(read(f[1]))), read_header(f[1], String), read_header(f[1])
            end
        end
        epoch = frame_epoch(hdr; timestamp_key, exptime_key)
        mask = fill_value_mask(image)
        wcs = try
            first(_header_or_solved_wcs(image, header, hdr, epoch, solve_cache; reference_stars))
        catch e
            e isa ErrorException || rethrow(e)
            @warn "frame left out of the light curve: $(sprint(showerror, e))" path
            return (frame=k, epoch=epoch, ra=NaN, dec=NaN, x=NaN, y=NaN, recentred=false, flux=NaN,
                    flux_err=NaN, own=(zeropoint=NaN, n_stars=0, scatter=NaN), v=Float64[], scale=NaN)
        end
        stars === nothing && (stars = _field_reference_stars(image, wcs, epoch))
        refined = refine_wcs(image, wcs, stars; mask)
        pointing = frame_pointing(hdr)
        if !refined.refined && pointing !== nothing
            # a header WCS can be wrong beyond refining: a mirrored axis
            solved = solve_wcs(image, pointing..., stars; mask)
            solved.wcs === nothing || (refined = solved)
        end
        refined.refined && (wcs = refined.wcs)
        own = photometric_zeropoint(image, wcs, stars; mask, aperture_radius)
        g = gain !== nothing ? gain : haskey(hdr, "GAIN") ? Float64(hdr["GAIN"]) : nothing

        x, y = world_to_pix(wcs, Float64[ra_at(epoch), dec_at(epoch)])
        m = _moving_aperture(image, x, y, aperture_radius, annulus, mask, g)
        recentred = false
        if m.flux / m.flux_err >= recentre_snr
            d1, d2 = _refine_centroid(image .- m.background, round(Int, y), round(Int, x), aperture_radius)
            if hypot(d2 - x, d1 - y) <= aperture_radius
                x, y = d2, d1
                m = _moving_aperture(image, x, y, aperture_radius, annulus, mask, g)
                recentred = true
            end
        end
        ra, dec = pix_to_world(wcs, Float64[x, y])
        # each star's G + 2.5 log10(flux): NaN off the frame, near the object or too faint
        ny, nx = size(image)
        edge = annulus[2] + 2
        v = map(stars) do s
            sx, sy = world_to_pix(wcs, Float64[s.ra, s.dec])
            edge <= sx <= nx - edge && edge <= sy <= ny - edge || return NaN
            hypot(sx - x, sy - y) > 3 * aperture_radius || return NaN
            f = _moving_aperture(image, sx, sy, aperture_radius, annulus, mask, g)
            f.flux > 0 && f.flux / f.flux_err >= _COMPARISON_MIN_SNR ? s.gmag + 2.5 * log10(f.flux) : NaN
        end
        (frame=k, epoch=epoch, ra=ra, dec=dec, x=x, y=y, recentred=recentred, flux=m.flux,
         flux_err=m.flux_err, own=own, v=v, scale=_pixel_scale_arcsec(wcs, size(image)))
    end
    # in order until a frame is measured (the field's stars and Gaia queries
    # fetched once), then the rest on every thread
    measured = Vector{Any}(undef, length(fits_paths))
    k = 0
    while k < length(fits_paths) && (k == 0 || isempty(measured[k].v))
        k += 1
        measured[k] = measure(k, fits_paths[k])
    end
    Threads.@threads for j in k+1:length(fits_paths)
        measured[j] = measure(j, fits_paths[j])
    end
    frames = [f for f in measured]

    own = [f.own.zeropoint for f in frames]
    own_err = [f.own.n_stars >= 2 ? 1.2533 * f.own.scatter / sqrt(f.own.n_stars) : NaN for f in frames]
    zp, zp_err = _comparison_zeropoints([f.v for f in frames])
    if zp === nothing
        zp, zp_err = own, own_err
    else
        # the level from every frame's own fit, many stars' colours averaged;
        # a few comparison stars' would differ from one set of frames to another
        d = filter(isfinite, own .- zp)
        isempty(d) || (zp = zp .+ median(d))
        # a frame too few comparison stars were measured in (cloud) keeps its own
        lost = .!isfinite.(zp)
        zp[lost] = own[lost]
        zp_err[lost] = own_err[lost]
    end
    rows = map(enumerate(frames)) do (k, f)
        good = isfinite(zp[k]) && f.flux > 0
        mag = good ? zp[k] - 2.5 * log10(f.flux) : NaN
        near = 2 * aperture_radius * f.scale
        star_nearby = any(s -> !(s.gmag > mag + 2.5) && _angular_distance_arcsec(f.ra, f.dec, s.ra, s.dec) <= near,
                          stars)
        phot_err = 2.5 / log(10) * f.flux_err / f.flux
        (frame=f.frame, epoch=f.epoch, ra=f.ra, dec=f.dec, x=f.x, y=f.y, recentred=f.recentred,
         star_nearby=star_nearby, flux=f.flux, flux_err=f.flux_err, mag=mag,
         mag_err=good ? hypot(phot_err, isfinite(zp_err[k]) ? zp_err[k] : 0.0) : NaN,
         zeropoint=zp[k], zeropoint_err=zp_err[k])
    end
    return Table(rows)
end

const _COMPARISON_MIN_SNR = 30.0

# Zero points from the stars measured in every frame (`v[k][i]` = G + 2.5 log10
# flux): each frame's median offset from the stars' means, and its error. Stars
# scattering more than 3 times the typical star are dropped first. `nothing`
# with fewer than 3 stars. A frame with fewer than 3 stars measured (none at
# all, or under cloud) is left out of the choice of stars, and gets NaN.
function _comparison_zeropoints(v)
    measured = findall(vk -> count(isfinite, vk) >= 3, v)
    isempty(measured) && return nothing, nothing
    n = length(v[first(measured)])
    common = [i for i in 1:n if all(isfinite(v[k][i]) for k in measured)]
    length(common) >= 3 || return nothing, nothing
    V = [v[k][i] for k in measured, i in common]
    offsets() = [median(V[k, :] .- vec(mean(V; dims=1))) for k in axes(V, 1)]
    spread = vec(std(V .- offsets(); dims=1))
    keep = spread .<= 3 * median(spread)
    count(keep) >= 3 && (V = V[:, keep])
    level = median(vec(mean(V; dims=1)))
    off = offsets()
    zp = fill(NaN, length(v))
    zp[measured] = level .+ off
    zp_err = fill(NaN, length(v))
    zp_err[measured] = map(axes(V, 1)) do k
        d = V[k, :] .- vec(mean(V; dims=1)) .- off[k]
        1.2533 * 1.4826 * median(abs.(d)) / sqrt(size(V, 2))
    end
    return zp, zp_err
end

# Aperture flux at (x, y) with its error: the annulus's robust scatter over the
# aperture's pixels, the error of the annulus median, and shot noise with `gain`.
function _moving_aperture(image, x, y, aperture_radius, annulus, mask, gain)
    px = _aperture_pixels(image, x, y, aperture_radius, annulus, mask)
    px === nothing && return (flux=NaN, flux_err=NaN, background=NaN)
    background = median(px.ring)
    n = length(px.inside)
    σ = 1.4826 * median(abs.(px.ring .- background))
    flux = sum(px.inside) - background * n
    variance = n * σ^2 + (π / 2) * n^2 * σ^2 / length(px.ring) +
               (gain === nothing ? 0.0 : max(flux, 0.0) / gain)
    return (flux=flux, flux_err=sqrt(variance), background=background)
end

"""
    recover_rotation_period(times, flux; minimum_period, maximum_period) ->
        (period, peak_period, power, false_alarm_probability, periodogram)

Recover a periodic signal (a rotating object's changing brightness) from a
light curve via the Lomb-Scargle periodogram (`LombScargle.jl`), searching
periods between `minimum_period` and `maximum_period` (same time unit as
`times`, e.g. days if `times` are Julian Dates as [`light_curve`](@ref)
returns).

`peak_period` and `power` are the periodogram's highest peak, and
`false_alarm_probability` (`LombScargle.fap`) the probability that noise
alone would produce a peak at least this strong; a small value is what
distinguishes a real periodic signal from a noise fluctuation.

A light curve with two maxima per cycle, as an elongated asteroid's, puts
the highest peak at half the period. `period` is twice `peak_period` where
alternate cycles differ: four harmonics of the double period fit
significantly better (F-test, 0.1%) than two of the peak's. It needs the
double period within `maximum_period` and the light curve to cover it;
otherwise, or where the two halves can't be told apart, `period` is
`peak_period`. [`fourier_period`](@ref) fits the period itself and gives
its uncertainty. `periodogram` is the full `LombScargle.Periodogram`, for
inspecting other peaks with `LombScargle.jl`'s own tools directly.
"""
function recover_rotation_period(times::AbstractVector{<:Real}, flux::AbstractVector{<:Real};
                                  minimum_period::Real, maximum_period::Real)
    minimum_period < maximum_period ||
        throw(ArgumentError("minimum_period must be less than maximum_period"))

    periodogram = lombscargle(times, flux;
                               minimum_frequency=1 / maximum_period,
                               maximum_frequency=1 / minimum_period)
    power = findmaxpower(periodogram)
    peak_period = only(findmaxperiod(periodogram))
    false_alarm_probability = fap(periodogram, power)
    period = _double_if_halves_differ(times, flux, peak_period, maximum_period)

    return (period=period, peak_period=peak_period, power=power,
            false_alarm_probability=false_alarm_probability, periodogram=periodogram)
end

# Two harmonics of P against four of 2P, which add only the difference between
# alternate cycles: twice the period where an F-test finds it at 0.1%.
function _double_if_halves_differ(times, flux, P, maximum_period)
    span = maximum(times) - minimum(times)
    (2P <= maximum_period && 2P <= 0.8span && length(times) >= 20) || return P
    window(Q) = (minimum_period=1 / (1 / Q + 1 / span), maximum_period=min(1 / (1 / Q - 1 / span), 0.8span))
    one = fourier_period(times, flux; window(P)..., order=2, trend=true)
    two = fourier_period(times, flux; window(2P)..., order=4, trend=true)
    d2 = length(times) - 10
    rss1, rss2 = one.rms^2 * (length(times) - 6), two.rms^2 * d2
    F = max((rss1 - rss2) / 4, 0.0) / (rss2 / d2)
    x, a = d2 / (d2 + 4F), d2 / 2
    return x^a * (1 + a * (1 - x)) < 0.001 ? two.period : P
end

"""
    fourier_period(times, mag; minimum_period, maximum_period, order=4, trend=false,
                   mag_err=nothing) -> (period, period_error, amplitude, rms, periods, chi2)

An asteroid's rotation period from its light curve, the way asteroid
photometry measures it: a Fourier series of the rotation phase to `order`
harmonics, fitted to the magnitudes by least squares at each trial period
between `minimum_period` and `maximum_period` (in the unit of `times`); the
period is where the fit's χ² is least.

An elongated asteroid brightens twice per rotation, and its two maxima
and two minima differ. One sinusoid, as in [`recover_rotation_period`](@ref)'s
periodogram, fits such a curve best at half the period; the second and
higher harmonics tell the two halves apart. 4 is usual for a night or a
few; a higher `order` fits more detail and needs more points. With
`trend`, a straight line in time is fitted with the series, for a drift
through the night that isn't the rotation (an unfiltered colour as the
airmass changes). `mag_err` weights the points.

`period_error` is where χ² grows by 1 from its least, the residuals scaled
to unit variance; `amplitude` is the fitted curve's peak to peak, and `rms`
the residuals'. `periods` and `chi2` (reduced) are the whole search, to see
its other minima: half or twice the period, or aliases of the gaps between
nights. When ``2P`` fits about as well as ``P``, compare the light curve
folded on each. A period longer than about 0.8 of the time the light curve
covers is warned of: a Fourier series that long fits any drift.
"""
function fourier_period(times::AbstractVector{<:Real}, mag::AbstractVector{<:Real};
                        minimum_period::Real, maximum_period::Real, order::Integer=4,
                        trend::Bool=false, mag_err::Union{Nothing,AbstractVector{<:Real}}=nothing)
    0 < minimum_period < maximum_period ||
        throw(ArgumentError("need 0 < minimum_period < maximum_period"))
    order >= 1 || throw(ArgumentError("order must be at least 1"))
    length(times) == length(mag) || throw(ArgumentError("times and mag must have the same length"))
    t = Float64.(times) .- mean(times)
    m = Float64.(mag)
    w = mag_err === nothing ? ones(length(m)) : 1 ./ Float64.(mag_err)
    nparam = 1 + trend + 2order
    dof = length(m) - nparam
    dof > 0 || throw(ArgumentError("$(length(m)) points can't fit $nparam parameters"))

    design(P) = hcat(ones(length(t)), (trend ? (t,) : ())...,
                     (f.(2π * k .* t ./ P) for k in 1:order for f in (sin, cos))...)
    function fit(P)
        X = design(P)
        β = (w .* X) \ (w .* m)
        return sum(abs2, w .* (m .- X * β)), β, X
    end
    # a frequency step that moves the highest harmonic by a tenth of a cycle
    # over the span; then a step 100 times finer around the least
    span = maximum(t) - minimum(t)
    step = 1 / (10 * order * span)
    n = max(2, ceil(Int, (1 / minimum_period - 1 / maximum_period) / step) + 1)
    freqs = collect(range(1 / maximum_period, 1 / minimum_period; length=n))
    chi = [first(fit(1 / f)) for f in freqs]
    i = argmin(chi)
    lo, hi = freqs[max(i - 1, 1)], freqs[min(i + 1, end)]
    local fine, fine_chi, j, inside
    # narrowed until the Δχ² = 1 interval spans several steps, so its width is measured
    for _ in 1:6
        fine = collect(range(lo, hi; length=201))
        fine_chi = [first(fit(1 / f)) for f in fine]
        j = argmin(fine_chi)
        inside = fine_chi ./ fine_chi[j] .* dof .<= dof + 1
        count(inside) >= 10 && break
        lo, hi = fine[max(j - 2, 1)], fine[min(j + 2, end)]
    end
    period = 1 / fine[j]
    within = 1 ./ fine[inside]
    period_error = (maximum(within) - minimum(within)) / 2

    # a period about as long as the series fits a drift as well as a rotation
    period > 0.8 * span &&
        @warn "the period found, $(round(period; sigdigits=4)), is $(round(Int, 100 * period / span))% of the " *
              "time the light curve covers: a drift through it fits as well as a rotation; observe longer, " *
              "fit `trend`, or lower `maximum_period`"
    best, β, X = fit(period)
    first_harmonic = 2 + trend
    phase = range(0, 1; length=1001)
    curve = [sum(β[first_harmonic+2(k-1)] * sinpi(2k * φ) + β[first_harmonic+2(k-1)+1] * cospi(2k * φ) for k in 1:order)
             for φ in phase]
    return (period=period, period_error=period_error, amplitude=maximum(curve) - minimum(curve),
            rms=sqrt(sum(abs2, m .- X * β) / dof), periods=1 ./ freqs, chi2=chi ./ dof)
end

"""
    object_geometry(object, epochs; observatory=nothing) -> Table

A solar-system object's geometry at each UTC Julian Date of `epochs`, from
JPL Horizons (a live request): `r` and `delta`, its distances from the Sun
and from the observer (AU), and `phase`, the Sun–object–observer angle
(degrees). `object` is looked up among small bodies: a number (`624`), a
designation (`"2019 UA185"`) or a name. The observer is the MPC site
`observatory`, or the Earth's centre.
"""
function object_geometry(object, epochs::AbstractVector{<:Real}; observatory=nothing)
    isempty(epochs) && return Table(r=Float64[], delta=Float64[], phase=Float64[])
    center = observatory === nothing ? "500@399" : "$(observatory)@399"
    jds = sort(unique(Float64.(epochs)))
    at = Dict{Float64,NamedTuple{(:r, :delta, :phase),NTuple{3,Float64}}}()
    query = Dict("COMMAND" => "'$(object);'", "EPHEM_TYPE" => "OBSERVER", "CENTER" => "'$center'",
                 "TIME_TYPE" => "UT", "QUANTITIES" => "'19,20,24'", "CSV_FORMAT" => "YES", "OBJ_DATA" => "NO")
    results = _horizons(query, [@sprintf("%.8f", jd) for jd in jds])
    for (chunk, result) in zip(Iterators.partition(jds, 50), results)
        rows = _horizons_geometry(result)
        length(rows) == length(chunk) || error("Horizons returned $(length(rows)) rows for $(length(chunk)) epochs")
        # Horizons answers in time order
        merge!(at, Dict(zip(chunk, rows)))
    end
    return Table([at[Float64(jd)] for jd in epochs])
end

function _horizons_geometry(result::AbstractString)
    m = match(r"([^\n]*)\n\*+\n\$\$SOE\n(.*)\$\$EOE"s, result)
    m === nothing && error("unexpected Horizons response: $(first(result, 300))")
    names = strip.(split(m[1], ','))
    col(name) = something(findfirst(==(name), names), 0)
    ir, id, ip = col("r"), col("delta"), col("S-T-O")
    min(ir, id, ip) > 0 || error("Horizons response lacks r, delta or S-T-O: $(m[1])")
    return map(split(strip(m[2]), '\n')) do line
        f = strip.(split(line, ','))
        (r=parse(Float64, f[ir]), delta=parse(Float64, f[id]), phase=parse(Float64, f[ip]))
    end
end

"""
    reduced_magnitudes(epochs, mags, geometry; G=0.15) -> (epoch, mag)

A moving object's light curve made comparable from night to night. Each
magnitude becomes the object's at 1 AU from the Sun and from the observer,
`mag - 5 log10(r * delta)`, and each epoch moves back by the light time, to
when the light left the object. `geometry` holds `r`, `delta` (AU) and
`phase` (degrees) for each epoch, as [`object_geometry`](@ref) returns.

With `G`, the brightening toward opposition is removed as well, by the H-G
phase function, which leaves the absolute magnitude H plus the rotation's
variation; 0.15 is the usual slope when the object's own isn't known.
`G=nothing` leaves it in, to be fitted together with the period.

Pass the result to [`recover_rotation_period`](@ref).
"""
function reduced_magnitudes(epochs::AbstractVector{<:Real}, mags::AbstractVector{<:Real}, geometry;
                            G::Union{Nothing,Real}=0.15)
    g = collect(geometry)
    length(epochs) == length(mags) == length(g) ||
        throw(ArgumentError("epochs, mags and geometry must have the same length"))
    epoch = [t - p.delta / _AU_PER_DAY_LIGHT for (t, p) in zip(epochs, g)]
    mag = [m - 5log10(p.r * p.delta) - (G === nothing ? 0.0 : _hg_phase_dimming(deg2rad(p.phase), G))
           for (m, p) in zip(mags, g)]
    return (epoch=epoch, mag=mag)
end
