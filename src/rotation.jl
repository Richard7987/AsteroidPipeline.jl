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

Each frame's WCS is refined and its zero point fitted against Gaia DR3 stars
(`reference_stars`, or fetched once for the first frame's field), as in
[`candidate_magnitudes`](@ref). Where the aperture's signal-to-noise reaches
`recentre_snr`, it is moved to the object's flux-weighted centroid; elsewhere
the flux is measured at the predicted position (forced photometry), which
keeps faint frames in the curve without biasing them upward. The background
is the median of `annulus`; `flux_err` combines its scatter, the error of
that median and, with `gain` (e-/ADU, else the header's `GAIN` if it has
one), the object's shot noise.

Returns a table, one row per frame: `frame`, `epoch`, `ra`, `dec` (where it
was measured), `x`, `y`, `recentred`, `star_nearby`, `flux`, `flux_err`, `mag`,
`mag_err` (Gaia G-equivalent; `NaN` where the frame has no zero point or the
flux isn't positive) and `zeropoint`. `star_nearby` marks a frame where a
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

    rows = map(enumerate(fits_paths)) do (k, path)
        image, header, hdr = FITS(path, "r") do f
            permutedims(Float64.(read(f[1]))), read_header(f[1], String), read_header(f[1])
        end
        epoch = frame_epoch(hdr; timestamp_key, exptime_key)
        mask = fill_value_mask(image)
        wcs = load_wcs(header)
        reference_stars === nothing && (reference_stars = _field_reference_stars(image, wcs, epoch))
        refined = refine_wcs(image, wcs, reference_stars; mask)
        refined.refined && (wcs = refined.wcs)
        zp = photometric_zeropoint(image, wcs, reference_stars; mask, aperture_radius).zeropoint
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
        good = isfinite(zp) && m.flux > 0
        mag = good ? zp - 2.5 * log10(m.flux) : NaN
        # a star counts if it could add a tenth of the object's flux
        near = 2 * aperture_radius * _pixel_scale_arcsec(wcs, size(image))
        star_nearby = any(s -> !(s.gmag > mag + 2.5) && _angular_distance_arcsec(ra, dec, s.ra, s.dec) <= near,
                          reference_stars)
        (frame=k, epoch=epoch, ra=ra, dec=dec, x=x, y=y, recentred=recentred, star_nearby=star_nearby,
         flux=m.flux, flux_err=m.flux_err,
         mag=mag,
         mag_err=good ? 2.5 / log(10) * m.flux_err / m.flux : NaN, zeropoint=zp)
    end
    return Table(rows)
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
        (period, power, false_alarm_probability, periodogram)

Recover a periodic signal (a rotating object's changing brightness) from a
light curve via the Lomb-Scargle periodogram (`LombScargle.jl`), searching
periods between `minimum_period` and `maximum_period` (same time unit as
`times`, e.g. days if `times` are Julian Dates as [`light_curve`](@ref)
returns).

Returns the period and power of the periodogram's highest peak, and its
false-alarm probability (`LombScargle.fap`) — the probability that noise
alone would produce a peak at least this strong; a small value is what
distinguishes a real periodic signal from a noise fluctuation. `periodogram`
is the full `LombScargle.Periodogram`, for inspecting other peaks (e.g. a
harmonic of the true period, a common ambiguity) with `LombScargle.jl`'s
own tools directly.
"""
function recover_rotation_period(times::AbstractVector{<:Real}, flux::AbstractVector{<:Real};
                                  minimum_period::Real, maximum_period::Real)
    minimum_period < maximum_period ||
        throw(ArgumentError("minimum_period must be less than maximum_period"))

    periodogram = lombscargle(times, flux;
                               minimum_frequency=1 / maximum_period,
                               maximum_frequency=1 / minimum_period)
    power = findmaxpower(periodogram)
    period = only(findmaxperiod(periodogram))
    false_alarm_probability = fap(periodogram, power)

    return (period=period, power=power, false_alarm_probability=false_alarm_probability,
            periodogram=periodogram)
end
