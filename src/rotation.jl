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

            # raw FITS (x, y) order matches world_to_pix's convention directly
            # (unlike detect_sources, which needs the permuted image) — see
            # the axis-order note in src/pipeline.jl.
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
