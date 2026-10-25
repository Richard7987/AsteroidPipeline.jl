"""
    variability_chi2(flux, flux_err; systematic_error_fraction::Real=0.03) -> (chi2, dof)

Chi-squared goodness of fit of `flux` against the constant-flux
(non-variable) null hypothesis, using the error-weighted mean as the
constant-flux estimate. `dof` is `length(flux) - 1` (one free parameter,
the mean).

The per-point error is `sqrt(flux_err^2 + (systematic_error_fraction * flux)^2)`:
a systematic floor, proportional to flux, added in quadrature. Without it a
bright star's tiny formal error makes the test flag flat-fielding and PSF
variations that `flux_err` doesn't model; on ZTF data the false positives
were the *brightest* stars. The 3% default keeps three confirmed variables
in three ZTF fields above `chi2_threshold=10` while flagging 0.6-2% of
constant stars; a 5% floor reaches 0% only by losing two of the variables.

A large `chi2 / dof` is evidence against the null hypothesis — i.e.
evidence of genuine flux variability beyond both statistical noise and
this systematic floor. Used by [`find_variable_sources`](@ref); exposed
separately since a caller may want the raw statistic rather than only a
threshold decision.
"""
function variability_chi2(flux::AbstractVector{<:Real}, flux_err::AbstractVector{<:Real};
                           systematic_error_fraction::Real=0.03)
    length(flux) == length(flux_err) ||
        throw(ArgumentError("flux and flux_err must have the same length"))
    eff_err = sqrt.(flux_err .^ 2 .+ (systematic_error_fraction .* flux) .^ 2)
    weights = 1.0 ./ eff_err .^ 2
    weighted_mean = sum(flux .* weights) / sum(weights)
    chi2 = sum(((flux .- weighted_mean) .^ 2) .* weights)
    return chi2, length(flux) - 1
end

"""
    _filter_high_snr(detections, max_relative_error) -> Vector

`detections` (a table with `x`, `y`, `flux`, `flux_err` columns) restricted
to rows with positive flux and `flux_err / flux < max_relative_error` — the
S/N floor shared by [`photometric_scale`](@ref) and
[`find_variable_sources`](@ref), factored out so both apply exactly the
same cut.
"""
_filter_high_snr(detections, max_relative_error::Real) =
    [row for row in detections if row.flux > 0 && row.flux_err / row.flux < max_relative_error]

"""
    photometric_scale(detections_per_frame; position_tolerance::Real=2.0,
                       max_relative_error::Real=0.10, min_stars::Integer=5) -> Vector{Float64}

Per-frame photometric scale factor, relative to frame 1, by ensemble
differential photometry, independent of any zeropoint header keyword.

For each frame `k`, stars are matched to frame 1 by position (within
`position_tolerance` pixels) and the factor is the **median** of
`flux_1 / flux_k` over matches with `flux_err / flux < max_relative_error`
in both frames, so one noisy or blended star can't swing it. It corrects
any uniform per-frame offset, such as a fixed aperture enclosing a
seeing-dependent fraction of the PSF (9-13% between ZTF frames).

If fewer than `min_stars` matches survive for some frame, its factor is
left at `1.0` with a warning.

Returns a vector the same length as `detections_per_frame`; index 1 is
always `1.0`.
"""
function photometric_scale(detections_per_frame; position_tolerance::Real=2.0,
                            max_relative_error::Real=0.10, min_stars::Integer=5)
    nframes = length(detections_per_frame)
    scales = ones(Float64, nframes)
    nframes < 2 && return scales

    ref = _filter_high_snr(detections_per_frame[1], max_relative_error)
    for k in 2:nframes
        candidates = _filter_high_snr(detections_per_frame[k], max_relative_error)
        ratios = Float64[]
        for d1 in ref
            best = _closest_detection(candidates, d1.x, d1.y, position_tolerance)
            best === nothing && continue
            push!(ratios, d1.flux / best.flux)
        end

        if length(ratios) < min_stars
            @warn "too few high-S/N matched stars to measure frame's photometric scale; leaving uncorrected" frame=k n_matches=length(ratios) min_stars max_relative_error
            continue
        end
        scales[k] = median(ratios)
    end

    return scales
end

"""
    find_variable_sources(detections_per_frame, timestamps;
                           position_tolerance::Real=2.0,
                           min_frames::Integer=length(detections_per_frame),
                           chi2_threshold::Real=10.0,
                           normalize::Bool=true,
                           max_relative_error::Real=0.10,
                           systematic_error_fraction::Real=0.03)

Match source detections across frames by consistent *position* (zero
assumed motion) rather than [`link_candidates`](@ref)'s linear-motion
model, then keep only sources whose flux is statistically inconsistent
with being constant — candidate variable stars or transients.

`detections_per_frame` holds per-frame detection tables (from
[`detect_sources`](@ref), with `x`, `y`, `flux`, `flux_err` columns);
`timestamps` gives each frame's time, as for [`link_candidates`](@ref).

Before the test:

- **S/N floor** (`max_relative_error`, 10%): detections with
  `flux_err / flux` at or above it are dropped in every frame.
- **Photometric normalization** (`normalize`): each frame's flux is
  rescaled by [`photometric_scale`](@ref). Set `normalize=false` if `flux`
  is already calibrated, or on the ZOGY path.
- **Systematic error floor** (`systematic_error_fraction`, 3%), passed to
  [`variability_chi2`](@ref), where its calibration is described.

Each frame-1 detection is matched to the closest detection within
`position_tolerance` pixels in every other frame. Groups of at least
`min_frames` points whose **reduced** chi-squared (`chi2 / dof`) exceeds
`chi2_threshold` are kept. A small fraction of constant stars still pass:
confirm a candidate independently (a VSX/SIMBAD match via
[`crossmatch_catalog`](@ref), or a period from
[`recover_rotation_period`](@ref)). With the default `min_frames`, a field
where few sources clear the S/N floor in every frame returns few or none.

On the ZOGY path (`S_corr` detections) a constant star isn't detected at
all, so the test only filters subtraction residuals; on the raw path it is
what separates a variable from a constant star.

Returns groups of `(frame, x, y, flux, flux_err)` named tuples sorted by
frame, the shape [`link_candidates`](@ref) returns (with `flux` normalized
when `normalize=true`), ready for [`astrometric_calibrate`](@ref) or, with
`timestamps`, for [`recover_rotation_period`](@ref).
"""
function find_variable_sources(detections_per_frame, timestamps;
                                position_tolerance::Real=2.0,
                                min_frames::Integer=length(detections_per_frame),
                                chi2_threshold::Real=10.0,
                                normalize::Bool=true,
                                max_relative_error::Real=0.10,
                                systematic_error_fraction::Real=0.03)
    nframes = length(detections_per_frame)
    nframes == length(timestamps) ||
        throw(ArgumentError("detections_per_frame and timestamps must have the same length"))

    Point = NamedTuple{(:frame, :x, :y, :flux, :flux_err),Tuple{Int,Float64,Float64,Float64,Float64}}
    groups = Vector{Point}[]
    nframes < 1 && return groups

    scales = normalize ?
             photometric_scale(detections_per_frame; position_tolerance, max_relative_error) :
             ones(Float64, nframes)
    filtered = [_filter_high_snr(d, max_relative_error) for d in detections_per_frame]

    for d1 in filtered[1]
        group = Point[(frame=1, x=Float64(d1.x), y=Float64(d1.y),
                        flux=Float64(d1.flux) * scales[1], flux_err=Float64(d1.flux_err) * scales[1])]

        for k in 2:nframes
            best = _closest_detection(filtered[k], d1.x, d1.y, position_tolerance)
            best === nothing && continue
            push!(group, (frame=k, x=Float64(best.x), y=Float64(best.y),
                           flux=Float64(best.flux) * scales[k], flux_err=Float64(best.flux_err) * scales[k]))
        end

        length(group) < min_frames && continue

        chi2, dof = variability_chi2([p.flux for p in group], [p.flux_err for p in group];
                                      systematic_error_fraction)
        dof > 0 && chi2 / dof > chi2_threshold && push!(groups, group)
    end

    return groups
end
