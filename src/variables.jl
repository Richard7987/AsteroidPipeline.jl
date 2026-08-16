"""
    variability_chi2(flux, flux_err) -> (chi2, dof)

Chi-squared goodness of fit of `flux` (with per-point uncertainty
`flux_err`) against the constant-flux (non-variable) null hypothesis,
using the error-weighted mean as the constant-flux estimate. `dof` is
`length(flux) - 1` (one free parameter, the mean).

A large `chi2 / dof` is evidence against the null hypothesis — i.e.
evidence of genuine flux variability rather than measurement noise.
Used by [`find_variable_sources`](@ref); exposed separately since a
caller may want the raw statistic rather than only a threshold decision.
"""
function variability_chi2(flux::AbstractVector{<:Real}, flux_err::AbstractVector{<:Real})
    length(flux) == length(flux_err) ||
        throw(ArgumentError("flux and flux_err must have the same length"))
    weights = 1.0 ./ flux_err .^ 2
    weighted_mean = sum(flux .* weights) / sum(weights)
    chi2 = sum(((flux .- weighted_mean) .^ 2) .* weights)
    return chi2, length(flux) - 1
end

"""
    find_variable_sources(detections_per_frame, timestamps;
                           position_tolerance::Real=2.0,
                           min_frames::Integer=length(detections_per_frame),
                           chi2_threshold::Real=3.0)

Match source detections across frames by consistent *position* (zero
assumed motion) rather than [`link_candidates`](@ref)'s linear-motion
model, then keep only sources whose flux is statistically inconsistent
with being constant — candidate variable stars or transients.

`detections_per_frame` is a vector of per-frame detection tables (as
returned by [`detect_sources`](@ref), each with `x`, `y`, `flux`,
`flux_err` columns); `timestamps` gives each frame's observation time,
matching [`link_candidates`](@ref)'s convention (used here only to keep
frame order meaningful for periodicity follow-up, not for any position
prediction).

For each detection in frame 1, the closest detection within
`position_tolerance` pixels of that *same* position in every other frame
is matched (via [`link_candidates`](@ref)'s `_closest_detection` — no
velocity model, since a variable star does not move between frames).
Groups reaching at least `min_frames` matched points are then tested
against [`variability_chi2`](@ref)'s constant-flux null hypothesis; a
group is kept only if its **reduced** chi-squared (`chi2 / dof`) exceeds
`chi2_threshold`.

What "keep only if variable" filters depends on how `detections_per_frame`
was produced:

- **ZOGY path** (`detections_per_frame` from `S_corr` — `run_pipeline` or
  `search_field` given a `reference`): a constant star differences to
  zero and is never detected in `S_corr` at all, so a matched stationary
  group is already a variable/transient by construction; the chi2 test
  here is a second filter against the known subtraction-residual failure
  mode at bright stars (see `INVESTIGATION_LOG.md`), not the primary
  discriminant.
- **Raw path** (no `reference`): every star above `threshold` is detected
  in every frame, so the chi2 test is the only thing separating a genuine
  variable from an ordinary constant star.

`chi2_threshold`'s right value is data-dependent (like `threshold` and
`quality_max_std` elsewhere in this package) — it assumes `flux_err` is
the propagated-aperture-error estimate `detect_sources` provides (from a
uniform per-pixel noise, not a full variance map), so treat the absolute
chi2 scale as approximate and tune per survey.

Returns a vector of matched-point groups, each a vector of
`(frame, x, y, flux, flux_err)` named tuples sorted by frame — the same
shape [`link_candidates`](@ref) returns (with two extra fields), so a
group can be passed directly to [`astrometric_calibrate`](@ref).

For periodicity, build `(times, flux)` from a returned group's own
`frame`/`flux` fields and `timestamps`, and pass to
[`recover_rotation_period`](@ref) directly.
"""
function find_variable_sources(detections_per_frame, timestamps;
                                position_tolerance::Real=2.0,
                                min_frames::Integer=length(detections_per_frame),
                                chi2_threshold::Real=3.0)
    nframes = length(detections_per_frame)
    nframes == length(timestamps) ||
        throw(ArgumentError("detections_per_frame and timestamps must have the same length"))

    Point = NamedTuple{(:frame, :x, :y, :flux, :flux_err),Tuple{Int,Float64,Float64,Float64,Float64}}
    groups = Vector{Point}[]
    nframes < 1 && return groups

    for d1 in detections_per_frame[1]
        group = Point[(frame=1, x=Float64(d1.x), y=Float64(d1.y),
                        flux=Float64(d1.flux), flux_err=Float64(d1.flux_err))]

        for k in 2:nframes
            best = _closest_detection(detections_per_frame[k], d1.x, d1.y, position_tolerance)
            best === nothing && continue
            push!(group, (frame=k, x=Float64(best.x), y=Float64(best.y),
                           flux=Float64(best.flux), flux_err=Float64(best.flux_err)))
        end

        length(group) < min_frames && continue

        chi2, dof = variability_chi2([p.flux for p in group], [p.flux_err for p in group])
        dof > 0 && chi2 / dof > chi2_threshold && push!(groups, group)
    end

    return groups
end
