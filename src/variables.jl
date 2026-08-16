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
differential photometry — the survey-agnostic way to put every frame's
flux on a common scale without depending on any zeropoint header keyword
(real IASC campaign headers are unverified — see the docs site's "Using
real IASC campaign data" section).

For each frame `k`, matches stars against frame 1 by position (within
`position_tolerance` pixels, via [`link_candidates`](@ref)'s
`_closest_detection`), and takes the **median** of `flux_1 / flux_k` over
those matches — restricted, in both frames, to stars with
`flux_err / flux < max_relative_error`, so one noisy or blended star can't
swing the whole frame's scale.

Measured directly on real ZTF data (field 451, 2019-10-23, 119 matched
stars): the ensemble ratio disagreed with the frames' own `MAGZP`
zeropoints by 9-13%, essentially *unchanged* by the S/N cut above
(-8.9% to -13.2% either way) — so this is not a low-S/N selection effect.
It tracks each frame's own `SEEING` instead (1.805-2.009 px across the
five frames): a fixed-radius aperture (`aperture_radius` in
[`detect_sources`](@ref)) encloses a seeing-dependent fraction of a star's
total PSF flux, the standard "aperture correction" effect in photometry.
`photometric_scale`'s ensemble approach corrects for exactly this kind of
uniform per-frame multiplicative offset, whatever its cause — it does not
need to know *why* frame `k`'s stars all read low or high, only that they
do, consistently.

If fewer than `min_stars` matches survive the S/N cut for some frame,
that frame's factor is left at `1.0` (no correction) with a `@warn`,
rather than trusting a scale derived from a handful of stars.

Returns a vector the same length as `detections_per_frame`; index 1 is
always `1.0` by construction (frame 1 is its own reference).
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
                           chi2_threshold::Real=50.0,
                           normalize::Bool=true,
                           max_relative_error::Real=0.10)

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

Two corrections are applied before any variability test:

- **Low-S/N floor** (`max_relative_error`, default 10%): detections with
  `flux_err / flux` at or above this are dropped before matching, in every
  frame — a safety net against reporting "variability" from measurements
  too imprecise to support the claim.
- **Photometric normalization** (`normalize`, default `true`): a uniform
  per-frame flux-scale offset (see [`photometric_scale`](@ref) — on real
  ZTF data this tracked each frame's own seeing, a fixed-aperture-radius
  effect, not a code defect) is corrected before the chi-squared test, so
  a constant star doesn't read as variable just because one exposure's
  aperture happened to enclose a different fraction of its PSF. Set
  `normalize=false` to skip this (e.g. if `flux` is already calibrated).

For each detection in frame 1 (after the S/N floor), the closest detection
within `position_tolerance` pixels of that *same* position in every other
frame (after its own S/N floor) is matched. Groups reaching at least
`min_frames` matched points are tested against [`variability_chi2`](@ref)'s
constant-flux null hypothesis; a group is kept only if its **reduced**
chi-squared (`chi2 / dof`) exceeds `chi2_threshold`.

`chi2_threshold`'s default is deliberately not close to the naive "3" a
formally-correct chi-squared test would suggest. Measured on real ZTF
field 451 (119 matched, high-S/N stationary stars, after normalization):
the reduced-chi2 distribution has a real, heavy tail — 23% of ordinary
stars exceed a threshold of 3, and even a threshold of 20 still flags
13% — almost certainly because `detect_sources` finds a peak-detected
integer pixel each frame, not a sub-pixel centroid, so a 1 px jitter
against a small `aperture_radius` (default 3 px, comparable to a typical
seeing FWHM) shows up as a real, if spurious, flux change. `50.0` cuts
that to 8% on the same data — still elevated above the true stellar
variable fraction, not a clean cut, and honestly reported as such rather
than tuned to look better than it is: treat any candidate here as
requiring independent confirmation (a VSX/SIMBAD match via
[`crossmatch_catalog`](@ref), or a period recovered by
[`recover_rotation_period`](@ref)), not as self-evidently real. Properly
fixing the underlying cause — forced, sub-pixel-centroided photometry
instead of peak-position aperture photometry — is future work, not
attempted here.

With the default `min_frames` (every frame must match at high S/N), a
field where few sources clear the S/N floor in every frame will correctly
return few or no candidates rather than a list built on unreliable
photometry — lower `max_relative_error` or `min_frames` deliberately if a
shallower search is wanted, rather than reading an empty result as "no
code path found here".

What the chi2 test is discriminating against depends on how
`detections_per_frame` was produced:

- **ZOGY path** (`detections_per_frame` from `S_corr` — `run_pipeline` or
  `search_field` given a `reference`): a constant star differences to
  zero and is never detected in `S_corr` at all, so a matched stationary
  group is already a variable/transient by construction; the chi2 test
  here is a second filter against the known subtraction-residual failure
  mode at bright stars (see `INVESTIGATION_LOG.md`), not the primary
  discriminant. `normalize`/`photometric_scale` do not apply meaningfully
  here, since `S_corr` isn't on a physical flux scale to begin with —
  pass `normalize=false` on this path.
- **Raw path** (no `reference`): every star above `threshold` is detected
  in every frame, so the chi2 test (with normalization and the S/N floor)
  is the only thing separating a genuine variable from an ordinary
  constant star.

Returns a vector of matched-point groups, each a vector of
`(frame, x, y, flux, flux_err)` named tuples sorted by frame — the same
shape [`link_candidates`](@ref) returns (with two extra fields, and with
`flux`/`flux_err` already normalized when `normalize=true`), so a group
can be passed directly to [`astrometric_calibrate`](@ref).

For periodicity, build `(times, flux)` from a returned group's own
`frame`/`flux` fields and `timestamps`, and pass to
[`recover_rotation_period`](@ref) directly.
"""
function find_variable_sources(detections_per_frame, timestamps;
                                position_tolerance::Real=2.0,
                                min_frames::Integer=length(detections_per_frame),
                                chi2_threshold::Real=50.0,
                                normalize::Bool=true,
                                max_relative_error::Real=0.10)
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

        chi2, dof = variability_chi2([p.flux for p in group], [p.flux_err for p in group])
        dof > 0 && chi2 / dof > chi2_threshold && push!(groups, group)
    end

    return groups
end
