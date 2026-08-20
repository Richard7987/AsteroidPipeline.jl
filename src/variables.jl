"""
    variability_chi2(flux, flux_err; systematic_error_fraction::Real=0.03) -> (chi2, dof)

Chi-squared goodness of fit of `flux` against the constant-flux
(non-variable) null hypothesis, using the error-weighted mean as the
constant-flux estimate. `dof` is `length(flux) - 1` (one free parameter,
the mean).

The per-point error used is not `flux_err` alone, but
`sqrt(flux_err^2 + (systematic_error_fraction * flux)^2)` — a systematic
error floor, proportional to flux, added in quadrature. Without it, a
bright star's tiny *formal* (Poisson/background) error makes this test
wildly oversensitive to real but non-astrophysical systematics (imperfect
flat-fielding, frame-to-frame PSF variation, ...) that never show up in
`flux_err` at all. This was found, not assumed: on real ZTF data (field
451), the false positives this chi2 test produced were *not* the
faintest, noisiest stars, as a naive read of the S/N floor would suggest
— they were the **brightest** ones, with the smallest formal errors (a
median 0.25% relative error among flagged stars vs. 2.93% among the
rest), the textbook signature of a systematic floor rather than
underestimated statistical noise. A floor of 1% (a standard value in
forced-photometry pipelines, tried first) cut the false-positive rate at
a reduced-chi2 threshold of 3 from 22% (no floor) to 6%, and at
threshold 20 from 13% to essentially 0%. Sweeping the floor further
against the same real, matched stationary stars — not guessed — found
more room, and was later cross-checked against two more real,
independently-confirmed variable stars in two more real fields, not just
one: 3% (the current default) keeps every one of three real confirmed
variables above `chi2_threshold=10.0` (reduced chi2 13.8/61.2/17.5 across
the three fields) while cutting the false-positive rate on real
stationary stars in those same three fields to 1.3%/2.04%/0.55% (vs.
2.0%/2.55%/0.55% at the previous 2% default — better or equal in every
field, never worse) — see [`find_variable_sources`](@ref)'s docstring
for the full before/after table across all three. A literal 0% false
positive rate is achievable (a 5% floor reaches it on two of the three
fields), but only by also erasing two of the three real confirmed
variables' own signal below `chi2_threshold` — not a real improvement,
just a different failure mode.

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
                           chi2_threshold::Real=10.0,
                           normalize::Bool=true,
                           max_relative_error::Real=0.10,
                           systematic_error_fraction::Real=0.03)

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

Three corrections are applied before any variability test:

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
- **Systematic error floor** (`systematic_error_fraction`, default 3%,
  forwarded to [`variability_chi2`](@ref)): the actual dominant fix for
  this function's real, measured false-positive floor — see below.

For each detection in frame 1 (after the S/N floor), the closest detection
within `position_tolerance` pixels of that *same* position in every other
frame (after its own S/N floor) is matched. Groups reaching at least
`min_frames` matched points are tested against [`variability_chi2`](@ref)'s
constant-flux null hypothesis; a group is kept only if its **reduced**
chi-squared (`chi2 / dof`) exceeds `chi2_threshold`.

`chi2_threshold`'s default and `systematic_error_fraction` (forwarded to
[`variability_chi2`](@ref)) were both calibrated against real ZTF data —
first against field 451 (119 matched, high-S/N stationary stars), later
cross-checked against two more independent real fields — and the
calibration story is worth reading, because the first hypothesis tried
here was wrong. `detect_sources` used to center its aperture on
`PeakMesh`'s raw *integer*-pixel peak; the working theory was that
per-frame pixel-grid jitter in that peak, against a small
`aperture_radius`, was producing spurious flux swings. Fixing that
(`detect_sources` now refines to a sub-pixel centroid — see its own
docstring) barely moved the false-positive rate at all (23% → 22% at a
threshold of 3; 13% → 13% at 20) — a real, measured non-result, not
swept under the rug. Checking *which* stars were actually being flagged
settled it: not the faintest ones, as low-S/N selection would predict,
but the **brightest** — median 0.25% relative error among flagged stars
vs. 2.93% among the rest, the signature of a systematic error floor
(imperfect flat-fielding, frame-to-frame PSF variation — nothing
`flux_err`'s pure Poisson/background model captures) rather than
underestimated statistical noise. Adding that floor
(`systematic_error_fraction`, in `variability_chi2`) is what actually
fixed it: a 1% floor cut it to 6% at a threshold of 3, ~0% at 20; a real,
large improvement over the 23-8% range measured before either fix.

Sweeping the floor further, first against field 451 alone, found more
room without giving up real sensitivity — checked against a real
confirmed variable (ASASSN-V J183620.31, from a different real field),
not just the stationary-star side of the tradeoff. That first pass
settled on 2%. Later, two more real, independently-confirmed variables
were found in two more real ZTF fields (V1012 Mon and ASASSN-V
J072906.85-090518.2 — dense fields near the galactic plane, each with
over 180 real matched stationary stars of its own), letting the same
sweep be repeated three times over, each with its own real confirmed
variable and its own real false-positive sample, instead of trusting one
field's result to generalize. That repeat is what moved the default from
2% to 3%: at 3%, all **three** confirmed variables (reduced chi2
13.8/61.2/17.5 across field 451/V1012 Mon/ASASSN-072906) still clear
`chi2_threshold=10.0`, while the false-positive rate in every one of the
three fields is the same or better than at 2% (1.3%/2.04%/0.55% at 3% vs.
2.0%/2.55%/0.55% at 2%). A literal 0% false-positive rate is reachable —
a 5% floor gets there on two of the three fields — but only by also
pushing two of the three confirmed variables' own reduced chi2 below
`chi2_threshold`, i.e. by going blind to real variability rather than
genuinely improving anything; 3% is the largest floor across all three
real fields that doesn't cost that. Still not zero: treat any candidate
here as requiring independent confirmation (a VSX/SIMBAD match via
[`crossmatch_catalog`](@ref), or a period recovered by
[`recover_rotation_period`](@ref)), not as self-evidently real.

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
  mode at bright stars (see the [Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#PSF-timing-and-astrometric-noise-bugs-behind-ZOGY's-excess-detections)), not the primary
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
