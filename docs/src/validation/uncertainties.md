# Uncertainties of positions and magnitudes

## Summary

[`analyze_set`](@ref) gives each measured point an uncertainty in position
(`pos_err`, arcsec per axis) and in magnitude (`mag_err`), and
[`ades_psv`](@ref) reports them to the MPC as `rmsRA`, `rmsDec` and `rmsMag`.
Against JPL's ephemerides of the 12 numbered asteroids identified on 11
Pan-STARRS1 campaign sets (45 points), the position errors divided by the
uncertainties have a standard deviation of 1.19 (robustly 1.01), 64% of them
within 1 and 91% within 2, and no offset: −0.007″ in RA, −0.019″ in Dec. On
synthetic asteroids the magnitude uncertainties are about right (robust
deviation 1.2); the position ones first came out 1.5 times too small with a
mean offset of 0.05″, which were the test's: placed through the same fit as
the search's, the errors are smaller than the uncertainties and the offset
is gone. A first model of the position uncertainty was 4 times too small and
was changed; every run is reported below.

## 1. What is given

**Magnitude.** The aperture's noise, from the sky's scatter and the source's
own photons through the header's `GAIN` (1 without one), with the zero
point's: the scatter of its Gaia stars over the square root of their number.

**Position, per axis.** Two parts in quadrature. The Cramér–Rao bound for
the centroid of a Gaussian of the frame's FWHM on that noise, what an ideal
centroid could reach. And the frame's astrometric fit: the scatter of its
Gaia stars about it, per axis. The fit is affine, so the field's distortion
stays in its residuals, and a target's position carries it as a star's does.
Where only an offset was fitted the uncertainty is unknown (`NaN`, an empty
field in the ADES report).

## 2. Synthetic asteroids

**Design**, declared before running (`examples/uncertainty_test.jl`): the 13
Pan-STARRS1 sets of IASC campaigns in `data/real/iasc`, 60 synthetic
asteroids added to each ([`inject_movers`](@ref), V 19.5–22.0, 5–100″/h),
searched by [`analyze_set`](@ref) with its defaults. Every point of a
tracklet matched to a synthetic asteroid, against its true position and
magnitude in that frame. If the uncertainties are right, the errors divided
by them have a standard deviation of 1, 68% within 1 and 95% within 2. Two
sets marked unreliable by `analyze_set` are left out.

**First run** (seeds 4000 + the set's place): the astrometric part was the
fit's error of its mean, its scatter over the square root of the number of
stars. That is far too small:

| | Points | Error, median | Uncertainty, median | Deviation (robust) | Within 1 / 2 |
|:--|--:|--:|--:|--:|--:|
| Position | 1,899 | 0.073″ | 0.018″ | 22.8 (5.5) | 17% / 32% |
| Magnitude | 1,899 | 0.042 | 0.048 | 2.89 (1.19) | 56% / 85% |

**Second run**, with the astrometric part changed to the stars' scatter, on
new seeds (5000 + the set's place) so that the model is not checked on the
data that changed it:

| | Points | Error, median | Uncertainty, median | Deviation (robust) | Within 1 / 2 |
|:--|--:|--:|--:|--:|--:|
| Position, all | 1,897 | 0.076″ | 0.072″ | 4.60 (1.51) | 49% / 71% |
| V 19.5–20.5 | 806 | 0.071″ | 0.070″ | 4.78 (1.44) | 51% / 71% |
| V 20.5–21.5 | 798 | 0.078″ | 0.074″ | 4.63 (1.52) | 49% / 71% |
| V 21.5–22.0 | 293 | 0.094″ | 0.080″ | 3.95 (1.65) | 46% / 70% |
| stacked tracklets | 92 | 0.114″ | 0.098″ | 4.06 (1.77) | 42% / 67% |
| Magnitude, all | 1,897 | 0.045 | 0.049 | 2.00 (1.23) | 54% / 84% |
| V 19.5–20.5 | 806 | 0.028 | 0.029 | 1.35 (1.17) | 49% / 81% |
| V 21.5–22.0 | 293 | 0.095 | 0.119 | 2.53 (1.21) | 58% / 89% |

The position errors had a mean of +0.051″ in RA and −0.047″ in Dec. The
synthetic asteroids are placed through each frame's own astrometric fit,
while `analyze_set` measures through the frames' fit made together, so the
offset could be the test's or the pipeline's. The deviations well above the
robust ones come from a few points far off, likely ones blended with stars.

**Third run**, after [`inject_movers`](@ref) came to fit the frames together
as `analyze_set` does (found on [the sets kept aside](holdout.md): its own
fits were up to 1″ from the joint one on 5 frames of these sets). Declared
before it: if the two fits' disagreement was the cause, a deviation closer to
1 than 1.5 and a smaller offset. Same seeds; the asteroids drawn differ, as
each one's noise depends on its frame's zero point.

| | Points | Error, median | Uncertainty, median | Deviation (robust) | Within 1 / 2 |
|:--|--:|--:|--:|--:|--:|
| Position, all | 1,968 | 0.026″ | 0.072″ | 0.95 (0.53) | 87% / 98% |
| V 19.5–20.5 | 829 | 0.019″ | 0.070″ | 0.84 (0.39) | 94% / 99% |
| V 20.5–21.5 | 827 | 0.032″ | 0.073″ | 0.97 (0.64) | 85% / 98% |
| V 21.5–22.0 | 312 | 0.047″ | 0.079″ | 1.12 (0.86) | 73% / 95% |
| stacked tracklets | 104 | 0.055″ | 0.098″ | 0.89 (0.79) | 81% / 97% |
| Magnitude, all | 1,968 | 0.044 | 0.051 | 2.35 (1.24) | 53% / 84% |

The mean error is −0.002″ in RA and +0.001″ in Dec. The offset and the
excess were the test's. The position uncertainties are now larger than the
errors on synthetic asteroids, as they should be: a synthetic asteroid is
placed through the same fit it is measured through, so it carries none of
the fit's error, which the uncertainty includes and a real asteroid carries
(section 3).

## 3. Catalogued asteroids

**Design**, declared before running (`examples/known_astrometry_test.jl`):
the same 13 sets as they are, nothing added, searched by `analyze_set` with
its defaults. Every point of every tracklet identified as a catalogued
asteroid, against JPL Horizons' astrometric position of that asteroid from
F51 at the point's epoch. Horizons' own uncertainty is left out, which for a
numbered asteroid is milliarcseconds; points of unnumbered ones reported
apart. No offset means the mean errors over the sets within 2 standard
errors of 0. One set is marked unreliable; one other has no identified
asteroid.

| | Points | Objects | Error, median | Uncertainty, median | Deviation (robust) | Within 1 / 2 | Mean error RA, Dec |
|:--|--:|--:|--:|--:|--:|--:|:--|
| numbered | 45 | 12 | 0.055″ | 0.080″ | 1.19 (1.01) | 64% / 91% | −0.007″, −0.019″ |
| unnumbered | 47 | 13 | 0.057″ | 0.067″ | 1.63 (1.27) | 56% / 86% | +0.057″, −0.003″ |

One unnumbered asteroid, 2019 UA185, is left out of the table: its 3 points
are some 10° from JPL's position, whose orbit rests on 11 observations over 8 days
in 2019 (the MPC's, which identified it, is a later one). Of the rest, 2016
SX110's 4 points are 0.46″ off in RA together, the way an orbit's error
shows. Over the 11 sets' means the error is +0.043 ± 0.044″ in RA and
−0.012 ± 0.012″ in Dec: no offset.

## 4. What it means

On catalogued numbered asteroids, where the truth is independent of the
pipeline, the positions have no offset at this level and their
uncertainties are the right size: the reported `rmsRA` and `rmsDec` can be
taken as they are. The sample is small, 45 points of 12 asteroids on one
telescope. The offset and the larger errors first measured on synthetic asteroids were
the test's: with the frames fitted together for both, they are gone. The magnitude uncertainties are about right for most points and too
small for a few: a magnitude off by more than 3 of its uncertainties is
more likely a blend than noise.
