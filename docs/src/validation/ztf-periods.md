# Rotation periods of eight asteroids

## Summary

The test of [The light curve of an asteroid](hektor.md), on eight more asteroids
chosen before looking at their data. Of the eight, two could not be measured
(one too bright for ZTF, one with too few exposures). Of the six measured, two
periods came out within 1% of the published ones, (588) Achilles and (3451) Mentor;
four did not. With ZTF's few visits a night, a period from one season of survey
data is a candidate to check, not a measurement.

## 1. Declared before running

The objects and their periods were taken from the Asteroid Lightcurve Database
(LCDB, [Warner et al. 2009](https://doi.org/10.1016/j.icarus.2009.02.003),
revision of October 2023, as JPL's Small-Body Database gives it), without
looking at their ZTF data:

| | |
|:--|:--|
| Objects | four Jupiter trojans, two main-belt and two near-Earth asteroids, with periods of 4–11 h |
| Images | every ZTF r exposure of 2022 in which the object is predicted at V 13.5–19: brighter saturates, fainter is noise |
| Too few | fewer than 30 measurements |
| Recovered | the Fourier period, doubled (two maxima per rotation), within 1% of the published one |

Every object is reported, recovered or not.

## 2. Methods

`examples/ztf_periods.jl` fetches the images from IRSA and runs the same steps
as for Hektor: [`moving_light_curve`](@ref) night by night,
[`reduced_magnitudes`](@ref) with [`object_geometry`](@ref), and
[`fourier_period`](@ref) with two harmonics between 2 and 20 h. IRSA refused
most downloads made in parallel (HTTP 503): they are made one at a time.

## 3. Results

| Object | Published (h) | Exposures | Measurements, nights | Doubled period (h) | Off | Result |
|:--|:--|:--|:--|:--|:--|:--|
| (588) Achilles, trojan | 7.306 | 88 | 78, 53 | 7.3060 ± 0.0001 | 0.00% | recovered |
| (3451) Mentor, trojan | 7.702 | 83 | 80, 62 | 7.6970 ± 0.0001 | 0.06% | recovered |
| (884) Priamus, trojan | 6.8605 | 106 | 95, 68 | 12.0032 ± 0.0007 | 75% | not recovered |
| (911) Agamemnon, trojan | 6.592 | 132 | 106, 63 | 23.9457 ± 0.0019 | 263% | not recovered |
| (1627) Ivar, near-Earth | 4.795 | 124 | 120, 35 | 7.9123 ± 0.0002 | 65% | not recovered |
| (1036) Ganymed, near-Earth | 10.297 | 31 | 31, 22 | 6.3308 ± 0.0001 | 39% | not recovered |
| (243) Ida, main belt | 4.634 | 21 | 21 | — | — | too few |
| (216) Kleopatra, main belt | 5.385 | 0 | — | — | — | too bright: V 10.0–12.3 all year |

Exposures are those available; IRSA had no image for 61 of the predicted ones.
The errors are the fit's formal ones, smaller than the real uncertainty.

## 4. Discussion

**The failures.** In frequency, Priamus's fitted period lies 2.998 cycles a day
from half the published one: an alias of the once-a-night sampling, which the
data can't tell from the true period. Agamemnon's, 11.97 h before doubling, is
half a day, again the sampling rather than the asteroid. For Ivar (3.94 cycles
a day) and Ganymed (2.92) the offsets are near whole numbers but not equal to
them; Ganymed has only 31 measurements.

**Recovered, but how precisely?** Mentor's published period has three decimals.
Over the 2022 season, about 1,100 rotations, a difference of 0.005 h adds up
to most of a rotation, so the published period itself folds these data
poorly: 91% of random periods fold them better. The measured 7.6970 h may be
closer to the true period than the catalogue's value, or not: these data can't
say which.

**What would help.** Several visits a night, or a dedicated light curve over a
few nights, as for Hektor's published period; or combining ZTF's g and r.

## 5. Conclusion

Two of six. A period from ZTF survey data alone is a candidate, to be checked
for daily aliases and confirmed with a dedicated light curve.
