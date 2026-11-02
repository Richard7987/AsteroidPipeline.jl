# The light curve of an asteroid

## Summary

[`moving_light_curve`](@ref) followed the Jupiter trojan (624) Hektor through
107 public ZTF exposures over six months. Within a night, the magnitudes scatter
as much as their reported errors. [`fourier_period`](@ref) gives
3.46017 ± 0.00004 h, which doubled is 6.92035 ± 0.00008 h, against the published
6.920509 h. These points alone can't tell the period from half of it: the
doubling is an assumption, not a measurement.

## 1. Introduction

The pipeline finds moving objects; a rotation period asks more of it: that one
object be measured on every frame, on the same scale, from one night to the
next. Hektor has a precisely known period, from shape models fitted to light
curves taken from 1957 to 1991
([Marchis et al. 2014](https://doi.org/10.1088/2041-8205/783/2/L37)), and it is
bright enough for ZTF to measure well.

The period was known before this test, so this isn't a blind one. What was
declared before running it: the photometry within a night should scatter as much
as its reported errors, and, since ZTF visits an object only a few times a night,
the period should fit better than random ones, not come out of the periodogram
by itself.

## 2. Data

| | |
|:--|:--|
| Target | (624) Hektor, published period 6.920509 ± 0.000002 h |
| Images | ZTF, r band, July–December 2022: 107 of 112 exposures available, 600″ cutouts over 49 nights |
| Positions | predicted for each exposure by IRSA's moving object search, then recentred on the frame |
| Used | 85 measurements over 168 days; the others not measured, or marked `star_nearby` |

## 3. Methods

`examples/moving_light_curve_demo.jl` fetches the images and runs every step:
[`moving_light_curve`](@ref) night by night, [`reduced_magnitudes`](@ref) with
[`object_geometry`](@ref) to correct for distance and phase angle, then
[`recover_rotation_period`](@ref) and [`fourier_period`](@ref) between 2 and 20 h.

## 4. Results

Within a night, 43 points scatter 0.009 mag about a straight line, consistent
with reported errors of 0.003–0.012 mag.

| Period search | Period |
|:--|:--|
| periodogram's highest peak | 3.460 h |
| Fourier series, 2 harmonics | 3.46017 ± 0.00004 h |
| doubled | **6.92035 ± 0.00008 h** |
| published | 6.920509 ± 0.000002 h |

The doubled period is 0.00016 h (0.6 s) below the published one, about 2σ. The
amplitude of the fit is 0.51 mag.

![The measurements folded on the published period, with a Fourier series of two harmonics](../assets/hektor-light-curve.png)

**Period or half of it?** Fits to the same data around each:

| Harmonics | RMS near 6.92 h | RMS near 3.46 h |
|:--|--:|--:|
| 2 | 0.0524 mag | 0.0454 mag |
| 4 | 0.0440 mag | 0.0426 mag |

Half the period fits as well or better, so [`recover_rotation_period`](@ref)
keeps 3.46 h: its test finds no significant difference between alternate cycles.

**Against other periods**, with 2 harmonics:

| Data | Published period | Daily alias, 8.0864 h | Random periods, median | Random periods fitting better |
|:--|--:|--:|--:|--:|
| July–December, 85 points | 0.051 mag | 0.076 mag | 0.162 mag | 0% |
| July–August, 48 points | 0.047 mag | 0.030 mag | 0.146 mag | 2.8% |

## 5. Discussion

Six months of sparse points pin the period down to about one part in 10⁵, and
two months don't: over July and August, the daily alias fits better than the
published period. A period from a survey like this needs the long baseline.

Three limits:

- **The doubling is assumed.** An elongated asteroid has two maxima per rotation,
  so the period is taken as twice the fit's. Telling the two apart needs a
  night or more of dense coverage.
- **The fit leaves 0.045 mag**, five times the scatter within a night. A fixed
  curve doesn't describe six months: most likely Hektor is seen from a changing
  angle, which isn't modelled here.
- **One object, with a period known beforehand.** Asteroids of unknown period,
  and other telescopes, are not tested here.
