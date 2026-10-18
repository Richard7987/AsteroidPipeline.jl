# ZTF field 451

## Summary

On one night of real Zwicky Transient Facility (ZTF) images, the pipeline
recovers both catalogued asteroids in the field, with and without
difference imaging against a deep reference (ZOGY). Running it on real data
found five errors no synthetic test had, all fixed; and a frame degraded by
thin cloud is recognised and left out.

## 1. Introduction

This was the pipeline's first test on real images, and its first use of
ZOGY difference imaging ([Zackay, Ofek & Gal-Yam 2016](https://doi.org/10.3847/0004-637X/830/1/27)),
which subtracts a deep reference image of the static sky and can find
objects fainter than a single exposure shows. The questions: does each stage
work on real data, and does ZOGY behave as its model says?

## 2. Data

| | |
|:--|:--|
| Field | ZTF field 451, CCD 1, quadrant 1, r band |
| Science images | 5 × 30 s, 2019-10-23, over 6.25 h |
| Reference | 30 images of the same field from other nights |
| Pixel scale | 1.01″ |
| Catalogued asteroids | 2002 UY45 and 1997 KO3, both main belt |

The data are public; `examples/fetch_data.sh` downloads them and
`examples/real_data_demo.jl` runs the test.

## 3. Methods

Two runs on the same five frames: undifferenced (detection on the raw
frames, 8σ) and ZOGY ([`build_reference`](@ref), [`zogy_subtract`](@ref)).
Tracklets were cross-matched against SkyBoT. ZOGY's significance map,
`S_corr`, should have unit standard deviation; that was checked per frame.

## 4. Results

| Run | Tracklets | Catalogued asteroids recovered |
|:--|--:|:--|
| undifferenced | 133 | 2 of 2 |
| ZOGY | 667 | 2 of 2, at the same positions |

Five errors surfaced and were fixed, each with a regression test:

| Error | Effect on real data |
|:--|:--|
| Julian Dates sent to SkyBoT in scientific notation | every cross-match silently empty |
| the images not background-subtracted before the FFT | >99.9% of a frame "above 6σ" |
| one NaN at a reprojected frame's border | an entire `S_corr` map NaN |
| that border filled with 0 instead of the sky level | a biased background in the difference |
| the PSF measured before reprojection, used after | residuals at every bright star: `std(S_corr)` 2.0–2.6 |

With them fixed, `std(S_corr)` came to 1.1–1.2 in four of the five frames.
The fifth stayed at 2.0: one-sided excess light around stars, confined to
one part of the detector, with a flat background — thin cloud during the
exposure. `quality_max_std=1.5` leaves it out; both asteroids are still
recovered.

## 5. Discussion

Both runs recover both asteroids, but this field doesn't test ZOGY's real
advantage: both objects are bright enough to find without it. The larger
tracklet count with ZOGY isn't a measure of its noise either: leaving out the
clouded frame meant accepting 4-of-5-frame tracklets, which alone doubled
the count (334 → 667). A fair comparison needs fainter objects and the same
frames in both runs.

The errors found are the main result: each was invisible on synthetic data,
and two of them made the output silently wrong rather than failing. Real
data has to be part of testing.
