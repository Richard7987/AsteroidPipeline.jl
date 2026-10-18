# Pan-STARRS1 image sets

## Summary

On ten real Pan-STARRS1 image sets, the pipeline's positions agree with an
independent reduction to 0.13–0.19″; it recovers 17 of the 30 catalogued
asteroids predicted on the frames (every miss explained) and an uncatalogued
one confirmed by hand, finds ~80% of injected asteroids down to V≈21.5, and
leaves 107 tracklets to check across the nine good sets. In a
blind test on four sets, it missed no real object, and each of its false
candidates led to a filter that now removes that kind.

## 1. Introduction

A search pipeline has to do three things well: measure positions precisely
enough for the MPC, find the moving objects that are there, and not bury a
person in false ones. This report measures all three on images of one
telescope, Pan-STARRS1, whose short image sets (four exposures of one field
over about an hour) are distributed to the public for asteroid searches by
the International Astronomical Search Collaboration (IASC).

## 2. Data

| Sets | Taken | Use |
|:--|:--|:--|
| 5 (`XY14_p10`, `XY15_p01`, `XY25_p10`, `XY26_p01`, `XY42_p11`) | 2019 | catalogued asteroids |
| 1 (`XY54_p10`) | 2025-09-16 | an independent reduction (Astrometrica) with an uncatalogued object; completeness |
| 4 (`XY72_p01`, `XY72_p11`, `XY73_p00`, `XY73_p01`) | 2026-10-04 | blind test |

Each set is four 45 s exposures of one 10′ × 10′ chip (about 2,430 × 2,430
pixels of 0.256″), about 16 minutes apart. The images have gaps between detector cells
filled with a constant, `BLANK` pixels, and header astrometry several
arcseconds off.

## 3. Methods

**Pipeline**: [`analyze_set`](@ref) with [`pipeline_config`](@ref) at the
images' pixel scale — astrometry refitted to Gaia DR3, median-differenced
images searched at 6σ, the [artifact filters](../concepts/artifacts.md), and
tracklets of three or more frames ([How it works](../concepts/pipeline.md)).

**Catalogued asteroids**: predicted by SkyBoT from the observatory (code
F51) at each frame's mid-exposure time. One counts as recovered if a single
tracklet lies within 3″ of its prediction in at least three frames.

**Astrometry**: the pipeline's report for one object compared, line by line,
with a report made independently in Astrometrica (also calibrated against
Gaia).

**Completeness**: 180 synthetic asteroids (V 18–23, 5–100″/h) and 80 fast
ones (100–1000″/h) injected into `XY54_p10`'s real frames
([Measure completeness](../guides/completeness.md)).

**Blind test**: the four 2026 sets were analyzed before being checked by
eye; every candidate was then examined in the images and against the MPC's
orbits.

## 4. Results

### Astrometry

| Frame | Pipeline (RA, Dec) | Astrometrica (RA, Dec) | Separation |
|:--|:--|:--|--:|
| 1 | 00 12 25.816 −10 43 50.18 | 00 12 25.824 −10 43 50.33 | 0.19″ |
| 2 | 00 12 25.279 −10 43 54.80 | 00 12 25.286 −10 43 54.70 | 0.14″ |
| 3 | 00 12 24.750 −10 43 59.11 | 00 12 24.753 −10 43 58.99 | 0.13″ |
| 4 | 00 12 24.206 −10 44 03.48 | 00 12 24.215 −10 44 03.43 | 0.14″ |

![Pipeline minus Astrometrica, four frames](../assets/astrometry-vs-astrometrica.png)

The header astrometry alone was 6.8″ off; refitted to Gaia stars it has
0.07–0.08″ residuals. Getting the time right mattered as much: the headers'
`MJD-OBS` is in TAI (37 s after the UTC `DATE-OBS`), and MPC reports need
the middle of the exposure, not its start (22.5 s later). Magnitudes
differed by 0.1–0.9 mag, at signal-to-noise 4–8 per frame.

### Catalogued asteroids

| Configuration | Sets | Tracklets to check | Objects recovered |
|:--|--:|--:|:--|
| raw frames, 8σ, header astrometry | 6 | 5,548 | 3 of 6 |
| difference images, refitted astrometry | 6 | 157 | 6 of 6 |
| + artifact filters | 10 | 125 | 14 catalogued + 1 uncatalogued |
| + line test scaled to the astrometry ([`pipeline_config`](@ref)) | 9 good + 1 damaged | 107 + 84 | 17 catalogued + 1 uncatalogued |

The six objects of the first two rows are the detectable ones in those sets:
2009 SG135, 2015 XJ232, 2019 PK5, 2019 NB9 (a Jupiter Trojan), 2018 LT and
the uncatalogued one. The artifact filters cost no catalogued asteroid
(14 recovered with or without them) and removed 8 false tracklets.

**A bright asteroid missed — and why.** (13014) Hasslacher, V 17.4, crossed
`XY72_p01` and was never reported, though it is one of the brightest objects
in the set. It was detected in three frames (the fourth puts it in a cell
gap) and linked — then dropped by the straight-line test: its three points
lay 1.5 px from a line, and partial tracklets had to be within 1 px. The
asteroid wasn't off: the frames' astrometry was. That field has few Gaia
stars (9–14 matched per frame) and its solutions have 0.08–0.23″ residuals,
so no line could be checked to 1 px (0.26″). The test now allows three times
the frames' astrometric residual, at most 2 px. Hasslacher and 2017 DY33
(V 21.7) are found; across the nine good sets four more tracklets appear,
two of them these asteroids. On the damaged set, which [`analyze_set`](@ref)
doesn't search, the count rises from 22 to 84.

Of the 13 catalogued objects still missed, 5 are in the damaged set, 7 are
fainter than V≈21 (the depth of these images), and one, 4311 T-1, crosses a
cell gap in one frame and a star in the next. Two more that a first count
placed in the field are 30–90 px off the chip once the astrometry is
corrected.

### Completeness

![Fraction of injected asteroids recovered](../assets/completeness.png)

| Magnitude | 18–21.5 | 21.5–22 | 22–22.5 | 22.5–23 |
|:--|--:|--:|--:|--:|
| recovered | 79% (92/117) | 72% | 16% | 0% |

| Rate (″/h) | 5–10 | 10–100 | 100–400 | 400–1000 |
|:--|--:|--:|--:|--:|
| recovered | 50% | 61% | 43% | 18% |

### A blind test

| Set | Checked by hand | The pipeline's candidates were |
|:--|:--|:--|
| `XY72_p01` | nothing new | 2019 UA185 (a known object SkyBoT misplaced); a hot pixel |
| `XY72_p11` | nothing new | a ghost of a bright star |
| `XY73_p00` | images badly damaged | 26 tracklets of detector defects |
| `XY73_p01` | nothing new | a diffraction spike |

No real object was missed. Each false candidate was traced to its cause and
measured against the real objects; four changes followed:

| False candidate | What separated it | Change |
|:--|:--|:--|
| 2019 UA185 | SkyBoT's orbit was outdated; the MPC's file of that day put it within 4″ | identification against the MPC's orbit file ([`read_mpcorb`](@ref)) |
| hot pixel | 0 px of movement on the detector, against ≥32 px for real objects | `detector_artifact_radius` |
| ghost, spike | 1.72σ and 2.44σ of background excess, against −0.07 to 0.12σ | `max_background_excess` |
| damaged set | 137–1126 sources per reference star, against 0.6–5.4 | [`set_quality`](@ref) |

![Set quality on the ten sets](../assets/set-quality.png)

Re-run with the changes, the four sets gave no false candidate.

## 5. Discussion

The positions are as good as an established manual reduction's, which is
what the MPC needs. Completeness falls off between V 21.5 and 22.5, the
depth of a single 45 s exposure; going deeper needs a reference built from
other nights ([Deep references](../guides/difference-imaging.md)). Slow
(≤10″/h) and fast (≥100″/h) objects are found less often, for the reasons in
[Limitations](../concepts/limitations.md).

The blind test is the most informative result: on images no one had looked
at, every false candidate had a physical cause that a measurable property
separated from real objects. The filters built from them are specific to
those causes, not tuned to these sets: they removed no catalogued asteroid
in any of the ten sets, and 2 of 180 injected ones. They are, however, measured on one telescope; on another,
check them with [an injection test](../guides/completeness.md) and a few
known objects.

Ten sets are few. The numbers here describe these images; the methods —
catalogued objects followed frame by frame, injection, a blind test — are
what to repeat on others.
