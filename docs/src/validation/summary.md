# Validation summary

How well does the pipeline work? Each claim below was tested on real
telescope data against something independent of the pipeline: catalogued
asteroids' predicted positions, another program's measurements of the same
images, or synthetic asteroids added to real images.

## Data

| Data set | Telescope | Images | Ground truth | Report |
|:--|:--|:--|:--|:--|
| 10 image sets, 2019–2026 | Pan-STARRS1 (1.8 m, 0.26″/px) | 4 × 45 s per set | 14 catalogued asteroids; an independent Astrometrica reduction; 260 injected asteroids | [Pan-STARRS1](pan-starrs1.md) |
| field 451, one night | ZTF (1.2 m Schmidt, 1.0″/px) | 5 + 30 reference | 2 catalogued asteroids | [ZTF](ztf.md) |
| 3 fields with known variables | ZTF | one night each | 3 catalogued variable stars | [A variable star](variable-star.md) |
| (624) Hektor, 6 months | ZTF | 107 exposures | published rotation period | [The light curve of an asteroid](hektor.md) |
| 4 fields, 80 injected asteroids each | Pan-STARRS1 | 4 × 45 s per set | 320 injected asteroids | [Fainter objects by stacking](stacking.md) |
| 15 known asteroids, 40 arcs since 2024 | MPC astrometry, every station | 3 nights to 30 days | JPL's orbits | [Preliminary orbits](orbits.md) |
| 8 asteroids, 2022 | ZTF | 21–132 exposures each | published rotation periods, chosen beforehand | [Rotation periods of eight asteroids](ztf-periods.md) |
| 30 asteroids chosen by a rule, 2018–2024 | ZTF | 32–687 exposures each | published rotation periods; a shuffled control | [Rotation periods of thirty asteroids](period-sample.md) |

## Results at a glance

| Question | Answer | Measured on |
|:--|:--|:--|
| Are the positions right? | within 0.12–0.17″ of an independent reduction | Pan-STARRS1, 4 frames |
| Are catalogued asteroids found? | 17 of the 30 predicted on the frames; of the 13 missed, 5 are in the one damaged set, 7 are fainter than V≈21, and one crosses a cell gap and a star; with the frames' astrometry fitted together, 28 identified on 20 sets (24 before) | 10, then 20 Pan-STARRS1 sets |
| How faint? | ~82% found to V≈21.5, 21% at 22–22.5 frame by frame; with stacking, 101 of 320 injected at V 21–24 instead of 82 | 180 + 320 injected asteroids |
| How many false candidates? | 107 tracklets to look at across the 9 good sets (5,548 for 6 at first); 91 not identified across 20 sets (150 before the frames were fitted together) | Pan-STARRS1 |
| Are damaged sets recognised? | the one damaged set flagged; none of the other nine | Pan-STARRS1 |
| Variable stars? | all 3 confirmed variables found; 0.6–2% of constant stars flagged | ZTF |
| Orbits? | semi-major axis within a median 4.8% from 3 nights, 0.03% from 30 days; 39 of 40 predictions 30 days ahead within 2σ | MPC astrometry |
| Rotation periods? | 6.92035 ± 0.00008 h for a published 6.920509 h, doubling the fit, on one asteroid; on eight chosen beforehand, two of the six measured within 1%, four not; on thirty chosen by a rule, the period or its half within 1% for 13 of 29, against 3 of 29 with the magnitudes shuffled | ZTF |

![Tracklets a person has to look at, by configuration](../assets/tracklets-to-vet.png)

## How the tests are run

- **Catalogued asteroids** count as recovered only if one tracklet follows
  the object's predicted position — from the observatory, at each frame's
  own time — in at least three frames. Matching a tracklet's first point
  alone counted static stars near an asteroid's path as recoveries: 7 of 9
  in an early check were stars.
- **Injected asteroids** are added to the real frames with each frame's own
  PSF, zero point, noise and exposure trailing, and counted the same way.
- **A blind test**: four image sets were analyzed before anyone looked at
  them, then checked by hand ([Pan-STARRS1](pan-starrs1.md#A-blind-test)).

The scripts are in `examples/`: `iasc_demo.jl` and `analyze_sets.jl` (image
sets), `injection_test.jl` (completeness), `real_data_demo.jl` (ZTF),
`variable_star_demo.jl` (variable stars), `moving_light_curve_demo.jl`, `ztf_periods.jl` and `rotation_period_sample.jl` (rotation periods), `stacking_test.jl` (stacking), `orbit_test.jl` (orbits).
