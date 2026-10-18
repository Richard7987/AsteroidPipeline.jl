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

## Results at a glance

| Question | Answer | Measured on |
|:--|:--|:--|
| Are the positions right? | within 0.13–0.19″ of an independent reduction | Pan-STARRS1, 4 frames |
| Are catalogued asteroids found? | 17 of the 30 predicted on the frames; of the 13 missed, 5 are in the one damaged set, 7 are fainter than V≈21, and one crosses a cell gap and a star | 10 Pan-STARRS1 sets |
| How faint? | ~80% found to V≈21.5, 16% at 22–22.5 | 180 injected asteroids |
| How many false candidates? | 107 tracklets to look at across the 9 good sets (5,548 for 6 at first) | Pan-STARRS1 |
| Are damaged sets recognised? | the one damaged set flagged; none of the other nine | Pan-STARRS1 |
| Variable stars? | all 3 confirmed variables found; 0.6–2% of constant stars flagged | ZTF |

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
`variable_star_demo.jl` (variable stars).
