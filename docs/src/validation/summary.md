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
| 11 sets of IASC campaigns, 100 injected asteroids each | Pan-STARRS1 | 4 × 45 s per set | 1,100 injected asteroids, V 21–23.5 | [How faint, on campaign sets](faint.md) |
| 15 known asteroids, 40 arcs since 2024 | MPC astrometry, every station | 3 nights to 30 days | JPL's orbits | [Preliminary orbits](orbits.md) |
| 73 asteroids near opposition, 12 nights of 2025 | MPC astrometry, one station (P07) | 323 tracklets, then 969 false ones added | the objects' designations | [Linking across nights](linking.md) |
| 19 sets of IASC's October 2026 campaign | Pan-STARRS1 | 4 × 45 s per set | the MPC's Isolated Tracklet File | [Already reported: the ITF](itf.md) |
| the same sets | Pan-STARRS1 | 4 × 45 s per set | the MPC's orbits of the asteroids identified | [Near-Earth or not: digest2](digest2.md) |
| 11 sets of IASC campaigns, 100 injected asteroids each | Pan-STARRS1 | 4 × 45 s per set | 1,100 injected asteroids, V 20–23, and the catalogued ones | [Real or false](real-bogus.md) |
| 11 sets of IASC campaigns, as they are and with 60 injected asteroids each | Pan-STARRS1 | 4 × 45 s per set | JPL's ephemerides of 25 catalogued asteroids; 660 injected asteroids | [Uncertainties](uncertainties.md) |
| 8 asteroids, 2022 | ZTF | 21–132 exposures each | published rotation periods, chosen beforehand | [Rotation periods of eight asteroids](ztf-periods.md) |
| 30 asteroids chosen by a rule, 2018–2024 | ZTF | 32–687 exposures each | published rotation periods; a shuffled control | [Rotation periods of thirty asteroids](period-sample.md) |
| 24 sets of minor planets (of 33 tried), raw and as the archive reduced them, 2024 | Las Cumbres Observatory (0.4 m, 1 m, 2 m MuSCAT) | 5 per set | the archive's own reduction (BANZAI); JPL's ephemerides | [Raw frames, reduced here](calibration.md) |
| 15 sets, nothing set for the instrument, 100 injected asteroids each | ZTF (1.01″/px), Las Cumbres Observatory (0.4 m, 1 m), DECam (Blanco 4 m, 0.263″/px) | 4–5 per set | 1,500 injected asteroids; JPL's ephemerides of the catalogued ones | [Other telescopes](telescopes.md) |

## Results at a glance

| Question | Answer | Measured on |
|:--|:--|:--|
| Are the positions right? | within 0.12–0.17″ of an independent reduction | Pan-STARRS1, 4 frames |
| Are catalogued asteroids found? | 17 of the 30 predicted on the frames; of the 13 missed, 5 are in the one damaged set, 7 are fainter than V≈21, and one crosses a cell gap and a star; with the frames' astrometry fitted together, 28 identified on 20 sets (24 before) | 10, then 20 Pan-STARRS1 sets |
| How faint? | ~82% found to V≈21.5, 21% at 22–22.5 frame by frame; with stacking, 94 of 320 injected at V 21–24 instead of 77; on 11 campaign sets, half found to V 21.7 frame by frame and V 21.8 with stacking, 23% at 22–22.5 with stacking | 180 + 320 + 1,100 injected asteroids |
| How many false candidates? | 31 tracklets to look at across the 9 good sets (5,548 for 6 at first, 98 before the dipole test); 12 not identified across 20 sets (150 before the frames were fitted together, 79 before the dipole test) | Pan-STARRS1 |
| Are damaged sets recognised? | the one damaged set flagged; none of the other nine | Pan-STARRS1 |
| Variable stars? | all 3 confirmed variables found; 0.6–2% of constant stars flagged | ZTF |
| Orbits? | semi-major axis within a median 4.8% from 3 nights, 0.03% from 30 days; 39 of 40 predictions 30 days ahead within 2σ | MPC astrometry |
| Tracklets linked across nights? | every tracklet of the 47 objects seen on 3 nights or more, none wrong, with or without 969 false tracklets; on 2 nights, 168 wrong among false tracklets | MPC astrometry |
| Were the candidates reported already? | 2 of the 5 unidentified on 19 campaign sets, by Pan-STARRS itself in the same exposures, 0.10–0.21″ from ours | Pan-STARRS1, the MPC's ITF |
| Near-Earth by digest2? | 25 of 26 catalogued main-belt and Mars-crossing asteroids below 65, the NEO Confirmation Page's threshold; one Mars-crosser at 98 | Pan-STARRS1, 18 campaign sets |
| False tracklets? | 80 tracklets neither injected nor catalogued on 11 sets, 12 with the dipole test, for 6 of 589 injected asteroids lost (89, 19 and 5 of 516 before the injection was fixed); on 12 sets kept aside, 532 of 533 kept; a classifier trained on the tracklets did no better than the filters | Pan-STARRS1, 11 campaign sets |
| Are the uncertainties right? | positions: for 12 numbered asteroids, errors over uncertainties with deviation 1.19 (robustly 1.01), 91% within 2, no offset (−0.007″, −0.019″); on injected asteroids larger than the errors, which carry no error of the astrometric fit (1.5 times too small before the injection was fixed); magnitudes, with what lies still taken out, deviation 1.06 (robustly 0.97; 2.35 before) | Pan-STARRS1, 11 campaign sets |
| Raw frames reduced here as well as the observatory's? | on 13 sets: the same astrometric fits as its pixels' from the same headers, zero points as uniform (a median 0.000 mag apart), limiting magnitudes within 0.1, positions against JPL as close (0.070″ against 0.102″); 14 of its 18 catalogued asteroids, the 4 missed 2 at the threshold, 1 with a frame since fixed, 1 identified only on a second search | Las Cumbres Observatory, 3 cameras |
| Other telescopes? | half found frame by frame 0.82–1.18 mag brighter than the frames' limit (0.89 on Pan-STARRS1), catalogued as often as injected, position uncertainties right on injected asteroids; against JPL right on DECam and Las Cumbres, 1.67 times too small on ZTF; after two new tests and whole Gaia answers, 0.3 unidentified tracklets per set on ZTF and 1.0 on Las Cumbres, but 128 on one DECam set in the galactic bulge, which found 2 of 100 injected | ZTF, Las Cumbres, DECam: 15 sets |
| Rotation periods? | 6.92035 ± 0.00008 h for a published 6.920509 h, doubling the fit, on one asteroid; on eight chosen beforehand, two of the six measured within 1%, four not; on thirty chosen by a rule, the period or its half within 1% for 13 of 29, against 3 of 29 with the magnitudes shuffled | ZTF |

![Tracklets a person has to look at, by configuration](../assets/tracklets-to-vet.png)

## How sure

Most results are counts from a few dozen objects. Their 95% intervals
(Wilson's; `examples/validation_stats.jl`):

| Measured | Count | 95% interval |
|:--|:--|:--|
| catalogued asteroids predicted on the frames, found ([Pan-STARRS1](pan-starrs1.md)) | 17 of 30 | 39–73% |
| synthetic asteroids V 18–21.5 found ([Pan-STARRS1](pan-starrs1.md)) | 96 of 117 | 74–88% |
| synthetic asteroids V 21–24 found, frame by frame and with stacking ([Stacking](stacking.md)) | 77 and 94 of 320 | 20–29% and 25–35% |
| synthetic asteroids V 22.0–22.5 found, frame by frame and with stacking ([How faint](faint.md)) | 26 and 53 of 231 | 8–16% and 18–29% |
| half the synthetic asteroids found, frame by frame and with stacking ([How faint](faint.md)) | V 21.72 and 21.84 | V 21.60–21.80 and 21.71–22.04 (sets resampled) |
| confirmed variable stars found ([A variable star](variable-star.md)) | 3 of 3 | 44–100% |
| predictions 30 days ahead within 2σ ([Orbits](orbits.md)) | 39 of 40 | 87–100% |
| objects on 3 nights or more linked whole, none wrong ([Linking](linking.md)) | 47 of 47 | 92–100% |
| unidentified tracklets already in the ITF ([ITF](itf.md)) | 2 of 5 | 12–77% |
| catalogued asteroids, not near-Earth, below digest2's 65 ([digest2](digest2.md)) | 25 of 26 | 81–99% |
| synthetic asteroids still found with the dipole test ([Real or false](real-bogus.md)) | 583 of 589 | 98–100% |
| the same, on 12 sets kept aside ([Sets kept aside](holdout.md)) | 532 of 533 | 99–100% |
| other tracklets left by the dipole test ([Real or false](real-bogus.md)) | 12 of 80 | 9–24% |
| half the synthetic asteroids found on 12 sets kept aside, frame by frame and with stacking ([Sets kept aside](holdout.md)) | V 21.57 and 21.67 | V 21.49–21.67 and 21.53–21.86 (sets resampled) |
| points of numbered asteroids within 2 uncertainties ([Uncertainties](uncertainties.md)) | 41 of 45 | 79–96%, too narrow: 12 asteroids, not 45 independent points |
| catalogued asteroids of the archive's frames found in raw frames reduced here ([Raw frames](calibration.md)) | 14 of 18 | 55–91% |
| rotation periods within 1%, of those measured ([Eight asteroids](ztf-periods.md)) | 2 of 6 | 10–70% |
| period or its half within 1%; magnitudes shuffled ([Thirty asteroids](period-sample.md)) | 13 of 29; 3 of 29 | 28–62%; 4–26% |

Synthetic asteroids added to the same few fields share those fields' stars
and seeing, so these intervals are narrower than a new field would show;
where sets were resampled whole, they are not. Every setting was chosen on
the same 11 Pan-STARRS1 sets; [Sets kept aside](holdout.md) measures again on
others.

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
