# Real or false: a classifier and the dipole test

## Summary

On 11 Pan-STARRS1 sets of IASC campaigns with 100 synthetic asteroids added
to each, a logistic regression trained to tell real tracklets from false
ones, validated on whole sets it had not seen, did no better than the
search's own filters: 534 real tracklets kept and 67 false ones let through,
against 540 and 89. Its weights pointed at one measurement the filters did
not use, a negative beside the detection in the difference image. Added as a
test at the detection threshold, it took the tracklets that are neither
synthetic nor catalogued from 89 to 19; 511 of the 516 synthetic asteroids
found before were still found, and all 25 catalogued ones.

## 1. Method

**Labels.** No set of tracklets labelled by eye is large enough, so the
labels come from known truth. [`inject_movers`](@ref) adds synthetic
asteroids to a set; a search with the tests that decide what is real relaxed
(no flux-ratio, bright-background or detector-artifact test, a straight-line
residual of up to 3 px) lets false tracklets through. A tracklet is real if
[`injection_recovery`](@ref) matches it to a synthetic asteroid or
[`analyze_set`](@ref) identifies it as a catalogued one, false otherwise.

**Features.** Ten per tracklet, read from the frames' normalized differences
([`stack_difference`](@ref)): the number of points, the significance at
each (minimum and median), the ratio of the largest to the smallest flux,
the straight-line residual, the rate, the fraction of masked pixels, two
measures of shape, and the deepest negative pixel within 4 px of a point.

**Classifier.** Logistic regression on the standardized features, with a
ridge penalty. One whole set is left out at a time, so nothing from its
field is seen in training. The threshold keeps, on the training sets, as
many real tracklets as the search's filters keep there.

**The dipole test.** A moving object is absent from the sequence's median,
so in a frame's difference it is only a positive source. A star that a
frame's misalignment or seeing subtracts badly leaves a positive and a
negative lobe side by side. [`run_pipeline`](@ref)'s `max_negative` drops a
detection with an unmasked pixel that many σ below zero within 4 px.

## 2. Data

Declared before running, in `examples/real_bogus_dataset.jl`,
`examples/real_bogus_train.jl` and `examples/dipole_test.jl`:

| | |
|:--|:--|
| Sets | the 13 local Pan-STARRS1 sets of IASC campaigns; the two [`analyze_set`](@ref) marks unreliable left out |
| Synthetic asteroids | 100 per set, V 20.0–23.0, 5–100″/h, other seeds than [How faint](faint.md)'s |
| Searches | frame by frame, no stacking, observatory F51, the MPC's orbits |
| Dipole test | `max_negative` at 8 and at 6, the detection threshold, against the search without it |

Two things were changed before any result was seen and one after, all three
in the dataset script's header: relaxing detection as well made linking run
for half an hour on one set without finishing; the relaxed search through
[`analyze_set`](@ref) stalled asking the catalogue services about hundreds
of tracklets; and matching a relaxed tracklet to a catalogued asteroid by a
single point labelled 36 tracklets as catalogued on a set with 4, so 3
points within 1″ are required.

## 3. Results

### The classifier

4,022 tracklets: 566 real (542 synthetic, 24 catalogued) and 3,456 false,
3,362 of the false on the three crowded fields.

| Set left out | Real | False | Filters: real kept | false passed | Classifier: real kept | false passed |
|:--|--:|--:|--:|--:|--:|--:|
| `XY15_p01` | 31 | 20 | 30 | 1 | 25 | 0 |
| `XY25_p10` | 82 | 2,296 | 76 | 49 | 79 | 29 |
| `XY26_p01` | 62 | 329 | 55 | 12 | 58 | 7 |
| `XY42_p11` | 63 | 737 | 62 | 20 | 59 | 10 |
| `XY45_p01` | 46 | 4 | 46 | 0 | 46 | 1 |
| `XY45_p10` | 49 | 4 | 48 | 0 | 49 | 0 |
| `XY45_p11` | 51 | 12 | 50 | 1 | 50 | 4 |
| `XY54_p10` | 53 | 6 | 51 | 1 | 52 | 2 |
| `XY72_p01` | 54 | 12 | 52 | 2 | 45 | 1 |
| `XY72_p11` | 32 | 22 | 31 | 1 | 30 | 7 |
| `XY73_p01` | 43 | 14 | 39 | 2 | 41 | 6 |
| All | 566 | 3,456 | 540 | 89 | 534 | 67 |

On the three crowded fields the classifier passed 46 false tracklets against
the filters' 81 and kept 196 real ones against 193; on the other eight it
passed 21 against 8 and kept 338 against 347. Fitted on every set, the
largest weights were the deepest negative (+2.51 per standard deviation),
the straight-line residual (−2.45) and the number of points (+2.12); the
others were between −0.66 and +0.38.

### The dipole test

The same sets and synthetic asteroids, searched with the filters as they
were and with the test added:

| Set | Without: synthetic | other | At 8σ: synthetic | other | At 6σ: synthetic | other |
|:--|--:|--:|--:|--:|--:|--:|
| `XY15_p01` | 30 | 1 | 30 | 1 | 30 | 1 |
| `XY25_p10` | 72 | 49 | 72 | 23 | 70 | 5 |
| `XY26_p01` | 54 | 12 | 53 | 5 | 52 | 1 |
| `XY42_p11` | 60 | 20 | 59 | 16 | 59 | 5 |
| `XY45_p01` | 44 | 0 | 44 | 0 | 44 | 0 |
| `XY45_p10` | 46 | 0 | 46 | 0 | 46 | 0 |
| `XY45_p11` | 49 | 1 | 49 | 1 | 49 | 1 |
| `XY54_p10` | 50 | 1 | 50 | 1 | 50 | 1 |
| `XY72_p01` | 45 | 2 | 45 | 2 | 45 | 2 |
| `XY72_p11` | 30 | 1 | 30 | 1 | 30 | 1 |
| `XY73_p01` | 36 | 2 | 36 | 2 | 36 | 2 |
| All | 516 | 89 | 514 | 52 | 511 | 19 |

"Other" is every tracklet neither synthetic nor catalogued. The 25
catalogued asteroids were identified in all three searches. The eight
uncrowded fields did not change: 8 other tracklets in each search.
[`pipeline_config`](@ref) now sets the test at 6σ.

### With stacking, and the other validations

[`stack_search`](@ref) takes the same `max_negative`, for a track with such a
pixel within 4 px in any frame, and [`analyze_set`](@ref) passes it the
configuration's. The other validations on Pan-STARRS1 were then repeated with
the test in both searches:

| | Before | With the test |
|:--|:--|:--|
| [digest2](digest2.md) and [ITF](itf.md), 19 campaign sets: catalogued asteroids identified | 26 | the same 26, with the same scores |
| unidentified tracklets | 4 | the same 4 |
| [How faint](faint.md), 11 sets, 1,100 synthetic asteroids: found frame by frame | 299 | 295 |
| found with stacking | 358 | 354 |
| 50% found at, frame by frame and with stacking | V 21.68, 21.81 | V 21.67, 21.79 |
| found at V 22.0–22.5 with stacking | 34 of 214 | 34 of 214 |
| stacked tracklets not synthetic | 6 | 5 |
| [Stacking](stacking.md), the 2 of its 4 fields run both ways, 160 synthetic asteroids: found frame by frame, then with stacking | 35, 41 | 35, 41 |
| [Pan-STARRS1 sets](pan-starrs1.md): 260 synthetic asteroids on one set, found | 140 | 140 |
| 9 good sets, frame by frame: catalogued asteroids identified | 19 | 19 |
| tracklets to look at | 98 | 31 |
| 20 sets, frame by frame: catalogued asteroids identified | 28 | 28 |
| tracklets not identified | 79 | 12 |

The 4 asteroids no longer found are between V 21.25 and 22.0, on two of the
crowded fields. The stacked tracklet no longer added was the one on
`XY25_p10`, at V 23.3 and barely visible.

## 4. Discussion

**The classifier.** With these features a linear model is no gain over the
filters: fewer false tracklets on crowded fields, more on clean ones, and 6
fewer real ones overall. Its use here was to show which measurement the
filters lacked.

**The cost of the test.** The 5 synthetic asteroids lost at 6σ, 1% of
those found, were all on the crowded fields, where a real object is more
often within 4 px of a star's residual. An asteroid passing a star is lost
in that frame with or without the test; with it, also when the star is
merely close.

**What wasn't measured.** Real tracklets are almost all synthetic: 24
catalogued among 566. An uncatalogued real asteroid counts as false, or as
"other"; the production search found 6 unidentified tracklets on 13 sets,
so few. The 6σ and 8σ values were compared on the same 11 sets that suggested
the test, all of one survey; no other telescope was tried. On the one other survey tried, [ZTF field 451](ztf.md)
(5 frames, 1.01″ pixels), [`analyze_set`](@ref) gave the same two catalogued
asteroids and no other tracklet with the test as without it; no crowded
field of another survey was tried.

## 5. Reproduce

```bash
julia --project=. examples/real_bogus_dataset.jl
julia --project=. examples/real_bogus_train.jl
julia --project=. examples/dipole_test.jl
```
