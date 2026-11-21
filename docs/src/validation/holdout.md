# Sets kept aside

## Summary

Every setting of the search was chosen, or checked, on the same 11
Pan-STARRS1 sets of IASC campaigns. Twelve other sets, of IASC's October
2026 campaign, were kept aside on 2026-11-19, before they were used for any
measurement of a setting, and the measurements to repeat on them were
declared then, with what counts as agreeing. The first run failed through
the measuring tool: [`inject_movers`](@ref) placed synthetic asteroids up to
1″ off where the headers were a few arcseconds off, on these sets and, less,
on the 11. Fixed, and run again (no longer a clean test, as the first had been
seen): the dipole test kept 532 of 533 synthetic asteroids and every
catalogued one, as on the 11 sets; half were found at V 21.57 frame by
frame (21.49–21.67), against 21.72 there, not met, and at V 21.67 with
stacking (21.53–21.86), against 21.84, met. These sets are shallower than
the 11, by a median 0.1 mag at their 5σ limit. The next campaign's sets are
the clean test.

## 1. The sets

The 12 Pan-STARRS1 sets (4 × 45 s each) of IASC's October 2026 campaign that
no measurement of a setting used, as IASC sent them:

| Night | Sets |
|:--|:--|
| 2026-10-12 | `XY02_p00`, `XY75_p01`, `XY75_p10` |
| 2026-10-14 | `XY04_p01`, `XY04_p10`, `XY04_p11`, `XY05_p00`, `XY05_p01` |
| 2026-10-19 | `XY55_p00`, `XY55_p01`, `XY55_p10`, `XY55_p11` |

The other 7 sets of that campaign (`XY45_p01`, `XY45_p10`, `XY45_p11`,
`XY72_p01`, `XY72_p11`, `XY73_p00`, `XY73_p01`) are among the 13 the
settings were measured on. Before this page, the 12 were searched only with
[`analyze_set`](@ref)'s defaults, for [the ITF](itf.md) and
[digest2](digest2.md), and three of their candidates were looked at by eye;
nothing was injected into them.

## 2. Rules

- Nothing in the pipeline is changed on what these sets show. A fault they
  reveal is reported here; if it is fixed, the fix is measured on other
  data, and these sets are not used to choose it.
- Each measurement declared is run once and reported in full.
- The sets of IASC's next campaign will be added to those kept aside, listed
  here before any of them is opened.

## 3. Declared

In `examples/holdout_test.jl`. Expected, from the 11 sets the settings were
measured on, with 95% intervals (Wilson's for fractions, resampling whole
sets for magnitudes; `examples/validation_stats.jl`):

| | There | Agrees if |
|:--|:--|:--|
| Synthetic asteroids found with the dipole test (6σ), over those found without it | 511 of 516, 99.0% | 99.0% is inside this interval |
| Other tracklets the test leaves (neither synthetic nor catalogued) | 19 of 89, 21% | 21% is inside this interval; not measurable with fewer than 10 without the test |
| Catalogued asteroids identified | the same with and without the test | the same |
| Half the synthetic asteroids found, frame by frame | V 21.67 (21.49–21.85) | 21.67 is inside this interval |
| The same, with stacking | V 21.79 (21.59–21.95) | 21.79 is inside this interval |

The searches are those of [Real or false](real-bogus.md) (100 synthetic
asteroids per set, V 20–23, no stacking) and [How faint](faint.md) (100
per set, V 21–23.5, the defaults), with other seeds.

## 4. First run: not met, through the measuring tool

Run on 2026-11-19 as declared. Two sets were left out by the rules:
[`inject_movers`](@ref) found no photometric zero point in any frame of
`XY04_p11` and `XY05_p00`.

| | Expected | Measured on 10 sets | |
|:--|:--|:--|:--|
| synthetic asteroids found with the dipole test, over without | 99.0% | 308 of 308, 100% (99–100%) | met |
| other tracklets the test leaves | 21% | 42 of 42, 100% (92–100%) | not met |
| catalogued asteroids identified | the same | the same 6 | met |
| half found, frame by frame | V 21.67 | V 21.30 (21.14–21.65) | not met |
| half found, with stacking | V 21.79 | V 21.30 (21.14–21.76) | not met |

| Set | Synthetic found, of 100 (V 20–23) | Other tracklets | Half found, frame by frame |
|:--|--:|--:|:--|
| `XY02_p00` | 17 | 1 | V 21.12 |
| `XY04_p01` | 15 | 15 | never half |
| `XY04_p10` | 19 | 8 | V 21.44 |
| `XY05_p01` | 17 | 12 | never half |
| `XY75_p01` | 12 | 5 | V 21.17 |
| `XY75_p10` | 37 | 1 | V 21.57 |
| `XY55_p00` | 41 | 0 | V 21.71 |
| `XY55_p01` | 49 | 0 | V 21.48 |
| `XY55_p10` | 53 | 0 | V 21.88 |
| `XY55_p11` | 48 | 0 | V 21.62 |

The four sets of 2026-10-19 behaved as the 11 sets had (41–53 found, no
other tracklets); the six of 10-12 and 10-14 found 12–37 and left 1–15 other
tracklets, which the dipole test did not touch.

**Why.** On those six sets the headers' WCS is 4–33″ off. [`analyze_set`](@ref)
fits the frames together ([`refine_wcs_set`](@ref)) and matched 207–338 Gaia
stars at 0.12″; [`inject_movers`](@ref) fitted each frame on its own, which
matched 5–12 stars and landed about 1″ from the joint fit, or failed: 33″ is
past one frame's search, and those two sets had no zero point. The synthetic
asteroids were placed through the wrong fit, up to 1″ from where the search
measures them, against a recovery radius of 1.5″: found ones were not
counted, and counted as other tracklets instead. What failed was the
measuring tool, not the search. The same comparison on the 13 sets the
settings were measured on found it there too, less: 0.8–1.05″ off on 5 of
their frames, and no fit on 1 (`XY72_p11`), so the measurements that inject
asteroids were biased on those sets as well.

## 5. Changed after the first run

[`inject_movers`](@ref) now fits the frames together, as
[`analyze_set`](@ref) does, checked on a synthetic set with headers 25″ off,
where a frame's own fit fails. The pipeline is unchanged. Declared before
running, in the scripts' headers: the measurements that inject asteroids on
the 11 sets are run again with the same seeds ([How faint](faint.md),
[Real or false](real-bogus.md), [Uncertainties](uncertainties.md)), and these
12 sets a second time, against the values the 11 give with the same
injector. That second run is not a clean test of sets kept aside: the first
was seen, and the tool was fixed on what it showed. The clean test is the
next campaign's sets.

## 6. Second run: the injection fixed

Run on 2026-11-20, all 12 sets: with the frames fitted together, every frame
of `XY04_p11` and `XY05_p00` had a zero point. Expected are the values the
11 sets give with the same injector ([Real or false](real-bogus.md),
[How faint](faint.md)), with the same rules for agreeing:

| | Expected | Measured on 12 sets | |
|:--|:--|:--|:--|
| synthetic asteroids found with the dipole test, over without | 583 of 589, 99.0% | 532 of 533, 100% (99–100%) | met |
| other tracklets the test leaves | 12 of 80, 15% | 3 of 3 | not measurable |
| catalogued asteroids identified | the same | the same 10 | met |
| half found, frame by frame | V 21.72 | V 21.57 (21.49–21.67) | not met |
| half found, with stacking | V 21.84 | V 21.67 (21.53–21.86) | met |

| Set | Synthetic found, of 100 (V 20–23) | Other tracklets | 5σ limit | Half found, frame by frame | with stacking |
|:--|--:|--:|:--|:--|:--|
| `XY02_p00` | 42 | 0 | 22.40 | V 21.48 | V 21.48 |
| `XY04_p01` | 45 | 0 | 22.53 | V 21.56 | V 21.56 |
| `XY04_p10` | 38 | 0 | 22.46 | V 21.43 | V 21.46 |
| `XY04_p11` | 36 | 0 | 22.43 | V 21.21 | V 21.21 |
| `XY05_p00` | 44 | 0 | 22.57 | V 21.67 | V 21.81 |
| `XY05_p01` | 47 | 0 | 22.55 | V 21.71 | V 21.76 |
| `XY55_p00` | 49 | 1 | 22.63 | V 21.68 | V 21.70 |
| `XY55_p01` | 51 | 0 | 22.60 | V 21.68 | V 21.88 |
| `XY55_p10` | 57 | 1 | 22.59 | V 21.70 | V 22.05 |
| `XY55_p11` | 41 | 0 | 22.60 | V 21.55 | V 21.59 |
| `XY75_p01` | 42 | 1 | 22.48 | V 21.72 | V 21.90 |
| `XY75_p10` | 41 | 0 | 22.49 | V 21.43 | V 22.02 |

"Synthetic found" and "other" are without the dipole test; the last three
columns are of the second search (V 21–23.5, the defaults).

| V | Injected | Frame by frame | With stacking |
|:--|:--|:--|:--|
| 21.00–21.25 | 117 | 91 (78%) | 91 (78%) |
| 21.25–21.50 | 124 | 87 (70%) | 87 (70%) |
| 21.50–21.75 | 109 | 49 (45%) | 57 (52%) |
| 21.75–22.00 | 123 | 20 (16%) | 50 (41%) |
| 22.00–22.25 | 114 | 2 (2%) | 30 (26%) |
| 22.25–22.50 | 129 | 0 (0%) | 24 (19%) |
| 22.50–22.75 | 119 | 0 (0%) | 6 (5%) |
| 22.75–23.50 | 365 | 0 | 0 |

The six sets that had shown the fault found 38–47 synthetic asteroids, and
1 other tracklet between them: those of the first run were the injector's.
Frame by frame, half were found 0.15 mag brighter than on the 11 sets.
These 12 sets are shallower: their 5σ limits are V 22.40–22.63 (median
22.54), those of the 11 V 22.23–23.12 (median 22.64), and none is as deep as
the four crowded fields (V 22.81–23.12), which have the faintest
half-found magnitudes of the 11. Half were found a median 0.94 mag brighter
than a set's limit here, 0.89 there. Frame by frame, no asteroid fainter than
V 22.25 was found here (0 of 613), against 15 of 568 on the 11 sets. What
makes the difference was not measured; nothing is changed on it, by the
rules.
