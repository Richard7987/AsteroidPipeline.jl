# Sets kept aside

## Summary

Every setting of the search was chosen, or checked, on the same 11
Pan-STARRS1 sets of IASC campaigns. Twelve other sets, of IASC's October
2026 campaign, were kept aside on 2026-11-19, before they were used for any
measurement of a setting, and the measurements to repeat on them were
declared then, with what counts as agreeing. This page holds that
declaration and the results.

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
