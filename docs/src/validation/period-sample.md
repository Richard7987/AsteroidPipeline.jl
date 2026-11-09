# Rotation periods of thirty asteroids

## Summary

The test of [Rotation periods of eight asteroids](ztf-periods.md), on thirty
asteroids chosen by a rule fixed before measuring any of them: ten near-Earth
asteroids, ten Jupiter trojans and ten main-belt asteroids from the Asteroid
Lightcurve Database, measured on public ZTF images with the method used for
(624) Hektor and nothing tuned. Twenty-nine could be measured. The fitted period
came out within 1% of the published one, or of half of it, for **13 of the 29**;
4 fell on aliases of the daily sampling and 12 elsewhere. With each object's
magnitudes shuffled among its epochs, 3 of the 29 land within the same 1% by
chance. The recoveries are well above that rate, but a period from one season
of survey data remains a candidate to check, and which of a period and its half
is the true one, these data don't say.

## 1. Declared before measuring

Written in `examples/rotation_period_sample.jl` before any object was measured:

| | |
|:--|:--|
| Catalogue | Asteroid Lightcurve Database (LCDB, [Warner et al. 2009](https://doi.org/10.1016/j.icarus.2009.02.003)), current public release |
| Eligible | numbered; period quality U 3 or 3−; period 2–20 h; largest amplitude ≥ 0.3 mag; not from a survey, not binary, not tumbling; (624) Hektor left out |
| Groups | LCDB family codes: near-Earth (9101), Jupiter trojan (9202), main belt (9104, 9105, 9106: inner, middle, outer) |
| Selection | each group shuffled with a fixed seed (20261105); the first 10 with at least 60 ZTF r exposures at predicted V 14–19 within 180 days. Only the number of exposures is looked at, never a measurement |
| Images | those exposures, in the 180-day window with the most of them |
| Recovered | the fitted period, 2–20 h with two harmonics, within 1% of the published period or of half of it. Twice the period and daily aliases (1 or 2 cycles a day from the period or its half) are reported as such and don't count |
| Control | each object's magnitudes shuffled among its epochs and searched the same way |
| Reported | every selected object, measured or not |

Half the period counts because a light curve with two maxima per rotation, the
usual case, is fitted best at half the period, as for Hektor.

## 2. Sample

| Group | Eligible in the LCDB | Examined to find 10 | Taken |
|:--|:--|:--|:--|
| near-Earth | 281 | 60 | 10 |
| Jupiter trojan | 62 | 11 | 10 |
| main belt | 1,246 | 13 | 10 |

Three biases follow from the rule:

- **The main belt is its background.** Codes 9104–9106 are the LCDB's
  background populations; members of collisional families (Flora alone has 142
  eligible) have their own codes and were not in the pool. This was set when the
  groups were written, before measuring; it doesn't bear on the method, but the
  ten aren't a sample of the whole main belt.
- **Near-Earth asteroids long in view.** 50 of the 60 examined had fewer than 60
  exposures in 180 days: the ten taken are those ZTF saw bright for months.
- **Easier periods.** U 3 and an amplitude of at least 0.3 mag select
  well-determined periods of elongated bodies: easier than an unknown asteroid.

## 3. Methods

`examples/rotation_period_sample.jl` selects the sample, fetches the images and
measures every object:

1. **Exposures.** IRSA's moving object search lists the ZTF exposures that
   contain the object. Queried over the whole 2018–2024 survey, it returns
   nothing for near-Earth asteroids after several minutes, so it is queried only
   over the spans in which JPL Horizons predicts the object within half a
   magnitude of V 14–19.
2. **Images.** A 600″ cutout of each exposure in the window, fetched one at a
   time: IRSA answers parallel requests, and at times sequential ones, with
   HTTP 503, and these are retried with a growing wait.
3. **Photometry.** [`moving_light_curve`](@ref) night by night, as for Hektor;
   measurements not measured or marked `star_nearby` are dropped.
4. **Period.** [`reduced_magnitudes`](@ref) with [`object_geometry`](@ref) from
   Palomar (I41), then [`fourier_period`](@ref) with two harmonics between 2 and
   20 h, on the measurements and on the shuffled control.

**What was lost, and redone.**

- The first run fetched the images without retrying 503s: 2,762 cutouts were
  lost to it, and some objects were fitted on a fraction of their frames. Its
  results were set aside and everything was run again with the retries; in that
  run no image was lost to a 503.
- 1,092 of the 4,639 exposures in the windows have no image at IRSA (HTTP 404,
  the same on every attempt). They are not spread evenly: (2329) Orthos has 44
  of 184, (11509) Thersilochos 56 of 246, (2893) Peiroos 51 of 208, (782)
  Montefiore 52 of 190, (2674) Pandarus 32 of 165.
- For (1572) Posnania, the Gaia archive answered HTTP 429 (too many requests) on
  8 nights, which went unmeasured. It was measured again once those requests
  went through: 71 measurements gave 16.2297 h (twice the period); all 117 give
  4.8387 h, a daily alias. The second is the one counted.

## 4. Results

Exposures are those available at IRSA; measurements, those used for the fit.
The span is between the first and last measurement, in days.

| Group | Object | Published (h) | U | Amplitude (mag) | Exposures | Measurements | Span (d) | Fitted (h) | Result | Control |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| near-Earth | (6569) Ondaatje | 5.958 | 3− | 1.12 | 90 | 82 | 121 | 2.9792 | P/2 | other |
| near-Earth | (2329) Orthos | 4.6473 | 3− | 0.67 | 44 | 23 | 180 | 7.8655 | other | daily alias |
| near-Earth | (5332) Davidaguilar | 5.803 | 3 | 1.2 | 57 | 51 | 173 | 5.8017 | P | other |
| near-Earth | (347813) 2002 NP1 | 5.915 | 3 | 0.86 | 56 | 25 | 171 | 3.4332 | other | other |
| near-Earth | (6053) 1993 BW3 | 2.573 | 3 | 0.5 | 104 | 88 | 101 | 9.5916 | other | other |
| near-Earth | (53435) 1999 VM40 | 5.189 | 3 | 0.38 | 178 | 165 | 177 | 2.5927 | P/2 | other |
| near-Earth | (11405) 1999 CV3 | 6.504 | 3 | 0.89 | 72 | 69 | 55 | 2.0539 | other | other |
| near-Earth | (495615) 2015 PQ291 | 13.7504 | 3 | 0.89 | 51 | 43 | 50 | 6.8748 | P/2 | other |
| near-Earth | (5189) 1990 UQ | 6.653 | 3 | 1.02 | 66 | 56 | 75 | 2.9217 | daily alias | other |
| near-Earth | (2100) Ra-Shalom | 19.797 | 3 | 0.61 | 94 | 76 | 131 | 11.9358 | other | other |
| trojan | (4834) Thoas | 18.192 | 3 | 0.39 | 106 | 94 | 180 | 3.1575 | other | other |
| trojan | (3793) Leonteus | 5.621 | 3 | 0.31 | 102 | 98 | 176 | 12.0074 | other | other |
| trojan | (11509) Thersilochos | 17.367 | 3 | 0.38 | 56 | 34 | 177 | 6.2054 | other | other |
| trojan | (2893) Peiroos | 8.945 | 3 | 0.35 | 51 | 28 | 149 | 4.4743 | P/2 | other |
| trojan | (23135) Pheidas | 8.69 | 3− | 0.31 | 128 | 119 | 178 | 5.1100 | other | P |
| trojan | (2674) Pandarus | 8.48 | 3 | 0.65 | 32 | 2 | — | — | not measured | — |
| trojan | (5476) Mulius | 5.78 | 3 | 0.31 | 79 | 49 | 178 | 2.8896 | P/2 | other |
| trojan | (17314) Aisakos | 9.67 | 3 | 0.34 | 42 | 30 | 161 | 4.0419 | daily alias | other |
| trojan | (9030) Othryoneus | 6.3 | 3− | 0.6 | 51 | 38 | 102 | 4.5510 | other | P/2 |
| trojan | (15527) 1999 YY2 | 6.9903 | 3 | 0.54 | 131 | 121 | 180 | 6.9921 | P | P/2 |
| main belt | (5740) Toutoumi | 11.68 | 3− | 0.43 | 348 | 247 | 179 | 5.8408 | P/2 | other |
| main belt | (782) Montefiore | 4.0728 | 3 | 0.54 | 52 | 51 | 179 | 4.0725 | P | daily alias |
| main belt | (17711) 1997 WA7 | 3.66176 | 3 | 0.52 | 128 | 108 | 177 | 3.6750 | P | other |
| main belt | (1081) Reseda | 7.3002 | 3 | 0.34 | 119 | 101 | 177 | 5.2468 | daily alias | other |
| main belt | (8080) Intel | 3.5335 | 3 | 0.37 | 46 | 42 | 161 | 3.5330 | P | other |
| main belt | (4952) Kibeshigemaro | 9.166 | 3 | 0.47 | 152 | 24 | 49 | 5.1066 | other | other |
| main belt | (41545) 2000 RA43 | 6.283 | 3 | 0.58 | 151 | 137 | 173 | 6.2806 | P | other |
| main belt | (10166) Takarajima | 5.3712 | 3 | 0.84 | 687 | 560 | 148 | 2.6855 | P/2 | other |
| main belt | (1631) Kopff | 6.683 | 3 | 0.41 | 135 | 121 | 176 | 7.7643 | other | other |
| main belt | (1572) Posnania | 8.048 | 3 | 0.4 | 139 | 117 | 171 | 4.8387 | daily alias | other |

Amplitudes are the LCDB's largest. Every column, with the fit's formal errors,
amplitude and scatter, is in [`rotation-sample.csv`](../assets/rotation-sample.csv).

| | Recovered (P or P/2) | Daily alias | Other | Not measured |
|:--|:--|:--|:--|:--|
| near-Earth | 4 of 10 (1 + 3) | 1 | 5 | 0 |
| Jupiter trojan | 3 of 9 (1 + 2) | 1 | 5 | 1 |
| main belt | 6 of 10 (4 + 2) | 2 | 2 | 0 |
| **All** | **13 of 29** (6 + 7) | **4** | **12** | **1** |
| Control, shuffled | 3 of 29 | 2 | 24 | |

(2674) Pandarus: of 32 exposures, 2 measurements were left, too few for the fit.

## 5. Discussion

**Against chance.** The control is the rate at which a period from these
epochs falls within 1% of the published one, or of its half, with no light
curve in the data: 3 of 29, about one in ten. 13 of 29 is well above it
(one-sided Fisher exact test, p = 0.004), but at that rate about three of the
twenty-nine would land within 1% by chance alone, so not every one of the
thirteen need be real. In the control, (23135) Pheidas lands
on its published period while its own data don't: its fitted amplitude, 0.09
mag, is barely larger than its scatter, 0.08.

**How close.** The recovered periods differ from the published ones (or their
halves) by a median 0.015%, 0.004–0.36%, against the 1% allowed; (17711)
1997 WA7 is the 0.36%. The fit's formal errors, 0.00005–0.0003 h, are smaller
than these differences, and smaller than the real uncertainty.

**A period or its half.** Seven of the thirteen are at half the published
period. The data can't tell a curve with two maxima at period P from one with a
single maximum at P/2; the rule counted both, and choosing between them takes an
assumption, as for Hektor, or a dense light curve.

**The failures.** Not declared beforehand, so only as a description: the three
longest periods, 17–20 h, all fail, where a night's visits cover a small part of
a rotation and the daily sampling dominates; five of the twelve "other" have
fewer than 40 measurements (Orthos, 2002 NP1, Thersilochos, Othryoneus,
Kibeshigemaro), and two a fitted amplitude under twice their scatter (Thoas,
Pheidas). Four landed on daily aliases, as (884) Priamus did in the
test of eight.

**Against the test of eight.** Two of six there, where only half the period
counted (the doubled fit); thirteen of twenty-nine here, where the period itself
counts too (six of the thirteen). The two samples also differ: these were chosen by their
number of ZTF exposures and with a larger amplitude.

## 6. Conclusion

From one season of ZTF survey data and nothing else, the period or its half
came out for 13 of 29 asteroids chosen beforehand, against 3 of 29 by chance.
Such a period is a candidate: to be checked for daily aliases and for the
factor of two, and confirmed with a dense light curve.

## 7. Reproducing

```sh
julia --project=. examples/rotation_period_sample.jl --select   # the sample only
julia --project=. examples/rotation_period_sample.jl            # every measurement
```

The images, about 1.3 MB each, and the results go to
`data/real/rotation_sample/`; a run that stops resumes where it left off.
