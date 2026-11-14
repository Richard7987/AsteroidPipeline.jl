# How faint, on campaign sets

## Summary

On 10 Pan-STARRS1 sets of IASC campaigns, with 100 synthetic asteroids
injected into each (V 21.0–23.5), [`analyze_set`](@ref) found half of them
down to V 21.68 frame by frame and V 21.87 with stacking. Between V 22.0 and
22.5, stacking found 34 of 192 where the frame-by-frame search found 8;
fainter than V 22.5, they found 3 and 4 of 369. Pan-STARRS's own pipeline reaches
about V 22 and IASC's searchers find objects 0.30 mag fainter on average
(Miller et al. 2024): the pipeline does not yet reach the range where they
work.

## 1. Method

[`inject_movers`](@ref) adds asteroids to the set's frames as point sources
with the frames' measured FWHM, moving in straight lines and smeared along
their motion during each exposure; [`analyze_set`](@ref)
searches the result with its defaults, and [`injection_recovery`](@ref) counts
an injected asteroid as recovered when a tracklet follows it. A tracklet
counts for the frame-by-frame search if it is not in the result's `stacked`.
The magnitude where half are found is interpolated between bins of 0.25 mag.

## 2. Data

Declared before running, in `examples/faint_completeness.jl`:

| | |
|:--|:--|
| Sets | the 13 local Pan-STARRS1 sets of IASC campaigns (4 × 45 s each) |
| Injected | 100 per set, V 21.0–23.5 and 5–100″/h drawn uniformly, seed 1000 plus the set's place in the list |
| Search | [`analyze_set`](@ref) with its defaults (stacking on), observatory F51, the MPC's orbits |

Changed after the first run stopped at `XY72_p11`, before any result was seen:
a set whose injection fails, or that [`analyze_set`](@ref) marks unreliable,
is left out. That left out three: `XY14_p10`, the damaged set, marked
unreliable; `XY72_p11` and `XY73_p00`, where [`inject_movers`](@ref) found no
photometric zero point in one frame and stopped, though
[`analyze_set`](@ref) analyses both sets.

## 3. Results

All 10 sets together, 1,000 injected:

| V | Injected | Frame by frame | With stacking |
|:--|:--|:--|:--|
| 21.00–21.25 | 114 | 92 (81%) | 94 (82%) |
| 21.25–21.50 | 114 | 85 (75%) | 87 (76%) |
| 21.50–21.75 | 90 | 50 (56%) | 57 (63%) |
| 21.75–22.00 | 121 | 39 (32%) | 60 (50%) |
| 22.00–22.25 | 92 | 8 (9%) | 17 (18%) |
| 22.25–22.50 | 100 | 0 | 17 (17%) |
| 22.50–22.75 | 96 | 2 (2%) | 3 (3%) |
| 22.75–23.50 | 273 | 1 | 1 |
| **50% found at** | | **V 21.68** | **V 21.87** |

Each set, with the 5σ limit [`analyze_set`](@ref) reports for it:

| Set | 5σ limit | 50% frame by frame | 50% with stacking | Stacked, not injected |
|:--|:--|:--|:--|:--|
| `XY15_p01` | 22.23 | 21.22 | 21.27 | 3 |
| `XY25_p10` | 23.12 | 22.19 | 22.19 | 1 |
| `XY26_p01` | 23.00 | 21.98 | 21.98 | 0 |
| `XY42_p11` | 22.82 | 21.95 | 22.00 | 0 |
| `XY45_p01` | 22.43 | 21.62 | 21.85 | 1 |
| `XY45_p10` | 22.39 | 21.46 | 21.93 | 1 |
| `XY45_p11` | 22.31 | 21.46 | 21.50 | 0 |
| `XY54_p10` | 22.81 | 21.78 | 21.85 | 0 |
| `XY72_p01` | 22.64 | 21.76 | 21.79 | 0 |
| `XY73_p01` | 22.66 | 21.55 | 21.55 | 0 |

The 6 stacked tracklets not injected are the ones stacking adds on the same
sets without injections ([Fainter objects by stacking](stacking.md)): two
catalogued asteroids and a diffuse patch on `XY15_p01`, and one candidate each
on `XY25_p10`, `XY45_p01` and `XY45_p10`.

## 4. Discussion

**Against IASC.** Miller et al. (2024) give Pan-STARRS's pipeline a limit of
about V 22, from 3 detections in 4 exposures, and find the objects IASC's
searchers report 0.30 ± 0.07 mag fainter on average than those it reports.
Neither is a measured completeness, so the comparison is loose; but here no
set reaches 50% at V 22.2, and between V 22.0 and 22.5 the pipeline finds
about one in six. The direct case is NHU0002, at V≈22.3: found on its set,
and not in the MPC's Isolated Tracklet File ([Already reported](itf.md)).

**Stacking.** It moves the 50% point by 0.2 mag on all sets together, from 0
on three sets to 0.47 on `XY45_p10`, and finds all of the 17 recovered at V
22.25–22.5. Past V 22.5 it finds almost nothing, though the 5σ limits of
six sets are fainter than that.

**The 5σ limit is not the depth.** The 50% point with stacking is 0.46–1.11
mag brighter than the limit [`analyze_set`](@ref) reports, 0.89 on the
median: that limit is a single source's in one frame, and a tracklet
needs the object in several.

**Against the earlier test.** The injection test on the
[Pan-STARRS1 sets](pan-starrs1.md#Completeness), with an earlier version,
found 21% at V 22–22.5 frame by frame; here it is 4% (8 of 192). The
difference was not looked into.

**What it does not cover.** Injected asteroids move in straight lines with a
constant brightness; real ones may not. Only Pan-STARRS1, only
sets of four frames. Why stacking stops near V 22.5 was not measured.

## 5. Reproducing

```sh
julia --project=. examples/faint_completeness.jl
```

It reads the sets from `data/real/iasc` and keeps each set's result in
`data/real/faint_completeness`; a second run only prints the tables.

## References

Miller, C., et al. (2024). PASP, 136, 024502.
