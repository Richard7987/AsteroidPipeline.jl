# Near-Earth or not: digest2

## Summary

On the 18 searchable Pan-STARRS1 sets of IASC's October 2026 campaign
searched by one team, the MPC's digest2 gave 25 of the 26 catalogued
asteroids found a NEO score below 65, the score at which the MPC posts a
tracklet to its NEO Confirmation Page; none of the 26 is a near-Earth
object. The one above, 2026 SE8 (perihelion 1.46 au), scored 98. Of the 4
unidentified tracklets, one scored 69 and three 1–2; looked at, the one at 69
is not in any frame.

## 1. Method

digest2 (Keys et al. 2019) finds the orbits that could produce a tracklet's
motion, a few observations from one night, and scores each orbit class 0–100
by how much of a model population of the Solar System lies on those orbits.
[`digest2_score`](@ref) runs it on [`analyze_set`](@ref)'s tracklets as the
MPC does: with each station's astrometric uncertainty from the program's
`MPC.config` and the same seed for every tracklet. The NEO class is
perihelion q < 1.3 au; the "MPC interest" class (`Int`) adds orbits that are
eccentric, inclined or reach far out.

## 2. Data

Declared before running, in `examples/digest2_test.jl`:

| | |
|:--|:--|
| Sets | the 19 Pan-STARRS1 sets of IASC's October 2026 campaign searched by the NEO Hunters team |
| Search | [`analyze_set`](@ref) with its defaults, observatory F51, the MPC's orbits |
| digest2 | version 1.0.0, from `github.com/Smithsonian/mpc-public`, with its `MPC.config` |
| Truth | each identified asteroid's perihelion, from its MPC orbit |

## 3. Results

`XY73_p00` was marked unreliable and not searched. On the other 18:

| | Tracklets | NEO score |
|:--|:--|:--|
| Catalogued, not near-Earth (q 1.46–2.92 au) | 26 | median 9; 24 below 50; 2000 KO3 54, 2026 SE8 98 |
| Catalogued, near-Earth | 0 | — |
| Unidentified | 4 | 69 (`XY55_p01`); 2 (`XY45_p01`, `XY45_p10`); 1 (`XY04_p01`) |

2026 SE8 also scored 100 in the `Int` class; it is a Mars-crosser, just
outside the NEO limit. The tracklet on `XY55_p01` scored 92 in `Int`. Run
twice, every score was the same.

`XY75_p01`'s tracklet, unidentified in the [ITF check](itf.md), is now
identified, as 2026 VL6, among the MPC's orbits since. It scored 1.

## 4. Discussion

**What it says.** digest2 judges the motion of one night, not the object:
a tracklet of a main-belt asteroid is unlikely to score high, and here 25 of
26 didn't. One that moves like a near-Earth object can still be a
Mars-crosser on a nearby orbit, as 2026 SE8 is.

**What it doesn't.** It doesn't tell a real object from a false tracklet:
a star's jitter linked into a tracklet scores as a fast NEO
([Follow-up tools](../guides/follow-up.md)). The tracklet on `XY55_p01`, at
69, was found only by stacking, at 9.8″/h, with magnitudes of V 23.4–24.1
where measured: 0.8–1.5 mag fainter than the set's 5σ limit (22.6), where
the search found 1 of 297 injected asteroids between V 22.75 and 23.5
([How faint](faint.md)). Looked at, it is in none of the four frames nor in
their mean along its motion, and its light within 7 × 7 pixels is 0.5–3.6σ:
noise, scored as a slow near-Earth object.

**What wasn't tested.** No near-Earth object was among the identified
asteroids, so how often digest2 scores a real one below 65 here wasn't
measured; only one night of one survey.

## 5. Reproduce

```bash
julia --project=. examples/digest2_test.jl <digest2 executable> <set folder or .zip> ...
```

## References

Keys, S., et al. (2019). The digest2 NEO Classification Code. PASP, 131,
064501.
