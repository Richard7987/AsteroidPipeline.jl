# MPC digest2 Scoring

Scores each tracklet in an [`astrometric_calibrate`](@ref)-shaped table
for its likelihood of belonging to several real solar-system orbit
classes (NEO, Main Belt, Jupiter Trojan, ...), via `digest2` — the
Minor Planet Center's own real, already-validated short-arc classifier
(statistical ranging against a real population model), not anything
trained or guessed at here. Explored alongside three other candidate
technologies (see [Design refinements](@ref) for why GPU-accelerated
reprojection and local plate-solving were investigated and *not*
recommended); `digest2` and [Cross-Night Linking](@ref) were the two
worth building.

## Setup (one-time, outside this package)

`digest2` is a separate C program, not a Julia dependency:

```sh
git clone https://github.com/Smithsonian/digest2.git
cd digest2/digest2
make
```

The `MPC.config`, `digest2.model.csv`, and `digest2.obscodes` files this
needs at runtime already ship in that same `digest2/` directory —
nothing to write or copy yourself. Either put the resulting `digest2`
executable on `PATH`, or pass its path directly:

```julia
digest2_score(candidates, "I41"; digest2_path="/path/to/digest2/digest2/digest2")
```

## Example

```julia
scores = digest2_score(candidates, "I41")
```

produces one row per tracklet:

| id | rms | int_score | neo_score | n22_score | n18_score |
|--:|--:|--:|--:|--:|--:|
| 1 | 0.02 | 100.0 | 100.0 | 23.0 | 0.0 |
| 2 | 0.0 | 100.0 | 98.0 | 26.0 | 1.0 |

## A real finding: this is an orbit-class classifier, not a real/bogus one

Run for real against `real_data_demo.jl`'s ZTF field 451 baseline (133
tracklets from real, already-downloaded ZTF data; setup: `examples/fetch_data.sh`
then `examples/real_data_demo.jl`'s baseline stage), the result was the
*opposite* of a naive expectation:

| | `neo_score` |
|:--|--:|
| The 2 real, SkyBoT-confirmed known objects (2002 UY45, 1997 KO3 — both Main Belt) | 5, 7 |
| The other 131 tracklets | 100 (130 of 131) |

That's `digest2` working correctly, not a bug: it correctly recognized
the 2 real objects as *not* NEOs (they're Main Belt — their score is
concentrated in `digest2`'s `MB1`/`MB2` classes, which this function
doesn't parse out). The other 131 were never real moving objects at
all — inspecting one directly, its 5 detections jitter by under 1"
across the full 6.25 h baseline (`~3"/day` implied rate): a real,
stationary star, re-detected each frame and linked into a bogus
tracklet only because [`link_candidates`](@ref)'s `match_radius` in
that demo (10", looser than ZTF's real sub-arcsec centroiding
precision — the same lesson [Validating against real IASC (Pan-STARRS1) campaign data](../iasc-campaign-validation.md)
already learned once from `PERROR`, resurfacing here) was loose enough
to accept a star's own centroid jitter as if it were motion. `digest2`
has no way to tell a stationary star's jitter from a real, very distant
object moving at a dynamically-consistent near-zero rate — it scored
the hypothesis it was given correctly; the hypothesis itself was never
a real tracklet.

**Practical takeaway**: don't sort by `neo_score` on raw
[`link_candidates`](@ref) output as a "most interesting first" filter.
Tighten `match_radius` to the survey's own real astrometric precision
first, or otherwise filter for genuinely consistent motion, before
scoring — `digest2_score` classifies an *already-plausible* tracklet's
dynamical class well; it is not a substitute for that upstream quality
control.

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["digest2.jl"]
```
