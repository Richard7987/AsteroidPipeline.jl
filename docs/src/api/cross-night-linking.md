# Cross-Night Linking

[`link_candidates`](@ref) only links detections *within* one observing
night (a flat sequence of frames sharing one linear-motion model over
minutes to hours). [`link_across_nights`](@ref) solves the same kind of
problem one level up: grouping different nights' own tracklets that are
consistent with one real object, over gaps of days where real orbital
curvature — not just the instantaneous rate — starts to matter. This was
entirely greenfield: nothing in this pipeline had any night/date/session
concept before this (confirmed by a direct search of `src/` — `timestamps`
was, and still is for a single night, one flat `Vector{Float64}` of
Julian Dates in lockstep with `fits_paths`).

## Example

```julia
night1 = search_field(night1_paths; timestamp_key="OBSMJD").movers
night2 = search_field(night2_paths; timestamp_key="OBSMJD").movers
night3 = search_field(night3_paths; timestamp_key="OBSMJD").movers

groups = link_across_nights([night1, night2, night3])
```

returns a `Vector` of groups, each a `Vector` of `(night, id)` pairs —
tracklets, from different nights, judged to be the same real object.
Only groups spanning at least `min_nights` (default `2`) distinct nights
are returned.

## Design

Each night's tracklets are independently fit to a linear rate in a
local tangent-plane projection (reusing `_linfit` — the exact
same fitting routine [`link_candidates`](@ref)'s own within-night refit
already uses, not reimplemented), then extrapolated forward to every
other night's tracklets' own mean epoch. A pair is accepted within
`match_radius_arcsec` of the extrapolated position, optionally loosened
by `max_accel_arcsec_per_day2 * Δt_days^2` to allow for real curvature
over longer gaps (off, `Inf`, by default — mirroring
[`link_candidates`](@ref)'s own `max_speed::Real=Inf` default: a single
generous radius does the real work unless the caller has independent
reason to tighten it). Accepted pairs are grouped transitively via
union-find, since one tracklet can plausibly pair with several other
nights' tracklets that must all collapse into one group.

## An honest gap: not yet validated against a real multi-night case

Every other real-data claim in this project's docs is backed by an
actual run against real data. This one isn't yet, and that's stated
plainly rather than skipped past.

The existing real IASC/PS1 dataset was checked first (5 fields, 9 known
objects — see [Validating against real IASC (Pan-STARRS1) campaign data](../iasc-campaign-validation.md)): none of
those 9 objects repeats across fields, so there was no ready-made real
cross-night case sitting in data already on hand. A further, bounded
real check — querying SkyBoT for ZTF field 451's own sky position at
its original epoch and again 3 and 7 days later — confirmed *why* this
is genuinely hard to find by chance: every object SkyBoT reports within
that field changes completely from one query to the next. Even a
"slow" Main Belt object here (2002 UY45, ~852"/day) crosses a ZTF
quadrant's own ~35' width in a few days — real multi-night, same-tile
recovery needs an object caught unusually close to its apparent
stationary point, not just "a slow-ish asteroid," and finding one that
also happens to fall within ZTF's actual public observing history for
one specific tile is a real needle-in-a-haystack search, not something
a couple of cone searches turns up.

So: [`link_across_nights`](@ref)'s correctness is currently established
by synthetic exact-recovery tests only (construct a real linear-motion
object across 3 fabricated nights, plus unrelated same-night "noise"
tracklets, and confirm the object's tracklets group correctly while
noise never does — mirroring [`link_candidates`](@ref)'s own synthetic
test style exactly). `match_radius_arcsec`'s default (5") and
`max_accel_arcsec_per_day2=Inf` are starting points to tune against a
real campaign's own cadence and target population, not values
calibrated against this project's own measured data the way e.g.
`variability_chi2`'s `systematic_error_fraction` was. A real multi-night
validation remains open for whenever a suitable real case turns up (a
live campaign's own multi-night data would be the natural source).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["linking_multinight.jl"]
```
