# Linking across nights

## Summary

[`link_nights`](@ref) was given the real astrometry of every catalogued
asteroid that one station observed in one patch of sky over two weeks: 323
tracklets of 73 asteroids on 12 nights, mixed and without their names. It
linked every tracklet of all 54 objects seen on 2 nights or more, leaving
out 5 points that did not fit, and made no wrong linkage. With 969 false
tracklets moving like asteroids added, the linkages of 3 nights or more were
again 47 right and none wrong; on 2 nights, 168 were wrong and 7 right. A
linkage of 2 nights is therefore not counted unless asked for
(`min_nights`).

## 1. Method

Two tracklets of different nights, at most 15 days apart, are first compared
by their motion: for any motion of constant acceleration, the mean angular
velocity between them equals the mean of their own two velocities. The
allowance is 3σ of the measured velocities plus a bound on the next term.
The observer's offset from the Earth's centre is taken out first, for a grid
of distances from 0.02 to 50 AU, since its turning with the Earth changes a tracklet's apparent
motion by up to 2.3″/h for an object 1 AU away; a pair passes at any one
distance. A pair that passes is fitted with a two-body orbit, by the method
of [`fit_orbit`](@ref) started from those distances; it is linked when every
observation lies within 3σ (1.5″ for the default σ = 0.5″) of it. If one
does not, the point farthest from the orbit is left out and the orbit
refitted: one point of a tracklet, and only of a tracklet of 3 points or
more. Linkages grow one tracklet at a time, the one fitting best, while the
orbit still fits by that rule; each tracklet ends in one linkage at most,
longer linkages first.

## 2. Data

Declared before running, in `examples/linking_test.jl`:

| | |
|:--|:--|
| Window | 2025-10-14 to 2025-10-29 UTC: the new moon of 2025-10-21 ± 7 days |
| Objects | every MPCORB orbit within 1° of the anti-solar point on 2025-10-21.5 UTC: 297 |
| Station | Pan-STARRS1 (F51) was declared; it did not observe these objects in the window, so the station with the most CCD observations of them there was used: P07 (Space Surveillance Telescope), 1,713 observations |
| Tracklets | the station's observations of one object on one night, 2 points or more: 323, of 73 objects; 4–10 points over a median 3.9 h |
| Linking | the defaults (σ = 0.5″, 15 days, 3σ), linkages of 2 nights or more |
| False tracklets | declared after the first run: 3 per real tracklet on each night (969), each with the epochs of one real tracklet of that night, the straight-line motion of another, a start anywhere in the box the night's real points cover, and 0.3″ of noise per point |

The first run found a fault in the package, not in the linker: JPL Horizons
answers in time order, and the observer's positions asked out of time order
were matched to the wrong epochs. With it fixed, the run was repeated
unchanged. The run was repeated again after the points left out (§4) were
allowed; the results below are of that version.

## 3. Results

| | Real tracklets only | With 969 false ones |
|:--|:--|:--|
| Linkages of 3 nights or more | 47 right, 0 wrong | 47 right, 0 wrong |
| Linkages of 2 nights | 7 right, 0 wrong | 7 right, 168 wrong |
| Objects on 3+ nights (47): all tracklets linked | 47 | 47 |
| Objects on 2 nights (7) linked | 7 | 7 |
| Points left out | 5, in 4 linkages | 101, in 87 linkages |
| False tracklets in a linkage | — | 328, all on 2 nights |
| Orbit rms of the right linkages | median 0.34″, max 0.49″ | the same |
| Time | 114 s | 526 s |

Before points could be left out, 4 of the 47 objects on 3 or more nights
were linked in part, one of the 7 on 2 nights was not linked, and on 2
nights 110 linkages were wrong and 215 false tracklets linked, in 68 s and
269 s.

## 4. Discussion

**The points left out.** Before points could be left out, five tracklets
were missing from four objects' linkages, each with one point between 1.6″
and 2.7″ from the orbit of the rest, beyond the 1.5″ allowed. Leaving out
one point per tracklet links them all. The rule stays strict otherwise: a
tracklet of 2 points loses none, and every other observation must fit. A
survey whose astrometry is better than 0.5″ should pass its own σ.

**Two nights.** Among false tracklets moving like asteroids, two
tracklets of different nights fit a two-body orbit to 1.5″ often: 168
times here, against 7 real objects; leaving out a point made it more
often (110 before). A third night ruled every one of them
out. [`link_nights`](@ref) therefore keeps linkages of 3 nights or more
by default; that default was chosen after this measurement.

**What the test is not.** All the real tracklets here are of catalogued
objects, bright enough for one station to measure on several nights. A
campaign's own candidates are fainter, its false tracklets are of stars,
noise and artifacts rather than straight lines at asteroid rates, and the
real objects here are the catalogued ones that station measured, not every
object as faint. The false tracklets were made to be hard: they move like
asteroids by construction.

**Speed.** 114 s for 323 tracklets, 526 s for four times as many, on one
core of a desktop computer; refitting without a point about doubled both.

## 5. Reproduce

```bash
julia --project=. examples/linking_test.jl
```

The first run fetches each object's observations from the MPC into
`data/real/linking`; JPL Horizons is asked live.
