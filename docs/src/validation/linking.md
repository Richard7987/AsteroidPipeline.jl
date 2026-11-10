# Linking across nights

## Summary

[`link_nights`](@ref) was given the real astrometry of every catalogued
asteroid that one station observed in one patch of sky over two weeks: 323
tracklets of 73 asteroids on 12 nights, mixed and without their names. It
linked all 47 objects seen on 3 nights or more, 43 of them with every
tracklet, and made no wrong linkage. With 969 false tracklets moving like
asteroids added, the 3-night linkages were again 47 right and none wrong; on
2 nights, 110 linkages were wrong and 6 right. A linkage of 2 nights is
therefore not counted unless asked for (`min_nights`).

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
observation lies within 3σ (1.5″ for the default σ = 0.5″) of it. Linkages
grow one tracklet at a time, the one fitting best, while the orbit still
fits every observation; each tracklet ends in one linkage at most, longer
linkages first.

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
unchanged.

## 3. Results

| | Real tracklets only | With 969 false ones |
|:--|:--|:--|
| Linkages of 3 nights or more | 47 right, 0 wrong | 47 right, 0 wrong |
| Linkages of 2 nights | 7 right, 0 wrong | 6 right, 110 wrong |
| Objects on 3+ nights (47): all tracklets linked | 43 | 43 |
| in part | 4 | 4 |
| not linked | 0 | 0 |
| Objects on 2 nights (7) linked | 7 | 6 |
| False tracklets in a linkage | — | 215, all on 2 nights |
| Orbit rms of the right linkages | median 0.34″, max 0.49″ | the same |
| Time | 68 s | 269 s |

Objects linked in part, by the number of nights they were seen on: one of 10
seen on 6 nights, two of 11 on 8, one of 2 on 9. Each lost one or two
tracklets.

## 4. Discussion

**The tracklets left out.** Each of the five tracklets missing from those
four objects' linkages has one point between 1.6″ and 2.7″ from the orbit
of the rest, beyond the 1.5″ allowed. The rule is strict on purpose: a
linkage is accepted only when every observation fits. A survey whose
astrometry is better than 0.5″ should pass its own σ.

**Two nights.** Among false tracklets moving like asteroids, two
tracklets of different nights fit a two-body orbit to 1.5″ often: 110
times here, against 6 real objects. A third night ruled every one of them
out. [`link_nights`](@ref) therefore keeps linkages of 3 nights or more
by default; that default was chosen after this measurement. Five real
tracklets ended in a wrong 2-night linkage: three of objects seen on one
night only, one left out of a linkage in part, and one of an object seen
on 2 nights, which was then not linked.

**What the test is not.** All the real tracklets here are of catalogued
objects, bright enough for one station to measure on several nights. A
campaign's own candidates are fainter, its false tracklets are of stars,
noise and artifacts rather than straight lines at asteroid rates, and the
real objects here are the catalogued ones that station measured, not every
object as faint. The false tracklets were made to be hard: they move like
asteroids by construction.

**Speed.** 68 s for 323 tracklets, 269 s for four times as many, on one
core of a desktop computer.

## 5. Reproduce

```bash
julia --project=. examples/linking_test.jl
```

The first run fetches each object's observations from the MPC into
`data/real/linking`; JPL Horizons is asked live.
