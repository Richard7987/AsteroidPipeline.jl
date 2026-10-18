# Is a candidate new?

Over a million asteroids are catalogued; most moving objects in a set are
among them. Identifying them is a comparison of positions — with two
complications.

## Positions are predictions

A catalogue holds orbits, not positions. The position at the moment of an
exposure, seen from the observatory, is computed: here, from the MPC's orbit
file (two-body motion from the orbit's epoch, light-time corrected, the
observatory's position from [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/)),
and from [SkyBoT](https://vo.imcce.fr/webservices/skybot/). A prediction is
as good as the orbit:

| Orbit | Error of the prediction |
|:--|:--|
| numbered asteroid, many oppositions | well under 1″ |
| a few oppositions | arcseconds |
| short arc, or a recently changed orbit | arcminutes, or degrees |

So the match radius (5″ for SkyBoT, 15″ for the MPC's orbits, in three
frames) catches the well-known objects, and misses the poorly known ones —
which then show up as candidates.

## Motion tells

An imprecise orbit is mostly wrong in **time**: it puts the object ahead of
or behind where it is, on the right path. So a known object predicted near a
candidate, at the same rate, in the same direction, displaced **along** the
track, is almost certainly that object; displaced **across** the track, it is
another one.

![Same motion, displaced along the track or across it](../assets/nearby-geometry.png)

[`nearby_known_objects`](@ref) applies this to every candidate, against every
known object within 10′:

| Verdict | When |
|:--|:--|
| `:known` | a known object is predicted within 15″ |
| `:possibly_known` | one moves at the same rate (±10%, or 2″/h) in the same direction (±5°), displaced mostly along the track (across it by less than 30″, or 15% of the offset) |
| `:new` | no known object fits |

The tests have to be this strict. Within 10′, main-belt asteroids move
alike: around a real candidate on a Pan-STARRS1 set (the table below), two
unrelated objects moved within 25% of its rate and 20° of its direction. With
looser tests both would have "explained" it — wrongly: it is neither.

| Object | Separation | Rate | Direction | Verdict |
|:--|--:|--:|--:|:--|
| the candidate | — | 34.2″/h | 241° | |
| 2018 LT | 4.8′ | 43.0″/h | 212° | different motion |
| (610588) 2006 CX81 | 5.6′ | 31.8″/h | 232° | different motion (9° off) |
| (683089) 2007 JO28 | 8.1′ | 31.6″/h | 256° | different motion |

## Still not certain

A candidate with no known object near it is new *to the published orbits*.
It may already be on the MPC's [NEO Confirmation Page](https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html)
or reported by someone else that night; the MPC decides when the report
arrives.
