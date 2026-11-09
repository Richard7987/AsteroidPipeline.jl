# Preliminary orbits

## Summary

[`fit_orbit`](@ref) fitted 40 short arcs of real MPC astrometry of 15 known
asteroids, cut as a new discovery would have them: 3 nights, 10 days and
30 days. Against JPL's orbits for the same objects, the semi-major axis was
off by a median 4.8% from 3 nights, 0.5% from 10 days and 0.03% from
30 days, and the position predicted 30 days later by a median 37″, 15″ and
1.2″. The predicted position fell within the fit's own 1σ in 36 of the 40
arcs and within 2σ in 39. The one beyond is an arc JPL's orbit fits no
better than the one found: its uncertainty is not a linear one.

## 1. Method

Gauss's method on the first, middle and last observation gives a first
orbit per root of its eighth-degree equation; systematic ranging (Virtanen,
Muinonen & Bowell 2001) gives others: a grid of distances at the first and
last observation, joined by Lambert's problem, of which the five fitting
best are kept. Each is refined by least squares on the heliocentric state at
the middle epoch, the observers' positions from JPL Horizons, light time
included, and the best kept. Observations off by more than 3σ are left out
one at a time.

The motion is two-body around the Sun: no planet perturbs it. The
uncertainties come from the fit's covariance, propagated linearly.

## 2. Data

Declared before running, in `examples/orbit_test.jl`:

| | |
|:--|:--|
| Objects | the numbered asteroids 1000, 2000, …, 10000, chosen by number only; the near-Earth asteroids (433) Eros, (1036) Ganymed, (1627) Ivar; the Jupiter trojans (588) Achilles, (624) Hektor |
| Observations | the MPC's CCD observations of each from 2024-01-01 on, every station, 0.5″ each |
| Arcs | from the first night: the first 3 nights; everything within 10 days; within 30 days |
| Compared | osculating a, e, i against JPL Horizons' at the same epoch; the position 30 days after the arc, from the geocentre, against Horizons' ephemeris |

Five arcs had fewer than 3 nights and were not fitted. For (10000), the
3-night and 10-day arcs are the same observations, counted twice.

## 3. Results

| Arc | Arcs fitted | Median rms | a off, median (max) | +30 d off, median (max) | Beyond 2σ |
|:--|:--|:--|:--|:--|:--|
| 3 nights | 15 | 0.10″ | 4.8% (18%) | 37″ (911″) | 0 |
| 10 days | 11 | 0.14″ | 0.5% (2.7%) | 15″ (100″) | 0 |
| 30 days | 14 | 0.17″ | 0.03% (13%) | 1.2″ (34″) | 1 |

Every arc, object by object:

| Object | 3 nights: a (JPL) | 10 days: a | 30 days: a | +30 d off / σ, the three arcs |
|:--|:--|:--|:--|:--|
| (1000) | 3.29 ± 0.29 (3.170) | 3.166 ± 0.007 | 3.1706 ± 0.0007 | 345″/438″, 16″/27″, 0.3″/2.9″ |
| (2000) | 2.382 ± 0.017 (2.379) | 2.379 ± 0.002 | 2.3794 ± 0.0001 | 27″/84″, 1.5″/11″, 0.1″/1.2″ |
| (3000) | 2.04 ± 0.09 (2.351) | 2.349 ± 0.004 | 2.352 ± 0.001 | 911″/1451″, 8.2″/14″, 4.4″/6.1″ |
| (4000) | 2.1 ± 0.6 (2.589) | — | 2.26 ± 0.03 | 37″/68″, —, 34″/8.0″ |
| (5000) | 2.3 ± 0.2 (2.537) | — | 2.4 ± 0.2 | 52″/36″, —, 22″/36″ |
| (6000) | 2.5 ± 0.4 (2.603) | 2.6 ± 0.2 | 2.601 ± 0.007 | 25″/93″, 4.8″/38″, 1.3″/4.3″ |
| (7000) | 2.47 ± 0.05 (2.464) | — | 2.47 ± 0.02 | 1.9″/19″, —, 2.4″/8.5″ |
| (8000) | 2.96 ± 0.14 (3.065) | — | — | 19″/17″, —, — |
| (9000) | 2.12 ± 0.11 (2.230) | 2.17 ± 0.06 | 2.230 ± 0.002 | 147″/174″, 100″/93″, 1.2″/10″ |
| (10000) | 2.57 ± 0.08 (2.586) | same arc | 2.5854 ± 0.0006 | 37″/241″, same, 1.5″/6.5″ |
| (433) Eros | 1.36 ± 0.19 (1.458) | 1.46 ± 0.07 | 1.457 ± 0.002 | 34″/266″, 12″/105″, 4.8″/5.6″ |
| (1036) Ganymed | 2.53 ± 0.13 (2.666) | 2.60 ± 0.08 | 2.668 ± 0.007 | 44″/42″, 24″/28″, 0.8″/5.0″ |
| (1627) Ivar | 1.8 ± 0.5 (1.864) | 1.8 ± 0.2 | 1.864 ± 0.005 | 50″/315″, 35″/132″, 0.4″/5.4″ |
| (588) Achilles | 5.5 ± 1.9 (5.211) | 5.18 ± 0.16 | 5.212 ± 0.004 | 21″/440″, 16″/74″, 0.8″/3.9″ |
| (624) Hektor | 5.3 ± 1.1 (5.273) | 5.28 ± 0.06 | 5.273 ± 0.002 | 158″/355″, 3.8″/40″, 0.8″/2.1″ |

## 4. Discussion

**The uncertainties.** 36 of 40 predictions fell within 1σ, more than the 68%
a Gaussian would give: with 0.5″ per observation against residuals of
0.1–0.2″, the stated uncertainty is conservative for these data. Pass the
observations' own uncertainties as `sigma` when they are known.

**The exception.** (4000)'s 30-day arc is 20 observations on 5 nights from two
ATLAS stations; its fit has a = 2.26 ± 0.03 AU, JPL's 2.589, 13σ away. Yet
JPL's orbit fits those 20 observations no better: rms 0.47″, against 0.44″ for
the orbit found. The data allow both; the linear covariance describes only
the neighbourhood of the minimum, and a short, sparse arc's uncertainty is
not a linear one. The same arc's own prediction for 30 days later was off by
34″ for a σ of 8″.

**What it is for.** A first orbit: is the object a near-Earth asteroid, a
trojan, a main-belt one; where to look for it in the coming nights; whether
two nights' tracklets can be the same object. Not an orbit for impact
monitoring or precise ephemerides, which need planetary perturbations and the
MPC's or JPL's processing of every observation.

**What wasn't tested.** Arcs of the pipeline's own measurements (no set here
has an object observed on several nights), comets, and objects passing
within a few hundredths of an AU of the Earth.

## 5. Reproducing

```sh
julia --project=. examples/orbit_test.jl
```

It fetches the observations from the MPC into `data/real/orbits` and asks JPL
Horizons for the observers' positions and the comparison orbits.
