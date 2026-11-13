# Already reported: the Isolated Tracklet File

## Summary

Of the 5 tracklets [`analyze_set`](@ref) left unidentified on the 19
Pan-STARRS1 sets of IASC's October 2026 campaign searched by one team, 2 were
already in the MPC's Isolated Tracklet File (ITF): reported by Pan-STARRS's
own pipeline, in 3 of the 4 exposures each, 0.10–0.21″ from our positions.
IASC rejects such a tracklet as a duplicate. The other 3 were not reported.

## 1. Method

The ITF holds the observations reported to the MPC that are not yet linked to
any object, Pan-STARRS's unlinked detections among them. [`download_itf`](@ref)
keeps a copy (gzipped, ~130 MB, ten million observations, refreshed daily);
[`read_itf`](@ref) reads the observations of the set's dates; and
[`itf_check`](@ref) finds, for each unidentified tracklet, those within 2″ of
its straight-line motion and within 30 minutes of it. [`analyze_set`](@ref)
does it with its `itf` keyword.

## 2. Data

Declared before running, in `examples/itf_test.jl`:

| | |
|:--|:--|
| Sets | the 19 Pan-STARRS1 sets of IASC's October 2026 campaign searched by the NEO Hunters team |
| Search | [`analyze_set`](@ref) with its defaults, observatory F51, the MPC's orbits |
| ITF | as downloaded on 2026-11-11 |

## 3. Results

| Set | Unidentified | In the ITF |
|:--|:--|:--|
| XY75_p01 | 1 | P12rfge (F51): 3 observations, 0.10–0.21″ from the motion |
| XY45_p10 | 1 | P12r043 (F51): 3 observations, 0.15–0.21″ |
| XY04_p01, XY45_p01, XY55_p01 | 1 each | — |
| the other 15 | 0 | — |

Every ITF observation matched is from the same exposure as one of the set's,
to the second: Pan-STARRS's own report of the same object. The set
XY73_p00 was marked unreliable and not searched.

## 4. Discussion

**Duplicates.** Pan-STARRS's pipeline reports what it detects in 3 of 4
exposures; IASC's searchers are credited only for what it missed. Two of our
five candidates were in that first group and would have been rejected. Reading
one night of the file takes about 11 s.

**What it is not.** An ITF match says the tracklet was reported; it does not
say it is real, or that what is not in the ITF is new. The ITF changes daily
as tracklets are linked and leave it, and those then identify as known
objects instead.

## 5. Reproduce

```bash
julia --project=. examples/itf_test.jl <set folder or .zip> ...
```

The ITF and the MPC's orbits are downloaded on the first run; SkyBoT, Gaia and
JPL Horizons are asked live.
