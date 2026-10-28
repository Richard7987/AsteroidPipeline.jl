# Roadmap

What is being worked on, what is missing before it can be released, and
what comes after. Items move to the [release notes](../release-notes.md)
when they ship.

## Branches

`main` holds what has been released or is ready to be; new work goes on its
own branch until it is finished:

| Branch | For |
|:--|:--|
| `feature/<name>` | a new function or capability |
| `fix/<name>` | a bug fix too large for a single commit on `main` |
| `release-<x.y.z>` | fixes to an already released version |

Names are lowercase words separated by hyphens. A branch starts from `main`,
takes `main`'s later changes by merging it (a pushed branch is never
rewritten), and its release notes go under "Unreleased". CI runs on every
pushed branch; a branch is merged once it passes.

## In progress: light curves of moving objects

Branch `feature/moving-light-curve`.

**Done**

- [`moving_light_curve`](@ref): photometry that follows the object from
  frame to frame, with a zero point per frame against Gaia and a flag for
  frames where a star bright enough to matter lies in the aperture.
- Synthetic tests: motion, recentring, a star on the track, missing frames.
- [`reduced_magnitudes`](@ref) and [`object_geometry`](@ref): distance, phase
  angle and light time, with the geometry from JPL Horizons.
- ZTF check on (624) Hektor, 107 exposures in r, July–December 2022:
  scatter within a night 0.009 mag, consistent with the reported errors of
  0.003–0.012 mag. Folded on the published 6.9205 h period, July–December
  leaves 0.049 mag, better than every random period tried; the
  periodogram's highest peak is 3.460 h, half the period, as for any
  double-peaked light curve.

**Missing before merging**

- [ ] A period recovered from real data: several nights of one object with
      a known period, from a Mexican telescope (data expected), and with it
      the validation page (`validation/moving-light-curve.md`).
      `examples/moving_light_curve_demo.jl` runs the same checks on ZTF
      frames of Hektor, whose period is already known; a team's own
      frames test the whole path, from its telescope to a period.

**After**

- [ ] Positions for a known object from its ephemeris (JPL Horizons, as
      `mpcorb.jl` already uses for the observer), instead of supplying them
      by hand for each night.
- [ ] The light curve in the UI, starting from a tracklet in the results.
- [ ] Release 1.3.0.

## Pending on `main`

- [ ] Shorter warning texts in [`analyze_set`](@ref) when SkyBoT or other
      services fail.
- [ ] Tests for `neocp.jl` and plate solving using recorded responses, so
      Codecov measures them without depending on the services.
- [ ] Observer positions computed locally with SPICE (SPICE.jl +
      SPICEKernels.jl), as suggested on Julia's Discourse, instead of a
      request to Horizons.
