# AsteroidPipeline.jl

An open-source Julia pipeline for asteroid search campaigns, developed for
use with the [International Astronomical Search Collaboration
(IASC)](https://iasc.cosmosearch.org/).

## Purpose

The pipeline processes sequences of FITS frames from a survey field to:

1. **Difference** each science frame against a deep, static-sky reference
   stack (`build_reference`, `estimate_psf`, `zogy_subtract`), extending
   detection below the single-frame noise floor — optional; skipped if no
   reference is supplied.
2. **Detect** point sources in each frame (`detect_sources`) — on the
   difference image if step 1 ran, on the raw science frame otherwise.
3. **Link** detections across frames by consistent linear motion to form
   asteroid candidate tracklets (`link_candidates`).
4. **Calibrate** tracklet pixel positions to sky coordinates using each
   frame's WCS astrometric solution (`load_wcs`, `pix_to_sky`,
   `astrometric_calibrate`).
5. **Cross-match** candidates against known-object catalogs — SkyBoT, VSX,
   SIMBAD — to separate previously cataloged objects from candidates that
   warrant human verification (`crossmatch_catalog`).

`run_pipeline` runs steps 1-4 end to end on a sequence of FITS file paths,
returning a candidate table ready for `crossmatch_catalog`.

For a confirmed discovery, `light_curve` (forced aperture photometry at a
fixed sky position across a dedicated follow-up sequence) and
`recover_rotation_period` (a Lomb-Scargle periodogram over that light
curve, via `LombScargle.jl`) recover a rotation period — a separate,
optional follow-up step, not part of `run_pipeline` itself.

Planned extension: variable-star and transient detection (distinguishing
these from asteroid candidates in `crossmatch_catalog`'s output).

## Status

Early development. `detect_sources`, `link_candidates`, the WCS
calibration step, and `crossmatch_catalog` are implemented and wired
together end to end in `run_pipeline`, validated against synthetic FITS
frames with a known injected source track, and exercised against real
public survey data (see `examples/real_data_demo.jl`). Not yet run on a
real IASC dataset.

ZOGY difference imaging (Zackay, Ofek & Gal-Yam 2016) is implemented —
`build_reference` stacks a deep reference from many epochs via
reprojection (`Reproject.jl`) onto a common pixel grid, `estimate_psf`
measures each frame's empirical PSF from its own bright stars, and
`zogy_subtract` produces a statistically normalized detection-significance
map (`S_corr`, unit variance by construction) in Fourier space via
`FFTW.jl`. Validated against three falsifiable synthetic checks (identical
images subtract to zero; pure noise gives `std(S_corr) ≈ 1`; an injected
source's peak significance matches the analytic matched-filter prediction)
and against real ZTF data — see `examples/real_data_demo.jl`, which runs
the pipeline with and without differencing on the same frames and reports
which known objects each recovers, rather than assuming differencing
helps.

Several real bugs surfaced only by testing against real ZTF data (not by
the synthetic tests or code review), all fixed with regression tests:

- An empty cross-match result crashing table construction.
- `crossmatch_catalog(...; :skybot)` silently returning zero matches on
  every real call because Julian Date epochs (~2.4e6) print in scientific
  notation by default, which SkyBoT's API rejects as an empty epoch.
- `zogy_subtract` never background-subtracted its inputs. A raw FFT's DC
  term is a *sum*, not a mean, over the whole image (>10^6 pixels here),
  so even a modest sky-brightness mismatch between a science frame and
  the reference stack (routine — depends on moon phase and airglow, not
  the photometric zeropoint `build_reference` already matches) blew up
  into a near-constant offset that swamped `S_corr` almost everywhere:
  on real data this meant >99.9% of a frame reading "above 6 sigma" and
  zero genuine detections surviving.
- Reprojected frames carry a thin, legitimate NaN border wherever a
  dithered footprint doesn't fully cover the target grid; left
  unsanitized before `zogy_subtract`, a single NaN poisoned its
  whole-array `median`/`fft` calls and turned an entire frame's `S_corr`
  to NaN.
- `run_pipeline` filled that same border with `0.0` rather than the valid
  region's own background level. Negligible for a real ~10^6 px ZTF frame
  (~0.05% invalid), but a regression test using a synthetic frame with a
  large dither (~13% invalid) showed this biases `zogy_subtract`'s own
  background-median computation toward 0, contaminating the whole
  difference image — the same bug, just large enough at that scale to
  turn from invisible into a wrong `run_pipeline` output.
- `estimate_psf` measured each frame's PSF *before* reprojection, but
  `zogy_subtract` differenced the frame *after* — interpolation reshapes a
  PSF slightly, and that mismatch showed up as systematic subtraction
  residuals at bright stars: on real data, 85-97% of the excess
  `S_corr` detections in 4 of 5 test frames fell within 15 px of one.
  Moving the PSF measurement after reprojection, and feeding both frames'
  detected stars to `zogy_subtract`'s astrometric-noise term (`V_ast`) so
  residual misregistration is priced into the significance map rather
  than ignored, brought `std(S_corr)` from ~2.0-2.6 down to ~1.1-1.2 in
  those 4 frames (nominally exactly 1 under the model's assumptions).

A further issue was a dataset-selection mistake rather than a code bug:
the first real-data demo used a field/night where the faintest catalogued
asteroid was 1.2 mag *below* that night's own detection limit — no
algorithm can recover a signal that was never above the noise, so the
demo now picks a night with a known object bright enough to serve as
actual ground truth.

With all of this fixed (including a per-frame quality gate, `run_pipeline`'s
`quality_max_std` — see Known limitations), `examples/real_data_demo.jl`'s
actual result on that night (field 451, 2019-10-23): the undifferenced
baseline finds 133 tracklets and recovers both known objects in the field
(2002 UY45, 1997 KO3); ZOGY also recovers both, at consistent sky offsets
(confirming the subtraction is correctly calibrated, not just "not
obviously broken"), but finds 667 tracklets total and no *additional*
known object. The likely reason ZOGY doesn't show its expected depth
advantage here: both known objects in this field (Mv 18.4, 19.6) are
bright enough that the undifferenced baseline already recovers them
trivially — this dataset doesn't happen to contain a known object faint
enough to sit below a single frame's noise floor but above the deep
reference's, which is the specific regime ZOGY is for.

The 667 vs 133 gap is *not* a clean measurement of ZOGY's noise
properties, and should not be read as one — see the last bullet under
Known limitations for why (the quality gate, once it correctly excludes
a bad frame, forces a `min_frames` reduction that itself loosens the
matching combinatorics, confounding a direct comparison).

## Known limitations

- **`plate_solve` has no live validation.** It implements the full
  nova.astrometry.net login/upload/poll/fetch cycle (see its docstring),
  and its request/response-parsing logic is unit-tested directly, but no
  API key was available while writing it, so the actual network round
  trip has never been exercised — set `ENV["ASTROMETRY_API_KEY"]` to run
  that test (`test/runtests.jl`, `@testset "plate_solve"`) and confirm it
  actually solves a real frame before relying on it.
- **Empirical PSF, not a fitted model.** `estimate_psf` stacks real star
  cutouts rather than fitting an analytic profile (Gaussian/Moffat), which
  keeps it survey-agnostic but means its quality depends on having enough
  bright, isolated, unsaturated stars in the frame.
- **`zogy_subtract`'s astrometric-noise term (`V_ast`) is opt-in at the
  `zogy_subtract` level** — it needs `n_sources`/`r_sources` passed
  explicitly, and is `0` without them. `run_pipeline` always supplies
  them, so this only matters when calling `zogy_subtract` directly.
- **One real-data frame likely had a passing cloud during the exposure —
  now caught by a quality gate, with a real combinatorial side effect.**
  On the `examples/real_data_demo.jl` field/night, `std(S_corr)` improved
  to ~1.1-1.2 (from ~2.0-2.6) in 4 of the 5 science frames after fixing
  the PSF-timing/astrometric-noise bug above, but the night's last
  exposure stayed at ~2.0 and produced roughly twice as many bright
  sources as the other four (232 vs 83-111). GAIN, SEEING, MAGZP,
  whole-frame background level, saturated-pixel count, and airmass were
  all checked against the other four frames and none stood out. Directly
  differencing this frame's *raw* pixels against another frame's (before
  any ZOGY machinery) found the actual cause: ~1750 pixels with a
  significant, one-sided (no matching negative lobe, so not a
  registration/dipole artifact) excess, ~98% confined to the bottom ~30%
  of the detector — a control pair of two normal frames differenced the
  same way found only ~80-110 such pixels, distributed in proportion to
  where the frame's stars actually are (not concentrated in one region).
  The whole-frame background stayed flat with no gradient across that
  region, ruling out amplifier glow or vignetting. This pattern — real
  stars showing localized excess halos in one part of the frame, with no
  diffuse background change and no positional offset — is the signature
  of thin cloud scattering starlight during the 30 s exposure, over only
  part of the field; consistent with this being the one frame among the
  five with an invalid (negative) `MOONILLF` value in its own archived
  metadata, a plausible sign of degraded weather telemetry at that time.
  Not independently confirmed (no all-sky camera or cloud-sensor log was
  checked), but every alternative explanation checked was ruled out by
  direct measurement.

  `run_pipeline`'s `quality_max_std` gate (default `1.5`, using plain
  `Statistics.std` — see its docstring for why a MAD-based spread,
  robust to exactly this kind of minority-of-pixels outlier, turned out
  blind to this specific anomaly and was rejected) does correctly
  exclude this frame on real data (`frame_std=2.025 > 1.5`, confirmed by
  the actual warning `run_pipeline` emits), and both known objects are
  still recovered afterward. But `link_candidates` requires every frame
  to match by default, so a gated frame — contributing zero detections —
  makes no tracklet reachable at all unless `min_frames` is lowered to
  account for it (`examples/real_data_demo.jl` uses
  `length(fits_paths) - 1`); loosening that constraint by one frame,
  independent of anything about data quality, roughly doubled the
  tracklet count again (334 → 667) purely through looser matching
  combinatorics. So the gate is confirmed working, but the net real-data
  effect on tracklet count is now confounded by that combinatorial
  change rather than cleanly attributable to data quality — a clean
  before/after comparison would need to hold `min_frames` fixed some
  other way (e.g. by literally excluding the bad frame from the input
  list rather than gating it internally), not attempted here.

## Example: real data

`examples/real_data_demo.jl` runs the pipeline against real ZTF (Zwicky
Transient Facility) frames both with and without ZOGY differencing, and
cross-matches both against SkyBoT — a controlled comparison, not just a
demonstration. Fetch the data first (public, no authentication required):

```
examples/fetch_data.sh
julia --project=. examples/real_data_demo.jl
```

Building the reference stack (30 frames, each individually reprojected)
is the slow part — tens of minutes on a laptop, one-time per run.

## Rotation period recovery

For a confirmed discovery, given a dedicated photometric follow-up
sequence (many exposures over hours, at a fixed sky position — the
target should barely move between them, unlike the original discovery
epochs):

```julia
times, flux, flux_err = light_curve(fits_paths, ra, dec)
result = recover_rotation_period(times, flux; minimum_period=0.02, maximum_period=1.0)
result.period, result.false_alarm_probability
```

`minimum_period`/`maximum_period` bound the search (same units as
`times`, i.e. days) and should bracket the rotation periods physically
plausible for the object's size class. A small `false_alarm_probability`
is what distinguishes a real periodic signal from a noise fluctuation —
see the function's docstring.

## Plate-solving

For a frame with no WCS already in its header, and a
[nova.astrometry.net](https://nova.astrometry.net/) API key (free
registration):

```julia
run_pipeline(fits_paths; reference=reference, plate_solve_api_key=key)
```

or directly: `plate_solve(fits_path; api_key=key)`. This is a live
network round trip — upload, then poll until the frame solves — so it is
slow and requires connectivity; see the **Known limitations** entry above
before relying on it, since this project has not yet run it against the
real service end to end.

## Using real IASC campaign data

Not attempted in this project — real campaign access needs the user's
own IASC registration, not something this pipeline can fetch on its own
(unlike the public ZTF demo data above). Once campaign FITS files are in
hand:

- Point `run_pipeline` (or `examples/real_data_demo.jl`'s pattern) at the
  local file paths directly; no fetch script is needed for files you
  already have.
- Check `timestamp_key` and whether the frames already carry a WCS before
  assuming the `"MJD-OBS"` default and `plate_solve_api_key=nothing`
  (unset) both apply — genuinely unknown until real files are in hand,
  not verified against this codebase.
- Re-tune `threshold`, `match_radius`, and (if using the ZOGY path)
  `quality_max_std` for the new data the same way `examples/real_data_demo.jl`
  did for ZTF, rather than assuming the current defaults — calibrated
  against one specific survey's noise characteristics — transfer.
- Run the existing test suite first (`Pkg.test()`) to confirm the
  environment itself is sound, then adapt `examples/real_data_demo.jl` as
  a validation template: known objects in the field (via `crossmatch_catalog(...; :skybot)`)
  are the same kind of ground truth used there.

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

A reproducible development environment is also provided via Nix
(`flake.nix`):

```
nix develop
```

## Dependencies

- [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) — FITS I/O
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) — source detection and photometry
- [WCS.jl](https://github.com/JuliaAstro/WCS.jl) — astrometric (pixel-to-sky) calibration
- [Reproject.jl](https://github.com/JuliaAstro/Reproject.jl) — resampling frames onto a common pixel grid for reference stacking
- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) — Fourier-domain ZOGY difference imaging
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) — sub-pixel PSF stamp alignment
- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) — rotation-period periodogram analysis
- [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl) — catalog cross-match queries
- [JSON.jl](https://github.com/JuliaIO/JSON.jl) — nova.astrometry.net API requests (plate-solving)

## License

MIT. See `LICENSE`.
