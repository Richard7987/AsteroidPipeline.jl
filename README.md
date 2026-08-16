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

On that real dataset (field 451, 2019-10-23), the undifferenced baseline
finds 133 tracklets and recovers both known objects in the field (2002
UY45, 1997 KO3); ZOGY also recovers both, at consistent sky offsets
(confirming the subtraction is correctly calibrated), but finds 667
tracklets total and no *additional* known object — both known objects
here are bright enough that the baseline already recovers them trivially,
so this dataset doesn't exercise ZOGY's actual advantage (recovering
objects below a single frame's noise floor). The raw tracklet-count gap
is not a clean read on ZOGY's noise properties; see
[`INVESTIGATION_LOG.md`](INVESTIGATION_LOG.md) for why, and for the full
record of every real bug this project's real-data testing surfaced (five
so far, all fixed with regression tests) and how each was diagnosed.

## Known limitations

- **Empirical PSF, not a fitted model.** `estimate_psf` stacks real star
  cutouts rather than fitting an analytic profile (Gaussian/Moffat), which
  keeps it survey-agnostic but means its quality depends on having enough
  bright, isolated, unsaturated stars in the frame.
- **`zogy_subtract`'s astrometric-noise term (`V_ast`) is opt-in at the
  `zogy_subtract` level** — it needs `n_sources`/`r_sources` passed
  explicitly, and is `0` without them. `run_pipeline` always supplies
  them, so this only matters when calling `zogy_subtract` directly.
- **A quality-gated frame silently tightens `link_candidates`.**
  `run_pipeline`'s `quality_max_std` (default `1.5`) excludes a frame
  whose `S_corr` spread is too high (confirmed against real data — see
  [`INVESTIGATION_LOG.md`](INVESTIGATION_LOG.md)), but a gated frame
  contributes zero detections, and `link_candidates` requires every frame
  to match by default. Pass a lower `min_frames` (e.g.
  `length(fits_paths) - 1`) when using the ZOGY path, or no tracklet will
  ever be reachable if any frame gets gated — `examples/real_data_demo.jl`
  does this.

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
slow and requires connectivity. Validated against the real service: see
[`INVESTIGATION_LOG.md`](INVESTIGATION_LOG.md).

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
