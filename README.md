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

Planned extensions: variable-star and transient detection, and rotation
period recovery via Lomb-Scargle periodograms for confirmed discoveries.

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

Four real bugs surfaced only by testing against real ZTF data (not by the
synthetic tests or code review), all now fixed with regression tests:

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

A fifth issue was a dataset-selection mistake rather than a code bug: the
first real-data demo used a field/night where the faintest catalogued
asteroid was 1.2 mag *below* that night's own detection limit — no
algorithm can recover a signal that was never above the noise, so the
demo now picks a night with a known object bright enough to serve as
actual ground truth.

With both fixed, `examples/real_data_demo.jl`'s actual result on that
night (field 451, 2019-10-23): the undifferenced baseline finds 129
tracklets and recovers both known objects in the field (2002 UY45, 1997
KO3); ZOGY also recovers both, at consistent sky offsets (confirming the
subtraction is correctly calibrated, not just "not obviously broken"),
but finds 274 tracklets total and no *additional* known object. The
likely reason ZOGY doesn't show its expected depth advantage here: both
known objects in this field (Mv 18.4, 19.6) are bright enough that the
undifferenced baseline already recovers them trivially — this dataset
doesn't happen to contain a known object faint enough to sit below a
single frame's noise floor but above the deep reference's, which is the
specific regime ZOGY is for. The near-doubled tracklet count is not yet
explained (more real faint candidates, or more differencing artifacts at
this threshold — most likely some of both) and is flagged here rather
than glossed over.

## Known limitations

- **No plate-solving fallback.** `load_wcs` only parses a WCS solution
  already present in the FITS header; it cannot derive one from an
  unsolved frame. IASC campaign frames are generally pre-solved, so this
  is not currently blocking, but it will need addressing before the
  pipeline can be used on frames from other sources (e.g. own
  blazar/exoplanet-timing imaging). See the `TODO` on `load_wcs` in
  `src/astrometry.jl` for candidate approaches.
- **Empirical PSF, not a fitted model.** `estimate_psf` stacks real star
  cutouts rather than fitting an analytic profile (Gaussian/Moffat), which
  keeps it survey-agnostic but means its quality depends on having enough
  bright, isolated, unsaturated stars in the frame.
- **`zogy_subtract`'s astrometric-noise term is opt-in.** Passing
  `n_sources`/`r_sources` estimates residual sub-pixel registration error
  from matched star positions; without them it is treated as zero, which
  is only accurate if reprojection has already removed essentially all of
  it (true here, since `build_reference` reprojects onto the exact
  science-frame grid, but not a safe assumption in general).
- **`link_candidates`'s velocity comes from frames 1 and 2 only** (see its
  docstring) — closely spaced first frames amplify velocity error when
  extrapolated across a longer baseline. A robust fit across all frames
  would remove the sensitivity to frame spacing and ordering.
- **ZOGY's real-data tracklet count is not yet understood.** On the
  `examples/real_data_demo.jl` field/night, ZOGY produces roughly twice as
  many tracklets as the undifferenced baseline (274 vs 129), with no
  increase in known-object recovery. Not yet root-caused — plausibly more
  genuine faint candidates, more differencing artifacts at the current
  threshold, or both.

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
- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) — periodogram analysis
- [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl) — catalog cross-match queries

## License

MIT. See `LICENSE`.
