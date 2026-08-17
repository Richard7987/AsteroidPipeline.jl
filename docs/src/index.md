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

`find_variable_sources` searches the same per-frame detections for
**stationary**, flux-varying sources — variable stars and transients,
as opposed to `link_candidates`'s moving-object search — matching
positions across frames instead of a linear motion model, then filtering
by `variability_chi2` against a constant-flux null hypothesis.
`search_field` runs both searches from a single shared detection pass
over `fits_paths` (avoiding the cost of detecting twice), returning
`(movers, variables)`; use `run_pipeline` alone when only asteroid
candidates are needed.

For a confirmed discovery, `light_curve` (forced aperture photometry at a
fixed sky position across a dedicated follow-up sequence) and
`recover_rotation_period` (a Lomb-Scargle periodogram over that light
curve, via `LombScargle.jl`) recover a rotation period — a separate,
optional follow-up step, not part of `run_pipeline` itself. The same
periodogram applies directly to a `find_variable_sources` candidate's own
`(frame, flux)` points, for periodic variables.

## Status

Early development. `detect_sources`, `link_candidates`, the WCS
calibration step, and `crossmatch_catalog` are implemented and wired
together end to end in `run_pipeline`, validated against synthetic FITS
frames with a known injected source track, and exercised against real
public survey data (see `examples/real_data_demo.jl`) and real IASC
practice campaign data (see `examples/iasc_demo.jl` and the "Using real
IASC campaign data" section below — 9 real, independently-catalogued
objects recovered across 5 real Pan-STARRS1 fields).

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

`find_variable_sources`/`search_field` (stationary, flux-varying source
detection) and `fit_moffat_psf` (`estimate_psf`'s analytic-PSF fallback)
are implemented, tested against synthetic data, and — for
`find_variable_sources`'s photometric normalization, S/N floor, and
`chi2_threshold` default, and for `fit_moffat_psf`'s recovered PSF
width — calibrated directly against real ZTF data (see the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#detect_sources's-flux-has-been-silently-wrong-since-the-beginning:-a-transposed-aperture), including a real bug in
`detect_sources`'s own flux measurement this calibration work found and
fixed). Now also validated end to end, via `search_field`, against a real
field with an independently-confirmed variable star as ground truth: ZTF
field 487/CCD 12/quadrant 1/zr, night 2019-06-10, containing ASASSN-V
J183620.31 (VSX: type EW, period 0.322 d) — see
`examples/variable_star_demo.jl` and the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#Validating-search_field-against-a-real,-independently-confirmed-variable-star)
for the full run. The target was recovered (crossmatched against VSX at a
1.7" offset) among 22 variable candidates from 638 detections; a
Lomb-Scargle period fit to that single partial night (only ~32% of one
0.322 d period) found a real, highly significant periodic signal
(FAP ≈ 0) but at 0.5 d, not the true period — an expected outcome of the
partial phase coverage, declared before running rather than adjusted
after seeing it, and not treated as a failure of `find_variable_sources`
itself (which is what the crossmatch recovery actually validates).

On that real dataset (field 451, 2019-10-23), the undifferenced baseline
finds 133 tracklets and recovers both known objects in the field (2002
UY45, 1997 KO3); ZOGY also recovers both, at consistent sky offsets
(confirming the subtraction is correctly calibrated), but finds 667
tracklets total and no *additional* known object — both known objects
here are bright enough that the baseline already recovers them trivially,
so this dataset doesn't exercise ZOGY's actual advantage (recovering
objects below a single frame's noise floor). The raw tracklet-count gap
is not a clean read on ZOGY's noise properties; see the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#The-quality-gate's-combinatorial-side-effect-on-tracklet-count) for why, and for the full
record of every real bug this project's real-data testing surfaced so
far — most recently three more from the real IASC campaign work below —
all fixed with regression tests, and how each was diagnosed.

## Known limitations

- **Empirical PSF's quality still depends on the field, though less than
  it used to.** `estimate_psf` stacks real star cutouts, which captures
  the true PSF shape (wings included) without fitting a model family per
  instrument, but needs enough bright, isolated, unsaturated stars to do
  it. It now retries at a progressively relaxed `min_separation` (halved,
  up to `relaxation_attempts` times, default 2) before giving up — a
  field with a few usable stars at a tighter isolation radius still gives
  the real PSF shape, which the analytic fallback never can. Only once
  even the most relaxed attempt finds nothing does it fall back (by
  default) to `fit_moffat_psf` — a parametric Moffat fit — rather than
  failing outright; the fallback trades exact PSF shape for robustness,
  and is not itself a substitute for a genuinely well-behaved field.
- **`zogy_subtract`'s astrometric-noise term (`V_ast`) is still opt-in at
  the `zogy_subtract` level** — it needs `n_sources`/`r_sources` passed
  explicitly, and is `0` without them; this is deliberate API layering
  (`run_pipeline` always supplies them, so a direct caller who doesn't
  need the extra `detect_sources` cost can skip it), not something
  planned to change. What did change: leaving both unset now emits a
  `@warn` (once per session), so omitting `V_ast` is a visible choice
  instead of a silent default a direct caller could miss.
- **`find_variable_sources` still has a real, measured false-positive
  floor on real single-epoch aperture photometry, even after fixing it
  once.** The first hypothesis tried — pixel-grid jitter in
  `detect_sources`'s peak-pixel aperture centering — turned out to be a
  real but secondary effect: refining the aperture to a sub-pixel
  centroid (see `detect_sources`'s docstring) barely moved the
  false-positive rate (23% → 22% at `chi2_threshold=3` on real ZTF field
  451). The actual dominant cause, found by checking *which* stars were
  flagged, was a systematic photometric error floor — bright stars'
  tiny formal errors made ordinary flat-fielding/PSF-variation
  systematics look like huge chi2 significance. Adding that floor
  (`variability_chi2`'s `systematic_error_fraction`) cut the rate: a 1%
  floor took it to 6% at threshold 3, ~0% at threshold 20; sweeping the
  floor further (against the same real stars, and checked against a real
  confirmed variable — ASASSN-V J183620.31 — to make sure real
  sensitivity wasn't sacrificed for it) found more room, without giving
  up real detections: the current default, 2%, cuts the
  `chi2_threshold=10.0` rate to 2.0% on this dataset (vs 1%'s 3.3%),
  while that confirmed variable still clears the threshold with a 3x
  margin — see `find_variable_sources`'s
  docstring and the [Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#The-centroid-fix-barely-moved-the-false-positive-floor-—-the-real-cause-was-a-systematic-error-floor) for the
  full before/after numbers. Still not zero: treat a candidate as needing
  independent confirmation (a catalog match or a recovered period), not
  as self-evidently real.

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
using AsteroidPipeline

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
using AsteroidPipeline

run_pipeline(fits_paths; reference=reference, plate_solve_api_key=key)
```

or directly: `plate_solve(fits_path; api_key=key)`. This is a live
network round trip — upload, then poll until the frame solves — so it is
slow and requires connectivity. Validated against the real service: see
the [Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#plate_solve-validated-end-to-end-against-the-live-service).

## Using real IASC campaign data

Now attempted, against 5 real Pan-STARRS1 (PS1) IASC practice sets
("Practice Image Sets", 2019-08-28/09-04/09-24), each 4 exposures of the
same field over ~40-70 min — see `examples/iasc_demo.jl` (point it at
your own local practice/campaign FITS; IASC material isn't public, so
unlike the ZTF demo above there is no fetch script). `run_pipeline`
recovered **9 real, independently-catalogued objects** across the 5
fields via `crossmatch_catalog(...; :skybot)` — including a Jupiter
Trojan (2019 NB9) — the first end-to-end validation of this pipeline
against real IASC-style data, not just ZTF.

Getting there surfaced four real, fixed issues — see the
[Investigation Log](https://richard7987.github.io/AsteroidPipeline.jl/dev/investigation-log#Validating-against-real-IASC-Pan-STARRS1-campaign-data)
for the full story of each:

- `load_wcs` raised "Linear transformation matrix is singular" on every
  one of these real headers — PS1's legacy `CNPIX1`/`CNPIX2` keywords
  make wcslib build a separate, degenerate implicit WCS alongside the
  header's real, valid one. Fixed in `load_wcs` itself.
- These FITS files mark invalid/masked pixels via the standard `BLANK`
  keyword rather than `NaN`, which `FITSIO.jl` doesn't auto-convert —
  left alone, `detect_sources` read one masked region as real flux and
  produced 158,443 spurious detections in a single frame. Handled as a
  preprocessing step in `examples/iasc_demo.jl` (not in `src/`, since
  this is a real-FITS-ingestion concern, not `detect_sources`'s job).
- `crossmatch_catalog(...; :skybot)` queried one candidate at a time,
  fully sequentially — real candidate lists here (hundreds to
  thousands of tracklets) took minutes to hours, and a multi-hour run
  eventually died to a transient connection error with no retry. Fixed
  in `_crossmatch_skybot` itself: concurrent requests (real, measured
  ~8x wall-clock speedup) and a retry on transient HTTP errors.
- `match_radius`, first converted from `real_data_demo.jl`'s ZTF value to
  keep the same ~10" angular tolerance, turned out far looser than PS1's
  own real astrometric precision (`PERROR`, in these headers: 0.20-0.23")
  — on the densest field this produced 10,422 tracklets, almost all
  spurious duplicates of the same real objects (distinct real stars
  within 10" of each other, or the same object matched by several
  near-identical trial velocities). Retuned to 2" (~10x `PERROR`,
  measured from the headers, not guessed) and confirmed directly: the
  same 9 distinct known objects are still recovered in every field,
  while total tracklets across all 5 fields drop from 16,158 to 4,960
  (-69%) — this was cleanup of spurious duplicates, not lost detections.

General guidance for pointing this pipeline at other real campaign data:

- Point `run_pipeline` (or `examples/iasc_demo.jl`'s pattern) at the
  local file paths directly; no fetch script is needed for files you
  already have.
- Check `timestamp_key` and whether the frames already carry a WCS before
  assuming the `"MJD-OBS"` default and `plate_solve_api_key=nothing`
  (unset) both apply — PS1's headers happened to match `"MJD-OBS"`
  exactly, but that's this survey, not a general guarantee.
- Re-tune `threshold`, `match_radius` (see the `PERROR`-based lesson
  above — scale to the *survey's own* astrometric precision, not another
  survey's pixel scale), and (if using the ZOGY path) `quality_max_std`
  for the new data, rather than assuming defaults calibrated against
  ZTF/PS1 transfer as-is.
- Run the existing test suite first (`Pkg.test()`) to confirm the
  environment itself is sound, then adapt `examples/iasc_demo.jl` or
  `examples/real_data_demo.jl` as a validation template: known objects in
  the field (via `crossmatch_catalog(...; :skybot)`) are the same kind of
  ground truth used there.

## Dependencies

- [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) — FITS I/O
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) — source detection and photometry
- [WCS.jl](https://github.com/JuliaAstro/WCS.jl) — astrometric (pixel-to-sky) calibration
- [Reproject.jl](https://github.com/JuliaAstro/Reproject.jl) — resampling frames onto a common pixel grid for reference stacking
- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) — Fourier-domain ZOGY difference imaging
- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) — sub-pixel PSF stamp alignment
- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) — rotation-period periodogram analysis
- [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) — analytic (Moffat) PSF fallback fitting
- [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl) — catalog cross-match queries
- [JSON.jl](https://github.com/JuliaIO/JSON.jl) — nova.astrometry.net API requests (plate-solving)

## License

MIT. See `LICENSE`.
