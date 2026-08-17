# Investigation log

Chronological record of what real-data testing found and how each finding
was diagnosed — kept separate from the main wiki pages so those stay
focused references rather than narratives. Cross-referenced from
`docs/src/index.md`'s Known Limitations section and from the relevant
docstrings.

This page covers the original ZTF field 451 investigation, chronologically
first. Later validation work against other real datasets got large enough
to split onto their own pages:

- [Validating `search_field` against a real, independently-confirmed variable star](variable-star-validation.md)
- [Validating against real IASC (Pan-STARRS1) campaign data](iasc-campaign-validation.md)
- [Revisiting the two intentional-design "limitations"](design-refinements.md)

## Real ZTF data surfaced four bugs no synthetic test or code review caught

All fixed, all with regression tests added afterward:

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

## Dataset-selection mistake: the first real-data demo picked an unrecoverable field

Not a code bug: the first real-data demo used a field/night where the
faintest catalogued asteroid was 1.2 mag *below* that night's own
detection limit. No algorithm can recover a signal that was never above
the noise. The demo now picks a night with a known object bright enough
to serve as actual ground truth (see `examples/real_data_demo.jl`'s
header comment for how that field/night/reference set was chosen).

## PSF-timing and astrometric-noise bugs behind ZOGY's excess detections

Investigating why ZOGY produced ~2x more tracklets than the undifferenced
baseline on real ZTF data (274 vs 129 at the time), with `std(S_corr)`
not actually ~1 as the model assumes and 85-97% of the excess detections
landing within 15 px of a bright reference star, traced it to two
compounding issues:

- `estimate_psf` measured each frame's PSF *before* reprojection, but
  `zogy_subtract` differenced the frame *after*. Reprojection's
  interpolation subtly reshapes a PSF, and that mismatch left systematic
  subtraction residuals at every bright star. Fixed by measuring the PSF
  (and each frame's own detected sources) on the same post-reprojection
  pixels `zogy_subtract` actually uses, and feeding those sources through
  as `n_sources`/`r_sources` so `zogy_subtract`'s astrometric-noise term
  (`V_ast`) can price in residual misregistration instead of ignoring it.
- The reprojected border's NaN was filled with `0.0` rather than the
  valid region's own background level (see the bug list above).

On real data this brought `std(S_corr)` from ~2.0-2.6 down to ~1.1-1.2 in
4 of 5 test frames. The 5th stayed elevated and produced *more* bright
sources than before rather than fewer — investigated separately below.

## Root-causing the anomalous 5th frame: a likely passing cloud

GAIN, SEEING, MAGZP, whole-frame background level, saturated-pixel count,
and airmass were all checked directly against the other four frames and
none stood out. Differencing this frame's *raw* pixels against another
frame's, bypassing all ZOGY machinery, found the actual cause: ~1750
pixels with a significant, one-sided excess (no matching negative lobe,
so not a registration/dipole artifact), ~98% confined to the bottom ~30%
of the detector. A control pair of two normal frames, differenced the
same way, found only ~80-110 such pixels, distributed in proportion to
where the frame's actual stars are — not concentrated in one region. The
whole-frame background stayed flat across that region (no gradient),
ruling out amplifier glow or vignetting.

This pattern — real stars showing localized excess halos in one part of
the frame, flat background, no positional offset — is the signature of
thin cloud scattering starlight during the 30 s exposure, over only part
of the field, not a code or calibration bug. It also happens to be the
one frame among the five with an invalid (negative) `MOONILLF` value in
its own archived metadata, consistent with degraded weather telemetry at
that time. Not independently confirmed against an all-sky camera or cloud
sensor, but every alternative explanation checked was ruled out by direct
measurement rather than left as a guess.

## Quality gate: MAD looked more rigorous, but was blind to the real anomaly

A per-frame quality gate (`run_pipeline`'s `quality_max_std`) was added to
catch frames like the one above. First attempt used a MAD-based
(median-absolute-deviation) robust spread instead of plain
`Statistics.std`, specifically to avoid a real bright source skewing
`std` in a small *synthetic test* image (confirmed: at 100x100 px, a
genuine ~60 sigma injected source pushed a clean frame's std to ~4).

Against real data this backfired: the actual anomalous frame's excess
showed up as ~232 point-like bright residuals, not a bulk noise increase
— exactly the kind of minority-of-pixels outlier MAD is *designed* to
ignore. Its MAD measured ~0.997, indistinguishable from the good frames'
~0.98-1.10, while its plain `std` (2.025) stood clearly apart from the
good frames' (1.10-1.18). Reverted to plain `std`; fixed the synthetic
test's false positive by using a larger image (600x600) instead of
changing the statistic, since on real ~10^6 px data a genuine source's
few dozen pixels are too small a fraction to skew `std` regardless.

Confirmed against real data: the gate (default `quality_max_std=1.5`)
correctly excludes the anomalous frame (`frame_std=2.025 > 1.5`), and
both known SkyBoT objects are still recovered afterward.

## The quality gate's combinatorial side effect on tracklet count

`link_candidates` requires every frame to match by default
(`min_frames`). A gated-out frame contributes zero detections, so once
the quality gate correctly excludes one, no tracklet is reachable at all
unless `min_frames` is lowered to account for it —
`examples/real_data_demo.jl` uses `length(fits_paths) - 1`.

That one-frame relaxation, independent of anything about data quality,
roughly doubled the tracklet count again (334 → 667) purely through
looser matching combinatorics: requiring 4-of-4 frames to match is a
weaker constraint than requiring 5-of-5, so more spurious coincidental
matches pass regardless of how clean the underlying data is. The gate
itself is confirmed working (see above), but this means the real-data
tracklet-count comparison in `README.md` is *not* a clean read on ZOGY's
noise properties — the `min_frames` change confounds it. A clean
before/after comparison would need to hold `min_frames` fixed some other
way (e.g. by literally excluding the bad frame from the input list rather
than gating it internally), not attempted here.

## `plate_solve` validated end-to-end against the live service

nova.astrometry.net's own API docs require a `Referer:
https://nova.astrometry.net/api/login` header on programmatic file
downloads, as an anti-scraper-bot check — `plate_solve` didn't set it.
Added it to every GET request in `src/platesolve.jl` (submission-status
poll, job-status poll, WCS fetch); also switched the base URLs from
`http://` to `https://` to match the documented endpoints.

With a real API key: a synthetic image containing one fabricated point
source in noise (no genuine star pattern) correctly failed to solve,
astrometry.net reporting `status: "failure"` rather than returning a
wrong answer — the expected outcome, since plate-solving works by
matching real star-pattern asterisms against a sky index, and there was
no real pattern to match. A real ZTF frame (`data/real/science/`, its own
existing WCS ignored) solved successfully: recovered
`CRVAL ≈ (36.3997, 2.0038)`, consistent with that field's known centre
(~36.5, ~2.1). One upload attempt hit a transient `503` from the
service; a retry succeeded, so `plate_solve` callers should be prepared
to retry on transient server errors rather than assume a single failed
upload means the service is unusable.

## `detect_sources`'s flux has been silently wrong since the beginning: a transposed aperture

Building `find_variable_sources` required trusting `detect_sources`'s
`flux` for the first time in a *quantitative* way — every earlier use
(`link_candidates`, tracklet building, WCS calibration, cross-matching)
only depends on its `x`/`y` positions. That new dependency surfaced a bug
that had been present, silently, since `detect_sources` was first written:
its measured flux was wrong by orders of magnitude for almost any real
source.

Root cause: `Photometry.jl` (v0.9.8) is internally inconsistent between
its own two pieces. `PeakMesh`'s `extract_sources` reports positions in
the standard Cartesian sense — `to_nt(ci) = (x=ci[2], y=ci[1], ...)`, i.e.
`x` is the array's *second* dimension (column), `y` its *first* (row).
`CircularAperture`, in the same package, does the opposite internally —
`bounds`/`overlap` treat its own `.x` field as indexing the array's
*first* dimension and `.y` the second. `detect_sources` built
`CircularAperture(row.x, row.y, aperture_radius)` directly from
`PeakMesh`'s output, so every aperture was centred at the transposed
pixel — correct only when a source's row and column indices happened to
coincide (the diagonal), or invisibly wrong on a square canvas at a
generic position, or (on a non-square canvas, or once truly out of the
transposed array's bounds) landing on pure background instead.

Found by comparing `detect_sources`'s measured flux against the closed-form
enclosed-energy integral for an isolated, noise-free synthetic Gaussian on
a non-square canvas at an asymmetric position: measured flux came back
~140x too small. Swapping the aperture's constructor arguments
(`CircularAperture(row.y, row.x, aperture_radius)`) brought it to within
~1.4% of the analytic value — geometric quantization error, not a
remaining bug. `light_curve` (`src/rotation.jl`) was checked the same way
and is *not* affected: it deliberately never permutes its image (see its
own docstring), and `WCS.world_to_pix`'s returned `(x, y)` already matches
raw FITS `(NAXIS1, NAXIS2)` order — which happens to be exactly what
`CircularAperture` expects internally, by coincidence of two conventions
cancelling out rather than by any intentional match.

This had zero effect on every result validated so far in this project —
`run_pipeline`'s real-data tracklet counts (133 baseline / 667 ZOGY, both
recovering the same 2 known objects) depend only on detected *positions*,
never on `flux` — but it fully invalidated the first real-data
measurements made *while building* `find_variable_sources` earlier the
same session, before this was found: a claimed 29.4% median
`flux_err/flux` (actual, post-fix: **2.56%**), a claimed "only 2 of 120
stars clear a 10% S/N floor" (actual: **119 of 119**), and a claimed
8-40% ensemble-vs-`MAGZP` mismatch explained by low-S/N selection bias
(the mismatch is real, but stable at 9-13% with or without an S/N cut —
not a selection effect at all; see below). Every constant and docstring
claim built on those numbers was rewritten against the corrected
measurements rather than left standing.

## What was actually driving the ~10 percent photometric-scale mismatch, once flux was measured correctly

With flux measured correctly, `photometric_scale`'s ensemble ratio against
each frame's own `MAGZP` zeropoint still disagreed by 9-13% — but now
*independent* of the S/N cut, ruling out selection bias as the cause.
Checked against each frame's `SEEING` header value (1.805-2.009 px across
the 5 real frames): the frames with better (smaller) seeing than frame 1
all showed `photometric_scale < 1`, consistent with the standard "aperture
correction" effect — a fixed-radius aperture (`aperture_radius`, default 3
px, comparable to ZTF's own seeing) encloses a larger fraction of a star's
total flux when the PSF is more concentrated. `photometric_scale`'s
ensemble-differential approach corrects for this automatically (it only
needs frame-to-frame consistency, not a causal model), so no code change
was needed here — only the docstring's explanation, which had cited the
now-debunked selection-bias story.

## `find_variable_sources`'s chi-squared test has a real, measured false-positive floor on real data

Even with correct flux and photometric normalization, real ZTF data's
reduced-chi2 distribution over 119 matched stationary stars has a heavy
tail: 23% exceed a threshold of 3 (the textbook-reasonable default),
13% still exceed 20. The likely cause: `detect_sources` positions each
detection at `PeakMesh`'s integer peak pixel, not a sub-pixel centroid, so
a 1-pixel jitter between frames — from noise, or a slightly different PSF
realization — against a small aperture (3 px, close to the PSF's own
scale) produces a real, but spurious, flux swing frame to frame. Raising
`chi2_threshold`'s default to 50 cuts the false-positive rate to 8% on
this dataset — better, but not clean, and documented as such in
`find_variable_sources`'s own docstring rather than presented as solved.
The actual fix (forced, sub-pixel-centroided photometry instead of
peak-position aperture photometry) is a larger change, not attempted here.

## The centroid fix barely moved the false-positive floor — the real cause was a systematic error floor

The sub-pixel-centroiding fix flagged as the likely cause above was
implemented (`detect_sources` now refines each `PeakMesh` integer peak to
a flux-weighted first-moment centroid — within a `radius`-sized window,
falling back to the untouched integer position whenever the refinement
isn't trustworthy — before centering the aperture there; the *returned*
`x`/`y` columns are deliberately left as the original integers, so no
downstream position-dependent logic or exact-equality test is affected).
Re-measuring the same reduced-chi2 distribution on the same 119 matched
stars: 23% → 22% at threshold 3, 13% → 13% at threshold 20. A real
change, but not a meaningful one — reported as such rather than declared
a fix.

Checking *which* stars were actually being flagged settled it: not the
faintest ones, as low-S/N selection would predict, but the **brightest**
— median relative flux error 0.25% among flagged stars vs. 2.93% among
the rest. That is the textbook signature of a systematic error floor
(imperfect flat-fielding, frame-to-frame PSF variation — real effects
that `flux_err`'s pure Poisson/background model never captures), not
underestimated statistical noise or pixel-grid jitter. Adding a
systematic floor in quadrature (`variability_chi2`'s
`systematic_error_fraction`, default 1% — a standard value in
forced-photometry pipelines, not tuned to this dataset) is what actually
fixed it: 23% → 6% at threshold 3, 13% → ~0% at threshold 20.
`chi2_threshold=10.0` with that floor gives a 4% false-positive rate on
this dataset — set as the new default, replacing the old
threshold-50/8%-floor compromise.

## `min_frames` now auto-adjusts for quality-gated frames

The combinatorial side effect documented above (a `quality_max_std`-gated
frame silently making no tracklet reachable at all under the default
`min_frames`) previously required `examples/real_data_demo.jl` to pass
`min_frames = length(fits_paths) - 1` by hand on the ZOGY path.
`_detect_all_frames` now reports how many frames it gated (`n_gated`,
always `0` on the raw path), and `run_pipeline`/`search_field`'s
`min_frames` (and `search_field`'s `variability_min_frames`) default to
`nothing`, resolved internally to `length(fits_paths) - n_gated` — an
explicitly-passed value still overrides this exactly as before. Verified
against the real dataset: `examples/real_data_demo.jl`'s manual override
was removed, and the pipeline still recovers the same 133/667 tracklets
and both known objects on that field, now with no caller-side workaround.

## Batching `crossmatch_catalog(...; :vsx)`/`(...; :simbad)`

Migrating off the CDS X-Match service (see above) left `_crossmatch_simbad`/
`_crossmatch_vsx` querying one candidate per TAP request — fine for a
handful of candidates, not for a real candidate list of hundreds. Two
obvious batching routes were tried against the live services and rejected
before settling on a third: a TAP `UPLOAD`-based batch (one request for
any N) — SIMBAD accepts the multipart upload but fails to resolve the
uploaded table's columns; VizieR rejects anything that isn't a full
VOTable document. ADQL `UNION` — rejected outright by SIMBAD's parser
("UNION is not supported in ADQL"). What does work, confirmed against a
real multi-candidate query: `OR`-chaining one `CONTAINS(...) = 1`
cone-search clause per candidate into a single request, then resolving
client-side (via a haversine great-circle distance, since a single
`OR`-chained query has no per-row candidate tag) which candidate each
returned row matches. Batched 50 candidates per request — conservative,
not measured at higher N on these free, shared, anonymous-use services.
Cuts N requests to `ceil(N / 50)`.

## Tuning `find_variable_sources`'s systematic error floor past the first value that worked

The 1% systematic error floor (see above) was the first value tried,
chosen because it's a standard number in forced-photometry pipelines, not
because it was shown to be optimal. With a real positive control now in
hand (ASASSN-V J183620.31, recovered by `find_variable_sources` earlier
this session — see its own
[validation page](variable-star-validation.md)), the floor could finally
be checked from *both* sides of the tradeoff at once, not just the
false-positive side: sweeping `systematic_error_fraction` against the
same 152 real, matched, high-S/N stationary stars from ZTF field 451 (a
slightly different count than the "119" quoted earlier — this sweep
additionally required `min_frames` at the full 5-frame count and
`normalize=true` together, narrowing the matched set), the
`chi2_threshold=10.0` false-positive rate dropped from 3.3% at a 1%
floor to 2.0% at 2%, 1.3% at 3%, and 0% at 5%. Naively, that argues for
as high a floor as possible — but a floor this large also suppresses
*real* variability, and that side had never been checked. Running the
same sweep against ASASSN-V J183620.31's own real forced-photometry
light curve: reduced chi2 falls from 123 (at 1%) to 31 (at 2%) to 13.8
(at 3%) to 5.0 (at 5%) — the last of which drops *below*
`chi2_threshold=10.0`, meaning a 5% floor would have made this exact,
real, independently-confirmed variable star invisible to
`find_variable_sources`. Settled on 2%: comfortably above 1%'s
false-positive rate, while leaving the real variable's signal a 3x margin
over threshold — the largest floor checked that doesn't cost real
detections, not the smallest false-positive rate achievable.

![False-positive rate on 152 real stationary stars, and reduced chi2 on a real confirmed variable, swept across systematic_error_fraction values](assets/systematic-error-floor-sweep.png)

## See also

- [Validating `search_field` against a real, independently-confirmed variable star](variable-star-validation.md)
- [Validating against real IASC (Pan-STARRS1) campaign data](iasc-campaign-validation.md)
- [Revisiting the two intentional-design "limitations"](design-refinements.md)
