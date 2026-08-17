# Investigation log

Chronological record of what real-data testing found and how each finding
was diagnosed — kept separate from `README.md` so the README stays a
focused reference rather than a narrative. Cross-referenced from
`README.md`'s Known limitations section and from the relevant docstrings.

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

## Validating `search_field` against a real, independently-confirmed variable star

Every real-data check so far had confirmed known *moving* objects (via
SkyBoT) but never a known *variable* — the one real dataset checked had
only one catalogued VSX variable in its footprint, too faint to serve as
a useful positive control. Fixed by choosing a real target with an actual
VSX catalog entry: ASASSN-V J183620.31 (type EW, a contact eclipsing
binary, period 0.322427 d), in ZTF field 487/CCD 12/quadrant 1/zr, night
2019-06-10 — a real high-cadence campaign with 144 exposures over 2.45 h,
thinned to 29 for `examples/variable_star_demo.jl`. Declared in advance:
2.45 h covers only ~32% of one period, so full period recovery from this
single night was not expected — the actual test was whether
`find_variable_sources` flags the star as variable at all.

First attempt did not run: `search_field`'s default `threshold=5.0` (used
at `8.0` in `real_data_demo.jl`) produced ~12,900 detections per frame
here, against field 451's ~130 — this field sits near the galactic plane.
`link_candidates`'s tracklet search is pairwise in detections/frame, and
did not finish in 35 minutes at that density; killed and confirmed via a
standalone check (`detect_sources` alone, one frame) that the counts were
real, not a hang. The target star is extremely bright at this field's
noise level (flux ≈ 87,000, ~1.7 px from its WCS-predicted position at
every threshold tested from 10 to 100), so raising the demo's threshold
to 60 — cutting detections/frame to ~1,260 — loses none of the signal
this run actually needs; a deliberate, field-specific tradeoff, not the
pipeline's own default. A second bug surfaced once linking finished:
`light_curve`'s default `timestamp_key="MJD-OBS"` doesn't exist in this
survey's headers (`search_field` was already correctly called with
`"OBSMJD"`) — a copy-paste omission in the demo script, not a pipeline
bug, fixed by passing the same key.

With both fixed, the run found 22 variable candidates from 638
detections; 2 matched a known VSX variable, including the target itself —
**ASASSN-V J183620.31, recovered at a 1.7" offset from its catalogued
position** — the positive-control criterion this whole exercise was
built to test, and it passed. A Lomb-Scargle fit to the recovered
candidate's own forced-photometry light curve (`light_curve` +
`recover_rotation_period`) found a real, highly significant periodic
signal (false-alarm probability ≈ 0) at 0.5 d, not the catalogued
0.322 d — the declared-in-advance outcome of fitting a period search to
a light curve covering less than a third of that period (almost
certainly an alias, not evidence against the true period), and not a
failure of `find_variable_sources` itself, which is what the crossmatch
recovery above actually validates.

## Validating against real IASC (Pan-STARRS1) campaign data

Every real-data check so far used ZTF. `docs/src/index.md`'s "Using real
IASC campaign data" section had stood as "not attempted" all session —
closed by running `examples/iasc_demo.jl` against 5 real Pan-STARRS1
(PS1) IASC practice sets (2019-08-28/09-04/09-24, 4 exposures each).
`run_pipeline` recovered 26 real, independently-catalogued objects
across the 5 fields via SkyBoT — including a Jupiter Trojan, 2019 NB9 —
but getting a clean run took three real, fixed bugs, found in this order.

**`load_wcs` failed on every one of these real headers.** Every PS1
header raised `"Linear transformation matrix is singular"` from wcslib.
Bisected a real header down to the exact cause (splitting it into halves,
testing each half in isolation, recursing into whichever half still
failed): `CNPIX1`/`CNPIX2` alone — a legacy IRAF/DSS plate-astrometry
keyword pair, present in these headers but with none of that convention's
other required keywords — was enough to reproduce it, even combined with
nothing but `SIMPLE`/`BITPIX`/`NAXIS`. wcslib reads that keyword pair as
the start of a *separate*, implicit DSS-style WCS description, and with
the rest of that convention absent builds an all-zero, degenerate linear
transform for it — a real wcslib parsing quirk, not anything wrong with
the header's own real, complete CTYPE/CRVAL/CRPIX/CDELT WCS, which parses
cleanly on its own. Fixed in `load_wcs`: on exactly this error, retry
after stripping just those two keyword's FITS cards and nothing else —
confirmed sufficient, not guessed. Regression test constructs a
synthetic header (a real WCS plus injected `CNPIX1`/`CNPIX2` cards) since
the real PS1 files can't be committed to the repo.

**`detect_sources` produced enormous numbers of spurious detections on
some frames.** One real frame: 176,165 "detections" at `threshold=8.0`
(field 451 on ZTF, by comparison, has ~130). Root cause: these FITS files
mark invalid/masked pixels using the standard `BLANK` header keyword
(scaled through `BZERO`/`BSCALE` like any other pixel value) rather than
`NaN`, and `FITSIO.jl` does not convert `BLANK` sentinels automatically.
One real frame had 158,443 pixels (2.7% of the image) pegged at exactly
that sentinel value (65535, from `BLANK=32767` + `BZERO=32768`) —
`detect_sources` read that as enormous real flux across a large masked
region and found a spurious "source" seemingly everywhere. Confirmed
directly: replacing those exact pixels with the frame's own valid-region
median (before any detection) dropped the same frame from 176,165 to 176
detections. Handled as a preprocessing step in `examples/iasc_demo.jl`
(`clean_blank_pixels`, writing cleaned copies preserving the original
header/WCS/timestamp exactly) rather than in `src/`, since BLANK-sentinel
handling is a real-FITS-ingestion concern specific to how a given survey
exports data, not something `detect_sources` itself should need to know
about.

**`crossmatch_catalog(...; :skybot)` was too slow, and not resilient, at
real scale.** Unlike `:vsx`/`:simbad` (batched via CDS TAP — see above),
SkyBoT has no batch mode, so `_crossmatch_skybot` queried one candidate
per request, fully sequentially. Fine for a handful of candidates; not
for hundreds to thousands of real tracklets. A full 5-field run of
`examples/iasc_demo.jl` took over two hours and then died outright, deep
into the fourth field's crossmatch (2,619 candidates), to
`"tls write failed: connection is closed"` — an uncaught, unretried
network error with no partial-progress recovery. Fixed two ways in
`_crossmatch_skybot`: concurrent requests (Julia `Task`s + a bounded
`Base.Semaphore`, not extra threads — this is a network-latency-bound
workload, and cooperative concurrency on however many threads Julia
already has is enough) and a single retry on any `HTTP.HTTPError`.
Benchmarked directly against the live SkyBoT service (20 real, identical
requests, repeated to isolate throughput from any one query's own
content): concurrency=8 gave a real, measured 2.6x speedup (12.5s vs
32.8s sequential); concurrency up to 60 ran clean with zero errors,
though gains flattened past ~20-40 (IMCCE's own server-side queueing, not
this code, by then). Settled on 20 — inside the tested-clean range, not
pushed to its edge, matching `_CDS_BATCH_SIZE`'s conservative philosophy.
The full rerun with both fixes completed end to end, no crash, in around
20 minutes total (all 5 fields) — down from a run that hadn't even
finished after two hours.

**`match_radius` was too loose, by a measured, corrected amount.**
`examples/iasc_demo.jl`'s `match_radius` was first converted from
`real_data_demo.jl`'s ZTF value to preserve the same ~10" angular
tolerance — a reasonable-looking choice that turned out to be looser than
PS1's own real astrometric precision. On the densest of the 5 fields,
this produced 10,422 tracklets from only ~500 detections/frame — almost
certainly distinct real stars within 10" of each other across frames
getting cross-linked into spurious tracklets, not 10,422 real moving
objects. The known SkyBoT objects were still correctly recovered in
every field regardless, but rather than guess at a tighter value, PS1's
own headers report the real number needed: `PERROR`, the astrometric
solution's per-star positional RMS residual, measured at 0.20-0.23"
across the fields checked here — not something assumed, read directly
from real data. Retuned `match_radius` to 2" (~10x `PERROR`, a
comfortable margin for real motion and centroiding noise, not the bare
residual) and reran all 5 fields: the same 9 distinct known objects were
recovered in every field (confirmed by name, not just by count — nothing
dropped out), while total tracklets across all 5 fields fell from 16,158
to 4,960 (-69%). The reduction is concentrated exactly where predicted:
the densest field (XY42_p11) went from 10,422 to 3,478; the two
previously "26 real objects" and "13/2619" style counts were actually
counting duplicate tracklet-rows around the same handful of real
objects, not 26 distinct discoveries — a reporting correction as much as
a code fix, worth noting since the inflated number was reported once,
here, before the retune caught it.

## Tuning `find_variable_sources`'s systematic error floor past the first value that worked

The 1% systematic error floor (see above) was the first value tried,
chosen because it's a standard number in forced-photometry pipelines, not
because it was shown to be optimal. With a real positive control now in
hand (ASASSN-V J183620.31, recovered by `find_variable_sources` earlier
this session), the floor could finally be checked from *both* sides of
the tradeoff at once, not just the false-positive side: sweeping
`systematic_error_fraction` against the same 152 real, matched,
high-S/N stationary stars from ZTF field 451 (a slightly different count
than the "119" quoted earlier — this sweep additionally required
`min_frames` at the full 5-frame count and `normalize=true` together,
narrowing the matched set), the `chi2_threshold=10.0` false-positive rate
dropped from 3.3% at a 1% floor to 2.0% at 2%, 1.3% at 3%, and 0% at 5%.
Naively, that argues for as high a floor as possible — but a floor this
large also suppresses *real* variability, and that side had never been
checked. Running the same sweep against ASASSN-V J183620.31's own real
forced-photometry light curve: reduced chi2 falls from 123 (at 1%) to 31
(at 2%) to 13.8 (at 3%) to 5.0 (at 5%) — the last of which drops *below*
`chi2_threshold=10.0`, meaning a 5% floor would have made this exact,
real, independently-confirmed variable star invisible to
`find_variable_sources`. Settled on 2%: comfortably above 1%'s
false-positive rate, while leaving the real variable's signal a 3x margin
over threshold — the largest floor checked that doesn't cost real
detections, not the smallest false-positive rate achievable.

## Revisiting the two intentional-design "limitations"

Two items in `docs/src/index.md`'s Known Limitations were design
decisions, not bugs — but "intentional" isn't the same as "as good as it
can be." Asked directly whether either could be genuinely improved
without abandoning the design choice behind it.

**`estimate_psf`'s empirical-vs-fallback split was a hard binary — a
field either had stars passing the isolation/saturation filter at exactly
the given `min_separation`, or it fell all the way back to an analytic
Moffat fit.** But a moderately (not severely) crowded field might have
real, usable stars at a *slightly* tighter isolation radius — the current
code was throwing that away and jumping straight to an approximation
instead of trying harder for the real thing. Added `relaxation_attempts`
(default 2): on an empty stamp list, halve `min_separation` and retry,
up to that many times, before falling back. A real empirical PSF from
fewer, closer stars still beats a parametric approximation, which is the
whole reason `estimate_psf` exists over just always using
`fit_moffat_psf`. Regression-tested with two synthetic Gaussian sources
25 px apart (fails the default `min_separation=40`, but 25 ≥ 20, the
first relaxed attempt) — recovers the real empirical PSF (correct FWHM)
instead of falling back, while `relaxation_attempts=0` on the same data
still fails cleanly, proving the relaxation is what does it. Writing that
test surfaced a second thing worth knowing: `fit_moffat_psf`'s own stamp
extraction (`stamp_size=25`, so a ±12 px half-window) doesn't check for
neighbor contamination the way `estimate_psf`'s isolation filter does —
stars closer than ~24 px apart corrupt each other's Moffat fit (tested
directly: clustering the existing fallback test's synthetic stars into a
6 px box made every fit fail to converge, where the original ~25 px
spacing fits cleanly) — not fixed here, since `fit_moffat_psf`'s own
docstring already documents that it deliberately skips the isolation
check as a defensible trade for a *fit* (a bad stamp shows up as a poor
residual, in principle), but the synthetic evidence says that trade has
a real limit worth knowing about.

**`zogy_subtract`'s `V_ast` opt-in was silent — a direct caller who
simply didn't pass `n_sources`/`r_sources` got `V_ast = 0` with no
signal anything was skipped**, unlike `run_pipeline`, which always
supplies them. The design itself (opt-in at the low-level function,
mandatory at the high-level one) is still the right call — computing
`n_sources`/`r_sources` costs two extra `detect_sources` passes, real
cost a direct caller might legitimately not want. What was missing was
visibility: added a `@warn` (once per session, via `maxlog=1`, so it
doesn't spam a caller who's already made an informed choice) when both
are left `nothing`. Regression test confirms it fires when omitted and
stays silent when both are supplied.
