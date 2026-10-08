# Validating against real IASC (Pan-STARRS1) campaign data

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since it's a self-contained validation story
against a different survey (Pan-STARRS1) than the ZTF narrative there.

Every real-data check so far used ZTF. `docs/src/index.md`'s "Using real
IASC campaign data" section had stood as "not attempted" all session —
closed by running `examples/iasc_demo.jl` against 5 real Pan-STARRS1
(PS1) IASC practice sets (2019-08-28/09-04/09-24, 4 exposures each).
`run_pipeline` was reported to recover 9 real, independently-catalogued
objects across the 5 fields via SkyBoT — including a Jupiter Trojan,
2019 NB9 — and getting a clean run took four real, fixed issues, found
in this order. (That 9 was later found to be 2: see
[Recovered objects, re-checked](@ref) below.)

## `load_wcs` failed on every one of these real headers

Every PS1 header raised `"Linear transformation matrix is singular"`
from wcslib. Bisected a real header down to the exact cause (splitting
it into halves, testing each half in isolation, recursing into whichever
half still failed): `CNPIX1`/`CNPIX2` alone — a legacy IRAF/DSS
plate-astrometry keyword pair, present in these headers but with none of
that convention's other required keywords — was enough to reproduce it,
even combined with nothing but `SIMPLE`/`BITPIX`/`NAXIS`. wcslib reads
that keyword pair as the start of a *separate*, implicit DSS-style WCS
description, and with the rest of that convention absent builds an
all-zero, degenerate linear transform for it — a real wcslib parsing
quirk, not anything wrong with the header's own real, complete
CTYPE/CRVAL/CRPIX/CDELT WCS, which parses cleanly on its own. Fixed in
`load_wcs`: on exactly this error, retry after stripping just those two
keyword's FITS cards and nothing else — confirmed sufficient, not
guessed. Regression test constructs a synthetic header (a real WCS plus
injected `CNPIX1`/`CNPIX2` cards) since the real PS1 files can't be
committed to the repo.

## `detect_sources` produced enormous numbers of spurious detections on some frames

One real frame: 176,165 "detections" at `threshold=8.0` (field 451 on
ZTF, by comparison, has ~130). Root cause: these FITS files mark
invalid/masked pixels using the standard `BLANK` header keyword (scaled
through `BZERO`/`BSCALE` like any other pixel value) rather than `NaN`,
and `FITSIO.jl` does not convert `BLANK` sentinels automatically. One
real frame had 158,443 pixels (2.7% of the image) pegged at exactly that
sentinel value (65535, from `BLANK=32767` + `BZERO=32768`) —
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

## `crossmatch_catalog(...; :skybot)` was too slow, and not resilient, at real scale

Unlike `:vsx`/`:simbad` (batched via CDS TAP — see the
[Investigation Log](investigation-log.md)), SkyBoT has no batch mode, so
`_crossmatch_skybot` queried one candidate per request, fully
sequentially. Fine for a handful of candidates; not for hundreds to
thousands of real tracklets. A full 5-field run of `examples/iasc_demo.jl`
took over two hours and then died outright, deep into the fourth field's
crossmatch (2,619 candidates), to `"tls write failed: connection is
closed"` — an uncaught, unretried network error with no partial-progress
recovery. Fixed two ways in `_crossmatch_skybot`: concurrent requests
(Julia `Task`s + a bounded `Base.Semaphore`, not extra threads — this is
a network-latency-bound workload, and cooperative concurrency on however
many threads Julia already has is enough) and a single retry on any
`HTTP.HTTPError`. Benchmarked directly against the live SkyBoT service
(20 real, identical requests, repeated to isolate throughput from any one
query's own content): concurrency=8 gave a real, measured 2.6x speedup
(12.5s vs 32.8s sequential); concurrency up to 60 ran clean with zero
errors, though gains flattened past ~20-40 (IMCCE's own server-side
queueing, not this code, by then). Settled on 20 — inside the
tested-clean range, not pushed to its edge, matching `_CDS_BATCH_SIZE`'s
conservative philosophy. The full rerun with both fixes completed end to
end, no crash, in around 20 minutes total (all 5 fields) — down from a
run that hadn't even finished after two hours.

## `match_radius` was too loose, by a measured, corrected amount

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
solution's per-star positional RMS residual, measured at 0.20-0.23
*pixels* (~0.06") across the fields checked here — not something
assumed, read directly from real data. (This page first reported it as
0.20-0.23" — a misreading: the header's own comment says "(pixels)",
and the separate `CERROR` keyword gives ~0.06" in arcsec, found while
validating the first 2025 set below.) Retuned `match_radius` to 2"
(~35x `PERROR`, a comfortable margin for real motion and centroiding noise, not the bare
residual) and reran all 5 fields: the same 9 distinct known objects were
recovered in every field (confirmed by name, not just by count — nothing
dropped out), while total tracklets across all 5 fields fell from 16,158
to 4,960 (-69%). The reduction is concentrated exactly where predicted:
the densest field (`XY42_p11`) went from 10,422 to 3,478; the two
previously "26 real objects" and "13/2619" style counts were actually
counting duplicate tracklet-rows around the same handful of real
objects, not 26 distinct discoveries — a reporting correction as much as
a code fix, worth noting since the inflated number was reported once,
here, before the retune caught it.

![Tracklet counts per field before and after retuning match_radius to PS1's real astrometric precision](assets/iasc-match-radius-retuning.png)

## Recovered objects (as first reported)

| Field | Tracklets (retuned) | Known objects recovered |
|:--|--:|:--|
| `XY14_p10` | 52 | — |
| `XY15_p01` | 132 | 2014 HO19 |
| `XY25_p10` | 567 | 4311 T-1, 2009 SG135, 2015 XJ232, 2019 PK5 |
| `XY26_p01` | 731 | 2001 SH320, 2011 SH185, 2019 NB9 (Jupiter Trojan) |
| `XY42_p11` | 3,478 | 2008 FA111 |

## The first set from a live campaign's practice round (2025 data)

The 2019 sets above are IASC's public "Practice Image Sets". In October
2026 the International Asteroid Search Campaign opened with a new
practice set, `ps1-NewPractice_3` (PS1 field `XY54_p10`, 4 × 45 s
exposures on 2025-09-16, ~48 min apart), which every team must measure
in Astrometrica and have checked before receiving real images. It was
measured in Astrometrica (running under Wine) independently of this
pipeline, giving an external ground truth for the first time — not just
SkyBoT's known objects:

- **NHU0001**, an object not in SkyBoT within 2', G ≈ 21, moving ~34"/h
  in a straight line at constant brightness — reported to IASC as the
  team's measurement.
- **2018 LT**, a known Mars-crosser at G ≈ 20. It first looked ~15"
  off its SkyBoT ephemeris — but that query was *geocentric*; this
  close-approaching object's parallax accounts for it, and queried for
  the observatory (`-loc F51`) SkyBoT agrees to <3". Astrometrica's own
  known-object box for it was also displaced (~60 px), so it had to be
  found by eye. **Always query SkyBoT topocentrically** for this kind of
  check: `crossmatch_catalog(rows, :skybot; radius, observatory="F51")`.

`examples/iasc_demo.jl` as it stood found **neither**: 117 tracklets,
0 known objects, and every position ~7" off. Each failure had a
separate cause, each confirmed directly before fixing.

### The header WCS was off by 6.8"

Matched against 38 Gaia stars, PS1's header solution was off by a
near-constant ~26.6 px (6.8") across the whole chip — although its own
`CERROR` keyword claims 0.06". Its `PCA*` polynomial distortion keywords
were ruled out as the explanation: evaluated at the chip's corners they
move a position by <0.5 px. Fixed by
[`refine_wcs`](@ref): detect stars in the frame, find the header's error
as one global offset by voting over all detection/catalog pairs, match,
and fit a linear TAN solution with sigma clipping —
Astrometrica's own "Data Reduction" approach. Against Gaia DR3
([`gaia_reference_stars`](@ref)), all four frames refined to a 0.07-0.08"
residual RMS, and NHU0001's positions agreed with Astrometrica's
independently Gaia-calibrated report to 0.08-0.18" (from 7.2").
`run_pipeline(...; refine_astrometry=true)` applies it per frame.

### Linking compared pixel positions across dithered frames

NHU0001 was detected in all four frames — at 5σ even on the raw
frames — and never linked. The raw (no-`reference`) path linked each
frame's *own* pixel positions, but PS1 dithers between exposures (up to
(6, 9) px here), bending a straight sky track by more than the 2"
(7.8 px) `match_radius`. Raw-path detections are now re-expressed in the
first frame's pixel grid before linking (`_to_common_grid`), so the
first frame's WCS applies to every row. This changes raw-path tracklet
counts — static stars now line up perfectly and link as zero-motion
"tracklets" (162 instead of 52 on `XY14_p10`), which is what `min_speed`
is for.

### Gaps between detector cells are filled with a constant, not flagged

~14% of every frame is the gaps between PS1's cells, filled with one
constant (160 ADU), marked neither by `BLANK` nor `NaN`.
[`fill_value_mask`](@ref) masks any constant-valued plateau (real sky
never repeats a value over a 3x3 block) for
[`detect_sources`](@ref)'s new `mask` keyword — survey-agnostic, no fill
value assumed.

### Single-frame detection at 8σ can't reach G ≈ 21

NHU0001 peaks at 7-9σ in each frame. Lowering `threshold` on raw frames
drowns in stars, so [`stack_difference`](@ref) adds the numerical
equivalent of a human's blink: each frame minus the median of the whole
(aligned, flux-scaled) sequence, with each bright star's core *and* halo
masked in proportion to its size — one saturated star's halo stayed 3σ
above sky out to ~40 px, and a fixed margin left its residuals in 14 of
16 surviving tracklets. Tracklets are then filtered by IASC's own "true
signature" tests (straight line, constant speed, ≤1 mag variation:
`min_speed`, `max_flux_ratio`).

### Network failures that used to stop or corrupt a run

Validating all this against the live services surfaced four more real
failure modes, all fixed:

- SkyBoT intermittently answered HTTP 200 with `# Flag: -1` and a
  crashed process's backtrace, which parsed as "no known objects" —
  every known object silently became a discovery. Now a
  [`SkyBoTServiceError`](@ref), retried, never an empty match.
- A catalog request hung for 30 minutes at 0% CPU with no timeout; all
  catalog requests now have an idle timeout.
- VizieR's TAP answered 503 for a sustained stretch, and later a TLS
  handshake broke mid-write; catalog queries retry with backoff on any
  service-side failure, and Gaia stars come from the ESA Gaia archive
  with VizieR as fallback.
- `load_wcs`'s `CNPIX1`/`CNPIX2` workaround only caught wcslib's
  "singular" message, but the same header failed "Invalid parameter
  value" inside a longer run — the implicit solution's fields are
  uninitialized memory, so the message varies; the retry now keys on
  the keywords being present, not on the message.

### More junk the difference images surfaced on the 2019 sets

Running the new path over the five 2019 sets as a regression check
turned up three more kinds of spurious detection, each measured before
being filtered:

- **Noisier detector cells.** One frame-wide noise figure let a couple
  of noisier cells produce ≈3,000 of a frame's ≈4,000 6σ detections
  (`XY15_p01`). [`stack_difference`](@ref) now divides each difference by
  a *local* (128 px box) noise estimate — floored at half the frame's
  own, after one box of near-constant pixels drove a frame's normalized
  values to about 1e15 and buried a 30σ asteroid (`XY25_p10`).
- **Cosmic rays and hot pixels.** 94-98% of the remaining detections
  were single-pixel spikes (neighbour/peak ratio ≈ 0, against 0.46-0.82
  for every real object measured). `detect_sources`'s new `sharpness`
  column and `min_sharpness=0.3` drop them.
- **Satellite trails.** Two broad trails in `XY25_p10` fragmented into
  rows of peaks that linked into about 240 3-frame tracklets; detections with
  two or more neighbours within 15 px are now dropped on the difference
  path.

## Recovered objects, re-checked

The table above counted a tracklet as "recovering" a known object if its
*first* point lay within 15" of SkyBoT's position. Re-checked properly —
SkyBoT queried topocentrically (`-loc F51`) at every frame's own epoch,
and a single tracklet required to sit near the object's predicted
position in at least 3 frames — most of those were static stars: the
original configuration truly recovered **2** objects in these five
sets, not 9. Some catalogued objects aren't recoverable from these
frames at all: fainter than about V 21.5, off the chip after all, or
(4311 T-1, V 19.3) in a cell gap in one frame and directly over a static
star in the next.

All six local sets, original configuration against the one now in
`examples/iasc_demo.jl` (3" match radius for the refined positions,
12" for the original's header-WCS ones; `XY54_p10`'s unknown object
checked against its Astrometrica measurement). Every recovery in the
new columns was confirmed to *move with* its ephemeris: within 1.3" of
it in every frame, at the ephemeris's own rate to 0.1"/h.

| Field | Real objects | Original (raw, threshold 8) | New, `min_frames=4` | New, `min_frames=3` (demo) |
|:--|:--|:--|:--|:--|
| `XY14_p10` | — | 162 tracklets | 0 | 0 |
| `XY15_p01` | — (2 known, V ≈ 21.5, not seen) | 234 | 0 | 0 |
| `XY25_p10` | 2009 SG135, 2015 XJ232, 2019 PK5 | 583: SG135² | 5: all 3 | 116: all 3 |
| `XY26_p01` | 2019 NB9 (Jupiter Trojan) | 743: ✓ | 1: ✓ | 13: ✓ |
| `XY42_p11` | — (3 known, V ≥ 20.6, not seen) | 3,528 | 1 | 25 |
| `XY54_p10` (2025) | 2018 LT, NHU0001 | 298: 2018 LT only¹ | 1: NHU0001 only | 3: both |
| **Total** | **6** | **5,548 tracklets, 3/6** | **8, 5/6** | **157, 6/6** |

¹ With the dithered-linking fix already applied; before it, the original
configuration found neither (117 tracklets).

These columns are the configuration as of that check; it has since been
improved to 7 real objects with 95 tracklets — see
[Completeness, measured by injecting synthetic asteroids](@ref).
² The evaluation also credits it with 4311 T-1, but that tracklet is the
static star the asteroid passes over (its 12" radius, sized for the
header WCS's error, lets a star near the middle of the track count).

2015 XJ232 (V 21.4) and 2019 PK5 (V 21.3) are new: the original
configuration's "recoveries" of them were static stars, and a first run
of this re-check missed them too — its scoring loop was written
`for name in names, id in ids ... break`, and in Julia that `break`
leaves *both* loops, so it stopped looking after the first object found
per set.

## Reports, epochs, photometry and completeness

### The pipeline's own MPC report, checked against Astrometrica's

`examples/iasc_demo.jl` now writes an 80-column MPC report for every
candidate that isn't a confirmed known object, with Gaia-calibrated
magnitudes ([`candidate_magnitudes`](@ref): each frame's WCS refined and
zero point fitted against Gaia DR3 — 23-27 stars per `XY54_p10` frame,
0.02-0.05 mag scatter), and compares it line by line with an
Astrometrica report of the same set when one is present. The first
comparison failed outright — no pipeline measurement at any of
Astrometrica's four epochs — and found two epoch errors:

- PS1's `MJD-OBS` is the shutter time in **TAI**: 37 s after
  `DATE-OBS`, "UTC start of exposure", in every IASC header checked
  (2019 and 2025), although the same header says `TIMESYS = 'UTC'`.
- Both PS1 and ZTF stamp the exposure's *start*; MPC reports need its
  *middle* (22.5 s later on PS1).

[`frame_epoch`](@ref) now reads `DATE-OBS` when present and adds half
of `EXPTIME`. Re-run, the pipeline's report for NHU0001 matched
Astrometrica's at all four epochs, **0.13-0.18" apart** in position;
magnitudes differed by 0.1-0.9 mag (21.4/21.4/21.7/21.8 against
21.2/21.0/21.6/20.9, at S/N 4-8 per frame). `digest2`, now built and on
`PATH` (see [MPC digest2 Scoring](@ref)), gives it NEO score 2: a
main-belt object, consistent with its ~34"/h motion.

### Completeness, measured by injecting synthetic asteroids

Six real moving objects can't trace a completeness curve, so
[`inject_movers`](@ref) adds synthetic ones to real frames — true sky
positions through each frame's refined WCS, true brightness through its
own zero point, the frame's measured PSF, Poisson noise, and the trail a
moving object leaves during the exposure — and
[`injection_recovery`](@ref) scores a run against them
(`examples/injection_test.jl`). On `XY54_p10` (180 main-belt-like
objects, 18-23 mag, 5-100"/h; 80 fast ones, 100-1000"/h):

| G mag | 18-22 | 22.0-22.5 | 22.5-23 |
|:--|:--|:--|:--|
| recovered | ~80% (plateau) | 16% | 0% |

The plateau's missing ~20% is the field itself: detector-cell gaps cover
~14% of every frame, and bright stars' cores and halos are masked. The
50% limit is G ≈ 22.1.

Injection found two real weaknesses, both fixed:

- **Fast objects' trails broke into several peaks**, which the cluster
  filter removed as junk: 300-700"/h objects were recovered 1 time in
  20. Compact groups are now merged into one detection at the trail's
  centre; only groups longer than a trail (satellite streaks) are
  dropped ([`_merge_clusters`](@ref AsteroidPipeline._merge_clusters)).
- **`link_candidates` only seeded from frames 1 and 2**, even with
  `min_frames=3`: an object missing from either — in a gap, over a star,
  or out of the field — could never be linked. Any frame pair now seeds
  when a frame may be missing.

Fast objects went from 8/80 to 22/80 recovered (44% up to 400"/h);
above ~500"/h most simply cross the 10' chip within one or two frames.
The extra seeds also produced more chance 3-frame alignments, so
partial tracklets must now lie on a straight line to within 1.0 px RMS
(`max_residual`; real objects measured 0.12-0.31 px, spurious ones a
median 2.6 px). Complete tracklets are exempt: faint real objects
reached 1.40 px.

Re-checked against SkyBoT ephemerides on all six sets, the configuration
in `examples/iasc_demo.jl` now recovers **7 real moving objects** — the
6 above plus 2008 FA111 (V 20.6, in `XY42_p11`), found for the first
time — with **95** tracklets to vet across the six sets.

Known limits, measured but not solved: slow objects (5-10"/h, 50%)
partly subtract themselves, since they overlap their own positions in
the median of a ~50 min set — only a deeper reference from other nights
(the ZOGY path) avoids that; and an object passing over a bright star is
masked out in that frame.

Sub-pixel alignment for [`stack_difference`](@ref) (`subpixel_alignment`
in [`run_pipeline`](@ref)) was measured the same way — 180 injected
objects in each of two sets — and is a trade-off, so it is left off:
it halved spurious tracklets on the star-heavy `XY25_p10` (356 tracklets
in all to 182) but recovered 5-9% fewer objects.
