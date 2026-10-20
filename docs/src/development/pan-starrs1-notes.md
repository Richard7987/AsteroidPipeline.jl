# Pan-STARRS1 working notes

The detailed record behind [Validation: Pan-STARRS1](../validation/pan-starrs1.md),
in the order things were found: every failure on real Pan-STARRS1 images, how it
was diagnosed and what fixed it. The images were the public practice sets and
campaign sets of the International Astronomical Search Collaboration (IASC).


The first run, on five 2019 practice sets (four exposures each), reported 9
catalogued objects recovered via SkyBoT — including a Jupiter Trojan,
2019 NB9. That 9 was later found to be 2 (see
[Recovered objects, re-checked](@ref) below). Getting a clean run took four
fixes, found in this order.

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
[Investigation Log](ztf-investigation.md)), SkyBoT has no batch mode, so
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

![Tracklet counts per field before and after retuning match_radius to PS1's real astrometric precision](../assets/iasc-match-radius-retuning.png)

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
`PATH` (see [Follow-up tools](../guides/follow-up.md)), gives it NEO score 2: a
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

## A blind test on four live campaign sets (2026-10-04)

The four sets of the 2026-10-04 delivery (`XY72_p01`, `XY72_p11`,
`XY73_p00`, `XY73_p01`) were run through the pipeline *before* being
vetted by hand in Astrometrica. The vetting found no new object in any of
them, and the pipeline missed nothing real: it linked 9 known asteroids,
plus one SkyBoT didn't place. But each of its "new" candidates was false,
for one of four reasons, and four changes now remove all of them.

### A known object SkyBoT didn't know where to find

One candidate in `XY72_p01` turned out to be **2019 UA185**. SkyBoT, and
JPL's Horizons, placed it far from where it was: both still used an older
orbit, while the MPC's own orbit file, updated the same day (MPEC
2026-T35), put it within ~4". Astrometrica, which reads that file
(`MPCORB.DAT`), identified it at once.

[`crossmatch_catalog`](@ref) now takes `:mpcorb`, with the orbits from
[`read_mpcorb`](@ref): each orbit is propagated (two-body, with
light-time) to the candidate's epoch, as seen from the observatory, whose
position comes from Horizons. A two-body orbit carried months from its
epoch is only good to arcseconds, and for archival data years away from
it, not even that — so `examples/iasc_demo.jl` asks *both* SkyBoT (5")
and `MPCORB.DAT` (15"), and either one identifies an object.

### A hot pixel, linked as a slow mover

Another candidate in `XY72_p01`, a "slow object", was a single hot pixel:
fixed on the detector, it moved across the sky with the telescope's small
dithers between exposures, and the linker saw a mover. Real objects
measured moved ≥32 px on the detector across a set; this one, 0.
`detector_artifact_radius` (2 px in the demo) drops any detection that
recurs at the same *detector* pixel in two or more other frames
([`_drop_detector_fixed`](@ref AsteroidPipeline._drop_detector_fixed)).

### Ghosts and diffraction spikes of bright stars

`XY72_p11` gave a candidate on a bright star's internal reflection (a
ghost), and `XY73_p01` one on a saturated star's diffraction spikes,
which rotate with the telescope from frame to frame and so survive the
difference stack. Both sit on light that a real faint object lacks: the
median of an 8-20 px annulus around them, in the raw frame, stood 1.72σ
(ghost) and 2.44σ (spike) above the sky, against −0.07 to 0.12σ for every
real object measured. `max_background_excess` (1.0σ in the demo) drops
those ([`_drop_on_bright_background`](@ref AsteroidPipeline._drop_on_bright_background)).

### A set that was not worth reducing

`XY73_p00` gave 26 tracklets, all junk from defective detector cells: in
its frames, 10σ detections outnumber the Gaia stars of the field
137-1126 to 1, where on every other set, good or merely difficult, the
ratio is 0.6-5.4. The astrometric refinement now warns when it exceeds
30, and [`set_quality`](@ref) tells a script whether to trust a set at
all; the demo then writes no candidate report for it, and says so.

### Effect

On the four new sets, the demo now writes no false candidate:
`XY72_p01`, `XY72_p11` and `XY73_p01` end with every tracklet identified
as a known object (2019 UA185 through `MPCORB.DAT`), and `XY73_p00` is
reported as unreliable. `XY54_p10`'s NHU0001 is still found.

Re-checked against SkyBoT on all ten sets, the two filters cost no real
object: 14 known objects recovered with or without them, with 125
tracklets instead of 133. On injected objects (`XY54_p10`, as above)
they cost 2 of 180 main-belt objects (110 → 108, both at V 19.5-20.5)
and none of the 80 fast ones (22 → 22).
`set_quality` flags none of the nine other sets.

## Eight sets of 2026-10-12 and 10-14: the frames' astrometry fitted together

The first analysis of eight new sets gave 60 candidates: 31 in `XY75_p01`,
21 in `XY05_p01`, 7 in `XY04_p01`, 1 in `XY04_p11`. Their rates gave them
away: nearly all moved 2.4–6″/h, all in one set the same way, by the same
per-frame offsets — stars, displaced by the frames' astrometry.

### Headers 32″ off

On `XY04_p11` and `XY05_p00` [`refine_wcs`](@ref) failed on every frame
(1–3 matches). Its vote searches offsets up to 60 px; rerun with 300 px,
the offset was (−112, −64) px on one and (111, 65) px on the other, the
same on all four frames, and 10–16 stars matched at 0.16–0.36″. The one
"candidate" in `XY04_p11` was (273964) 2007 KS7, 32″ off its predicted
position — the header's error — and so not identified.

### Frames out of step

On the other sets the refinement succeeded, yet frames disagreed. The
headers were compared with the refitted solutions on the sources detected
in two frames (120–800 of them), each mapped from one frame into the other:
through the headers they lined up to 0.1–0.35 px; through the refitted
solutions, off by up to 9 px. Each frame had been fitted to its own 6–13
stars of ~250, and 6 parameters on so few, mostly in one part of the chip,
extrapolate poorly. The single voted offset was itself the limit: it
differed between frames of one set (12.8–20.8 px on `XY04_p01`) — the
header's error varies across the chip, so one offset matches only where it
holds. Matched again from the first fit (8 px search), 70–90 stars matched
per frame at ~0.1″; a frame stuck at its first few starts again from the
best frame's solution carried through the headers.

Pooling all frames' stars into one fit per frame, through the headers'
relative alignment, then needs a check that the headers really do agree:
the pooled fit is repeated with a free offset for each frame's stars, and
rejected if any offset is significant. A first check — each frame's mean
residual — was too weak: a synthetic header 8 px out of line passed it, the
6-parameter fit absorbing the offset when that frame's stars sat in one
part of the chip.

### The flux test, per point

Aligned correctly, 2019 UA185 (`XY72_p01`) was detected in all four frames
within 0.5–1.1 px of Astrometrica's positions — closer than before (0.7–4.1
px) — and still not reported: fluxes 121, 114, 47 and 95, a ratio of 2.57
against the 2.5 allowed. The tracklet is now kept without its most
discordant point when the rest passes every test, the partial-tracklet
straight-line test included.

### Effect

The eight sets: no candidate, 2007 KS7 identified, residuals 0.10–0.13″.
On 20 sets, 176 → 119 tracklets, 150 → 91 not identified, 24 → 28 catalogued
asteroids identified; 2020 HW140 (V 21.9, by a cell gap) lost. Injection on
`XY54_p10`: 108 → 113 of 180, 22 → 27 of 80. NHU0001 against Astrometrica:
0.12–0.17″ (0.13–0.18″ before).
