# Validating against real IASC (Pan-STARRS1) campaign data

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since it's a self-contained validation story
against a different survey (Pan-STARRS1) than the ZTF narrative there.

Every real-data check so far used ZTF. `docs/src/index.md`'s "Using real
IASC campaign data" section had stood as "not attempted" all session —
closed by running `examples/iasc_demo.jl` against 5 real Pan-STARRS1
(PS1) IASC practice sets (2019-08-28/09-04/09-24, 4 exposures each).
`run_pipeline` recovered 9 real, independently-catalogued objects
across the 5 fields via SkyBoT — including a Jupiter Trojan, 2019 NB9 —
but getting a clean run took four real, fixed issues, found in this
order.

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

![Tracklet counts per field before and after retuning match_radius to PS1's real astrometric precision](assets/iasc-match-radius-retuning.png)

## Recovered objects

| Field | Tracklets (retuned) | Known objects recovered |
|:--|--:|:--|
| XY14_p10 | 52 | — |
| XY15_p01 | 132 | 2014 HO19 |
| XY25_p10 | 567 | 4311 T-1, 2009 SG135, 2015 XJ232, 2019 PK5 |
| XY26_p01 | 731 | 2001 SH320, 2011 SH185, 2019 NB9 (Jupiter Trojan) |
| XY42_p11 | 3,478 | 2008 FA111 |
