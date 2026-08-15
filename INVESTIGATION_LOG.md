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
