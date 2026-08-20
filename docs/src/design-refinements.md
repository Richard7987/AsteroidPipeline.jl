# Design refinements

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since these aren't bugs found on one specific
real dataset the way the other pages here are: two are design decisions
revisited to see if they could be improved without abandoning the choice
behind them, and the rest are validation/performance investigations that
don't belong to any single dataset-specific page.

## `estimate_psf`'s empirical-vs-fallback split was a hard binary

A field either had stars passing the isolation/saturation filter at
exactly the given `min_separation`, or it fell all the way back to an
analytic Moffat fit. But a moderately (not severely) crowded field might
have real, usable stars at a *slightly* tighter isolation radius — the
current code was throwing that away and jumping straight to an
approximation instead of trying harder for the real thing. Added
`relaxation_attempts` (default 2): on an empty stamp list, halve
`min_separation` and retry, up to that many times, before falling back.
A real empirical PSF from fewer, closer stars still beats a parametric
approximation, which is the whole reason `estimate_psf` exists over just
always using `fit_moffat_psf`. Regression-tested with two synthetic
Gaussian sources 25 px apart (fails the default `min_separation=40`, but
25 ≥ 20, the first relaxed attempt) — recovers the real empirical PSF
(correct FWHM) instead of falling back, while `relaxation_attempts=0` on
the same data still fails cleanly, proving the relaxation is what does
it. Writing that test surfaced a second thing worth knowing:
`fit_moffat_psf`'s own stamp extraction (`stamp_size=25`, so a ±12 px
half-window) doesn't check for neighbor contamination the way
`estimate_psf`'s isolation filter does — stars closer than ~24 px apart
corrupt each other's Moffat fit (tested directly: clustering the
existing fallback test's synthetic stars into a 6 px box made every fit
fail to converge, where the original ~25 px spacing fits cleanly) — not
fixed here, since `fit_moffat_psf`'s own docstring already documents
that it deliberately skips the isolation check as a defensible trade for
a *fit* (a bad stamp shows up as a poor residual, in principle), but the
synthetic evidence says that trade has a real limit worth knowing about.

## `zogy_subtract`'s `V_ast` opt-in was silent

A direct caller who simply didn't pass `n_sources`/`r_sources` got
`V_ast = 0` with no signal anything was skipped, unlike `run_pipeline`,
which always supplies them. The design itself (opt-in at the low-level
function, mandatory at the high-level one) is still the right call —
computing `n_sources`/`r_sources` costs two extra `detect_sources`
passes, real cost a direct caller might legitimately not want. What was
missing was visibility: added a `@warn` (once per session, via
`maxlog=1`, so it doesn't spam a caller who's already made an informed
choice) when both are left `nothing`. Regression test confirms it fires
when omitted and stays silent when both are supplied.

## `photometric_outlier_threshold` was untested, but not for lack of a real anomaly to test it against

`run_pipeline`'s raw-path `photometric_outlier_threshold` docstring used
to say the check had "not been validated against a real,
independently-confirmed anomaly" — technically true, but the *reason*
given (every frame's raw-path photometric scale on the field-451 dataset
stayed within ~3-13% of the median, "so this signal did not, and on this
dataset could not have, flagged that frame") turned out to conflate two
different things once actually measured, one frame at a time rather than
summarized as a range. The one real, confirmed anomaly available is a
specific frame — the likely passing cloud `quality_max_std` catches, via
its elevated `S_corr` standard deviation (~2.0 vs. ~1.1-1.2) and ~232
excess bright residuals. Measuring `photometric_scale` for *that exact
frame* directly (not the dataset's full spread): 0.64% deviation from
the field median — an order of magnitude under even a strict threshold,
let alone the default 20%.

That's not "the check is untested" — it's a real, measured negative
result with a real explanation: `quality_max_std` and
`photometric_outlier_threshold` are sensitive to different failure
modes, not two safety nets for the same one. A cloud during ZOGY
differencing inflates per-pixel residual noise in `S_corr` directly —
exactly what `quality_max_std` watches — without necessarily causing a
*uniform* per-frame flux-scale shift across every star, which is the
only thing `photometric_scale`'s ensemble-ratio approach can see. The
check remains genuinely unvalidated for the failure mode it's actually
meant to catch (real transparency loss, guiding/focus problems — a
uniform brightness shift), since no real example of *that* has turned up
in this project's data yet — but "never tested against any real
anomaly" was no longer an accurate way to describe it, and now isn't.

## `build_reference`'s real bottleneck, a real 2x win, and a real crash found and reverted

`real_data_demo.jl`'s documented "tens of minutes" reference-build time
had never been profiled — just described. Measured directly, on the
real 30-frame field-451 reference set:

**The per-pixel median-combine loop had a real, fixable cost.** The
original code built a fresh `[stack[ci, k] for k in ... if valid[ci,k]]`
array comprehension for *every pixel* — one heap allocation each, over
a million times for a real frame. Replacing it with a single reusable
buffer, filled in place per pixel instead of reallocated, gave a real
2x speedup with byte-identical output (verified directly, not assumed
from the allocation count alone), measured on the same real 30-frame
field-451 combine step:

| | time | allocations | GC time |
|:--|--:|--:|--:|
| Original (fresh array per pixel) | 1.39 s | 9.46 M | 31% |
| Reusable buffer | 0.68 s | 4.20 M | 5% |

**But the combine step was never the real bottleneck.** The same
benchmark that measured the 2x win also measured `Reproject.reproject`
itself: ~24s per frame, ~two orders of magnitude more than the *entire*
combine step for all 30 frames combined (~1s). The real "tens of
minutes" cost is almost entirely reprojection, not combination.

**Parallelizing reprojection across frames — obviously safe on paper,
genuinely unsafe in practice.** Each frame's reprojection is
independent of every other frame's, and `Reproject.jl`'s own source has
no shared mutable state (checked directly). Wrapped the per-frame loop
in `Threads.@threads` on that basis — and a real multi-threaded run on
real data segfaulted inside `WCS.jl`'s `pix_to_world!`, which wraps
`wcslib` (a C library) via `ccall`. Checking a Julia package's own
source for thread-safety isn't enough when it calls into a C library:
the transitive dependency needs the same scrutiny, and `wcslib`
apparently doesn't tolerate concurrent calls. Reverted immediately back
to a sequential loop; kept as a documented, real finding (in
`build_reference`'s own docstring) so the same "obviously parallelizable"
mistake isn't attempted again the same way.

`Reproject.reproject` does expose an `order` keyword (interpolation
order — `0` for nearest-neighbor instead of the default bilinear), but
the segfault's own stack trace showed the real per-pixel cost is in
`pix_to_world!` itself — the WCS coordinate transform each output pixel
needs before any interpolation happens at all — which `order` has no
effect on, so it wasn't pursued: a real accuracy cost (nearest-neighbor
reprojection reintroduces the same kind of sub-pixel registration error
the reference stack exists to average out) for a speedup the evidence
said wouldn't materialize.

**Processes succeed where threads crashed — with a real, different
pitfall of their own.** `Distributed.jl` sidesteps the specific hazard
above: each worker process has its own independent `wcslib` state, so
concurrent calls from different processes aren't the same shared-state
problem concurrent calls from different *threads* are. But the obvious
first attempt — reproject each frame inside a `pmap` closure that
receives the frame's already-built `WCSTransform` — segfaulted too, for
a different reason: `WCSTransform` holds pointers into `wcslib`-allocated
C memory that's only valid in the process that created it, and Julia's
generic `Serialization` doesn't reconstruct that state on the receiving
worker — confirmed via a real crash inside `WCS.jl`'s `getproperty`/
`convert_string`, deserializing a `WCSTransform` sent from the main
process. The fix: never let a `WCSTransform` cross the wire at all.
`WCS.jl` provides `to_header`/`from_header`, an exact round trip through
a plain FITS header string (just a `String`, no pointers, serializes
safely) — send that instead, and have each worker rebuild its own local
`WCSTransform` via [`load_wcs`](@ref) before calling
`Reproject.reproject`. Confirmed end to end on `build_reference` itself
(not just the technique standalone), on the same real 30-frame field-451
reference set, 8 worker processes, sequential and distributed measured
back to back with identical output: no crash, 295.92s vs. 951.31s
sequential — a real **3.21x**, not full 8x core-count scaling, since
each worker still re-parses its own WCS header per call and `pmap`'s own
scheduling/serialization isn't free. Landed as
`build_reference`'s `workers` keyword (opt-in; the caller supplies
already-running worker processes, since spawning and managing a process
pool is an environment concern, not something a data-processing function
should own as a side effect).

The "tens of minutes" sequential cost is, as far as this investigation
could establish, a real, currently-irreducible property of reprojecting
one frame through `wcslib` at a time — but it is no longer irreducible
*in total*: spreading that same per-frame cost across independent
processes is a real, measured 3.21x win, not something left unexplored
for lack of trying.

## Cross-validating `find_variable_sources`'s systematic error floor against two more real, independent variable stars

The 2% `systematic_error_fraction` default (see the Investigation Log's
own floor-tuning story) rested on a single real field: false positives
measured against 152 real ZTF field-451 stars, sensitivity checked
against one real confirmed variable, ASASSN-V J183620.31, from a
*different* field entirely. One confirmed variable is a thin evidence
base for a default that trades away real sensitivity — asked directly
whether more real variables could be found to check it against, instead
of trusting one field to generalize.

**Finding two more real, densely-sampled confirmed variables.** Queried
VizieR's VSX table (the same real TAP service `crossmatch_catalog(...;
:vsx)` already uses) for bright (mag 12-16), short-period (0.2-0.6 d),
real eclipsing binaries, then checked each candidate's actual ZTF
exposure history via IRSA's metadata API for real same-night,
same-quadrant density — most real candidates only had 5-7 exposures on
their best night (ordinary ZTF cadence; the earlier 144-exposure night
was a special high-cadence campaign, not typical). Out of 150 real
candidates checked this way, two stood out sharing the same field and
night: **V1012 Mon** (333 exposures) and **ASASSN-V
J072906.85-090518.2** (332 exposures), both field 360, night
2019-01-08 — denser real coverage than the original 144-exposure case,
and both in a real, unusually dense stellar field (Monoceros straddles
the galactic plane): ~30,000 raw detections/frame at the pipeline's
default `threshold=5`, needing `threshold=100` to bring that down to a
workable ~200-580/frame, still denser than field 451's own ~130.

**Both real targets recovered cleanly.** Run through the same
`search_field`-style matching field 451 used: V1012 Mon recovered 0.63"
from its VSX position, ASASSN-V J072906.85 at 1.55" — both well inside
typical match tolerances — with 197 and 182 real, full-coverage matched
stationary stars respectively (more than field 451's 152, in each
field individually).

**Repeating the exact same floor sweep three times over, not once.**

| Floor | Field 451 FP% / var χ² | V1012 Mon FP% / var χ² | ASASSN-072906 FP% / var χ² |
|:--|--:|--:|--:|
| 1% | 3.3% / 123.4 | 6.12% / 542.3 | 2.76% / 151.2 |
| 2% (old default) | 2.0% / 31.0 | 2.55% / 137.4 | 0.55% / 39.2 |
| **3% (new default)** | **1.3% / 13.8** | **2.04% / 61.2** | **0.55% / 17.5** |
| 4% | (not measured) | 1.53% / 34.5 | 0.55% / **9.9** ← below threshold |
| 5% | 0.0% / **5.0** ← below threshold | 1.02% / 22.1 | 0.55% / 6.3 ← below threshold |

At exactly 3% — a floor already measured on field 451 in the original
sweep, so no extrapolation needed there — **all three** real confirmed
variables stay above `chi2_threshold=10.0` (margins 1.38x, 6.1x, 1.75x),
while the false-positive rate is the same or better than 2% in every
single field. At 4%, ASASSN-V J072906.85's own signal already drops
below threshold; at 5%, two of the three do. A literal 0% false-positive
rate is real and reachable (field 451 and, separately, ASASSN-072906
both hit it at 5%) — but as a byproduct of losing real sensitivity, not
as a genuine improvement, exactly the tradeoff the original single-field
sweep already warned about, now confirmed rather than assumed on two
more independent real datasets. `systematic_error_fraction`'s default
moved from 2% to 3% on this evidence; 0% remains achievable only by
accepting that cost, which is why it isn't the target.

![False-positive rate and each field's own real confirmed variable's reduced chi², swept across systematic_error_fraction, for three independent real ZTF fields](assets/systematic-error-floor-sweep.png)

## Looking for another `build_reference`-sized win — profiled, and mostly not found

Asked directly whether anything else in the pipeline was worth
parallelizing or otherwise speeding up, the way `build_reference`'s
reprojection loop was. Profiled the real stages of `_detect_all_frames`,
`find_variable_sources`, `photometric_scale`, and `link_candidates`
directly (`@time`, not guessed) against the V1012 Mon dataset above —
useful for this specifically because its field is unusually dense
(~200-580 detections/frame even at `threshold=100`, denser than field
451's own ~130), making it a real stress case rather than a typical one.

**Most stages are already fast.** `find_variable_sources` (0.20s) and
`photometric_scale` (0.08s) across all 34 frames were not remotely close
to being a bottleneck, despite an initial worry that their per-star
position-matching loop might scale badly with detection count.
`detect_sources` itself: 0.06s for a single dense frame (523 sources).

**One real, safe fix found: `_detect_all_frames`'s `detections_per_frame`
was an untyped `Vector{Any}`.** `detect_sources` always returns the exact
same concrete `TypedTables.Table` type on every call — confirmed
directly (`typeof(...)` compared for both its empty and non-empty return
paths, not assumed from reading the source alone) — so there was no
reason for the container holding those results to be `Any`-typed. Fixed
by typing it explicitly from that same concrete type. Measured before
and after on this same dense dataset: no significant wall-clock change
here (within noise, ~5.4-5.6s either way) — `_detect_all_frames`'s real
cost is FITS I/O and `detect_sources`'s own internal allocations, not
this container's typing. Kept anyway as a real, zero-risk type-stability
correctness improvement (every downstream consumer of
`detections_per_frame` — `find_variable_sources`, `link_candidates`,
`photometric_scale` — now sees a concretely-typed vector instead of
`Any`), not reported as a performance win it didn't measurably produce.

**`link_candidates`'s real cost on this dense field (0.57s, 52.6% GC
time) is a known, already-mitigated design tradeoff, not a hidden
bug.** Its seed-tracklet search is a pairwise loop over every
frame-1/frame-2 detection pair — confirmed directly: 523 × 365 = 190,895
pairs checked here, of which only 658 survive the `max_speed` filter and
go on to the more expensive per-candidate work (a linear
`_closest_detection` scan across every other frame, plus a
least-squares refit that itself allocates three small temporary arrays
per candidate). On a typical field (like field 451's own ~130
detections/frame) the pair count drops by roughly (130/523)² ≈ 6%,
putting this well under 50ms — this 0.57s figure is specific to an
unusually dense field, and `threshold` (already tuned up for exactly
this kind of dense field in both this dataset and the earlier field-487
case) is the pipeline's existing, deliberate answer to it, not something
this profiling pass found reason to change.

No second `build_reference`-sized win turned up. The pipeline's real
remaining cost, at typical real field densities, is dominated by what
was already found and fixed.
