# Design refinements

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since these aren't bugs found on one specific
real dataset the way the other pages here are: two are design decisions
revisited to see if they could be improved without abandoning the choice
behind them, and two are validation/performance investigations that
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

**No further speedup found.** `Reproject.reproject` does expose an
`order` keyword (interpolation order — `0` for nearest-neighbor instead
of the default bilinear), but the segfault's own stack trace shows the
real per-pixel cost is in `pix_to_world!` itself — the WCS coordinate
transform each output pixel needs before any interpolation happens at
all — which `order` has no effect on, so it wasn't pursued: a real
accuracy cost (nearest-neighbor reprojection reintroduces the same kind
of sub-pixel registration error the reference stack exists to average
out) for a speedup that the evidence says wouldn't materialize. The
"tens of minutes" cost is, as far as this investigation could establish,
a real, currently-irreducible property of reprojecting real frames
through `wcslib` one at a time — not something left unoptimized for lack
of trying.
