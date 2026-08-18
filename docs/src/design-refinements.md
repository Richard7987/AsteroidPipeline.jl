# Revisiting the two intentional-design "limitations"

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since it's a design review, not a bug found on a
specific real dataset the way the other pages here are.

Two items in `docs/src/index.md`'s Known Limitations were design
decisions, not bugs — but "intentional" isn't the same as "as good as it
can be." Asked directly whether either could be genuinely improved
without abandoning the design choice behind it.

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

## GitHub Pages was never actually serving the wiki, despite a working deploy pipeline

The GitHub mirror's docs site returned 404 even though
`.github/workflows/Documenter.yml` reported `success` on every run, and
`deploydocs` printed a clean `Deploying: ✔` with all five criteria
checked. Diagnosed with `gh` (installed and authenticated specifically
for this, not guessed from log snippets) rather than assumed from the
green checkmark: `gh api repos/.../pages` returned a plain 404
("Not Found") — GitHub Pages itself was disabled at the repository
level, Source stuck on "None". The workflow's own job log, read in
full, confirmed the deploy step really had pushed a working site:
`fatal: 'upstream/gh-pages' is not a commit` (the expected first-deploy
message), followed by a clean orphan-branch creation and a real push —
git's own "Create a pull request for 'gh-pages'" hint, which only
appears when a *new* branch is actually accepted by the remote. So the
branch existed with real content; nothing was being served from it
because Pages was off, a setting independent of whether gh-pages exists.
Fixed directly: `gh api repos/.../pages -X POST` with
`source.branch=gh-pages`, confirmed via the API's own `status: "built"`
and a live `curl` returning 200 on multiple pages.

A related, real mistake made and caught in the same investigation: a
first attempt at adding a `.nojekyll` file (see above) touched it into
`docs/build/<n>/` before `deploydocs` ran, reasoning that this local
directory was "the site root." It isn't — `deploydocs` copies that
directory's contents into a nested `dev/` subfolder of the actual
gh-pages branch, generating `versions.js` and a redirect `index.html`
separately at the *real* root, which no local directory corresponds to.
Confirmed by listing the live branch's root contents via `gh api`
(`dev/`, `index.html`, `versions.js` — no `.nojekyll` anywhere) after
the first attempt had already been deployed. Corrected by adding
`.nojekyll` directly to the gh-pages root via the API instead, a
one-time fix rather than build-time logic: `deploydocs`'s own push
mechanism checks out the *existing* branch and only touches specific
known paths (each version's subfolder, `versions.js`, the root
`index.html`) rather than wiping it, so a root file added once persists
across every future automated deploy.

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
`Reproject.reproject`. Confirmed on the same real 30-frame field-451
reference set, 8 worker processes: no crash, and a real 331.8s vs. the
~720s sequential baseline above — about **2.2x**, not full 8x
core-count scaling, since each worker still re-parses its own WCS header
per call and `pmap`'s own scheduling/serialization isn't free. Landed as
`build_reference`'s `workers` keyword (opt-in; the caller supplies
already-running worker processes, since spawning and managing a process
pool is an environment concern, not something a data-processing function
should own as a side effect).

The "tens of minutes" sequential cost is, as far as this investigation
could establish, a real, currently-irreducible property of reprojecting
one frame through `wcslib` at a time — but it is no longer irreducible
*in total*: spreading that same per-frame cost across independent
processes is a real, measured ~2x win, not something left unexplored for
lack of trying.
