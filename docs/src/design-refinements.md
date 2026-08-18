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
