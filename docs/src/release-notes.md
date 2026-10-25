# Release notes

What changed in each version. Each release's documentation stays available
from the version picker at the top right.

## 1.1.0 — 2026-10-25

**Added**

- [`analyze_set`](@ref) takes `reference_stars`: a catalogue of one's own in
  place of the Gaia DR3 query, or work without the network.

**Changed**

- Tests: [`analyze_set`](@ref) run end to end on a synthetic set, and
  digest2 built in CI so [`digest2_score`](@ref) is tested there.
- Tests: the package is checked with Aqua.jl (method ambiguities, stale
  dependencies, compat bounds, type piracy).
- Docstrings and comments say what the code does and why, without the
  history of how each value was found; that lives in the validation pages.

**Fixed**

- SkyBoT moved to `ssp.imcce.fr`; its old address answers 403, so known
  objects were identified from the MPC's orbit file alone. The new service
  is used, with the observatory passed as `-observer` (it ignores `-loc` and
  computes geocentric positions) and an empty field (HTTP 204) read as no
  objects.
- [`julian_date_to_iso8601`](@ref) gave day 32 when rounding to the
  millisecond reached midnight on a month's last day; it now moves to the
  next month. The MPC date field gave `.000001` instead of `.000000` in the
  same case.
- The terminal UI kept the engine's log file open after starting it.

## 1.0.0 — 2026-10-23

The first stable release: from here on, a change that breaks the API raises
the major version. The images of any telescope are handled without
settings beyond the observatory code — checked on Pan-STARRS1's and ZTF's.

**Fixed**

- Images whose header astrometry is far off, or whose frames hold few
  matching Gaia stars, gave dozens of false candidates — every star seemed
  to move ~3″/h — and could report a known asteroid as new, 32″ from where it
  was. The frames' astrometry is now fitted together
  ([`refine_wcs_set`](@ref), used by [`run_pipeline`](@ref) and
  [`analyze_set`](@ref)): a wider search when a header is far off, stars
  matched again from each fit, and one solution per frame from the stars of
  all of them. On eight 2026-10-12/14 Pan-STARRS1 sets, 60 false candidates
  fell to none; on 20 sets, catalogued asteroids identified rose from 24 to
  28 ([Validation](validation/pan-starrs1.md)).
- Images from telescopes other than Pan-STARRS1, checked on ZTF's: the
  observation time is read from any of the usual header keywords
  (`DATE-OBS`, `MJD-OBS`, `OBSMJD`, `JD`…; [`frame_epoch`](@ref)) —
  [`analyze_set`](@ref) failed on ZTF images, which have no `DATE-OBS`; the
  astrometric searches are in arcseconds, for any pixel scale; and the MPC's
  orbit file is used only for images within a year of its orbits' epoch —
  propagated to 2019 images it put a known asteroid 155″ off, so older
  images use SkyBoT alone.
- Known asteroids could go unidentified, and pass for new ones: SkyBoT,
  asked for a 5″ cone, returned nothing for an asteroid 0.94″ from its
  centre in two of four frames. It is now always asked for at least 2′,
  and the distance is applied here.
- A SkyBoT outage stopped [`analyze_set`](@ref) with an error, instead of
  a warning: the failed request came wrapped in the error of the parallel
  task that made it. Known objects then come from the MPC's orbit file alone.
- A tracklet with one point too faint or too bright for the brightness test
  is kept without that point, if the rest passes every test: a faint known
  asteroid (2019 UA185, V ≈ 21) was being dropped for one frame at
  signal-to-noise ~5.

**Added**

- [`refine_wcs`](@ref) returns the matched stars (`matches`).
- Test coverage is measured in CI and published on Codecov.

**Changed**

- The package is in Julia's General registry (0.2.2, since 2026-10-22):
  `Pkg.add("AsteroidPipeline")`.
- Releases: a version is registered in Julia's General registry before its
  tag is created ([Releasing a version](development/releasing.md)).

## 0.2.2 — 2026-10-18

**Fixed**

- The package can be installed on Julia 1.11, the oldest it supports: the
  compat bound of `p7zip_jll` (which ships with Julia, at a fixed version per
  release) allowed only Julia 1.13's. The terminal app's analysis engine
  could not be set up on Julia 1.11 either.

**Changed**

- Tests run on Julia 1.11 and on the newest release.

## 0.2.1 — 2026-10-18

**Added**

- Your own temporary designations in reports — a team's numbering,
  `NHU0002` — with `designations` in [`mpc_report`](@ref),
  [`mpc80_report`](@ref) and [`ades_psv`](@ref).
- [`neocp_check`](@ref): candidates compared with the objects awaiting
  confirmation on the MPC's NEO Confirmation Page; [`analyze_set`](@ref) runs
  it for images up to 30 days old.
- The documentation is published per version: `stable`, each release, `dev`.
- Releases are published to PyPI (`pipx install asteroidpipeline`).

**Fixed**

- A bright asteroid could be dropped by the straight-line test when a
  field's astrometry is imprecise (few Gaia stars): the test now allows three
  times the frames' astrometric residual, at most 2 px. On ten real
  Pan-STARRS1 sets, two more catalogued asteroids are found
  ([Validation](validation/pan-starrs1.md)).

## 0.2.0 — 2026-10-18

The first release.

- [`analyze_set`](@ref): an image set — a `.zip`, a folder or FITS files —
  to an MPC report: astrometry refitted to Gaia DR3, damaged frames left
  out, the search, known objects identified against SkyBoT and the MPC's own
  orbit file ([`download_mpcorb`](@ref)), every known object in the field
  listed, and the report written.
- [`nearby_known_objects`](@ref): each candidate compared with every known
  object within 10′, for one moving the same way (a known object with an
  imprecise orbit).
- [`check_mpc_report`](@ref): a report checked before it is submitted.
- [`frame_site`](@ref), [`observatory_candidates`](@ref): the observatory's
  MPC code from the frames' FITS headers.
- Filters for hot pixels and bright-star light; [`set_quality`](@ref) for
  damaged image sets.
- A terminal app, `asteroidpipeline`, built on the library (see the
  repository's README).
- The documentation rewritten: tutorial, guides, explanation, validation
  reports.
