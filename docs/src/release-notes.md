# Release notes

What changed in each version. Each release's documentation stays available
from the version picker at the top right.

## Unreleased

**Fixed**

- The observer's positions from JPL Horizons were matched to the wrong epochs when
  asked out of time order: Horizons answers in time order. [`fit_orbit`](@ref) sorts
  its observations first and was not affected; [`orbit_ephemeris`](@ref) at unsorted
  epochs was.
- Documentation: the [limitations](concepts/limitations.md) and the
  [completeness guide](guides/completeness.md) quoted an earlier version's
  completeness (44% at 100–400″/h, 16% at V 22–22.5); the measured values are 57%
  and 21%. The limitations add rotation periods from surveys, and the Pan-STARRS1
  report says which astrometry its table of positions comes from.

## 1.5.1 — 2026-11-09

**Added**

- [Validation: rotation periods of thirty asteroids](validation/period-sample.md), chosen
  from the LCDB by a rule fixed before measuring (10 near-Earth, 10 Jupiter trojans,
  10 main belt): the period or its half within 1% for 13 of the 29 measured, against
  3 of 29 with each object's magnitudes shuffled.

**Fixed**

- Every request to JPL Horizons for a list of times failed on 2026-11-09, and with
  it [`analyze_set`](@ref) (the observer's position, to identify known asteroids),
  [`object_geometry`](@ref) and [`fit_orbit`](@ref). A list that gets no table back
  is now asked one time at a time, and the observer's position is read from
  Horizons' vector tables in text form, which kept their numbers when the CSV ones
  lost them.

## 1.5.0 — 2026-11-08

**Added**

- [`fit_orbit`](@ref) and [`orbit_ephemeris`](@ref): a preliminary two-body orbit
  from astrometry over several nights (Gauss's method and systematic ranging, then
  least squares), with its uncertainty, and where it puts the object later. On 40
  arcs of MPC astrometry of 15 known asteroids, a within a median 4.8% of JPL's from
  3 nights and 0.03% from 30 days ([Preliminary orbits](validation/orbits.md)).
- [`stack_search`](@ref) finds moving objects too faint for any one frame, by
  adding the frames along trial motions, with a threshold measured on the set's
  own noise; [`analyze_set`](@ref) adds what it finds (`stacked`). On four
  Pan-STARRS1 fields with injected asteroids: 82 recovered → 101, none false
  ([Fainter objects by stacking](validation/stacking.md)).
- [Validation: rotation periods of eight asteroids](validation/ztf-periods.md), chosen
  before looking at their ZTF data: two of the six measured within 1% of the
  published period, four on aliases of the nightly sampling or without a clear period.

**Changed**

- [`analyze_set`](@ref), when only 3 frames are searched, leaves out candidates
  fainter than the set's `limiting_mag` or with no magnitude, with a warning: nothing can confirm them,
  and on dense fields they ran to over a thousand, each asked of the catalogues.

**Fixed**

- [`analyze_set`](@ref) measured candidates with no magnitude on sets whose
  header is further off than one frame's refinement searches (`NaN`, and no
  magnitude in the report): it now uses the frames' joint fit.
- `limiting_mag` used a fixed 5-pixel aperture, too bright on telescopes with
  large pixels, where [`below_limit`](@ref) could reach real objects: it now
  uses an aperture of the frames' measured FWHM.

## 1.4.0 — 2026-11-01

**Added**

- [Validation: the light curve of an asteroid](validation/hektor.md): (624) Hektor
  on public ZTF images, 6.92035 ± 0.00008 h for a published 6.920509 h.
- [`analyze_set`](@ref) returns `limiting_mag`, from the frames' own noise and
  zero point, and [`below_limit`](@ref) the candidates fainter than it, most
  likely noise. The terminal app marks them and lists them last; none is discarded.
- A DOI, through Zenodo: [How to cite](citing.md) has it, and each release gets its own.
- [How to cite](citing.md) and `CITATION.cff`: the H-G model, the Moffat profile and VSX.

## 1.3.1 — 2026-11-01

**Added**

- `CITATION.cff` and [How to cite](citing.md): the reference for the software
  and for what it builds on.

## 1.3.0 — 2026-11-01

**Added**

- [`fourier_period`](@ref): an asteroid's rotation period from a Fourier
  series of its rotation phase, the way asteroid photometry measures it,
  with its error and amplitude. A period about as long as the light curve
  is warned of.
- [`gaia_reference_stars`](@ref) keeps every answer on disk
  ([`gaia_cache_dir`](@ref)), so a field analyzed or observed again needs
  no query, and falls back to a third copy of Gaia DR3, NOIRLab's Astro Data
  Lab in Tucson, when the ESA archive and VizieR (both in Europe) are down,
  as they were together for days in October 2026.
- [`moving_light_curve`](@ref): the light curve of a moving asteroid, measured
  where it is in each frame (its motion fitted from a few positions), with
  frames where it crosses a star marked, and the zero points from the same
  comparison stars in every frame. A frame with no astrometric solution is
  left out with a warning instead of stopping the curve. Frames are measured
  on every thread Julia has (`julia -t auto`); the terminal app starts its
  engine on every core.
- [`reduced_magnitudes`](@ref) and [`object_geometry`](@ref): magnitudes
  corrected for distance, phase angle and light time with the geometry
  from JPL Horizons, so several nights make one light curve.
- [`solve_wcs`](@ref): an astrometric solution for a frame whose header has
  none, from the telescope's pointing alone ([`frame_pointing`](@ref)) —
  pixel scale, orientation and a mirrored axis found against Gaia DR3,
  without a web service. [`run_pipeline`](@ref) and [`analyze_set`](@ref)
  use it for every frame without a WCS, before nova.astrometry.net; in a
  sparse field, a frame it can't solve starts from the solution of the frame
  nearest in time. Gaia is queried once per field: a header giving the
  pointing in two equinoxes, frame by frame, cost a query at every change
  (a 400-frame night took 10–15 min, now about 2). With `reference_stars`,
  the stars are taken from them, without the network.
- [`analyze_set`](@ref) warns when the frames span too little time for
  slow objects: below about 4 × seeing ÷ the span, an object is in the
  sequence's median and subtracts itself away.

**Changed**

- [`crossmatch_catalog`](@ref): while the CDS is down, `:simbad` asks
  SIMBAD's copy at Harvard and `:vsx` AAVSO's own VSX, with a warning.
  A service that failed is tried last for the next 15 minutes, as for Gaia.
- [`recover_rotation_period`](@ref): `period` is twice the periodogram's
  peak (now `peak_period`) where alternate cycles of the light curve
  differ, as an elongated asteroid's two maxima per rotation do: the peak
  is at half the period. Where the halves can't be told apart it is the
  peak, as before.
- Tests: Gaia, SIMBAD, VizieR, astrometry.net and the NEO Confirmation
  Page answer from recorded answers (`test/fixtures`), and the suite runs
  without waiting on them; [`plate_solve`](@ref) and [`neocp_check`](@ref)
  are tested end to end without an API key. The newest Julia in CI still
  queries them live ([Development setup](development/setup.md)).
- Shorter warnings.
- The Nix flake is no longer in the repository: the app and the library
  install the same way on every system ([Installation](installation.md)).
- [`pipeline_config`](@ref) takes `n_frames`: a tracklet needs half of
  them, and at least 3. [`analyze_set`](@ref) passes the set's; a 4-frame
  set is unchanged.

**Fixed**

- [`analyze_set`](@ref) warns when frames are left out as bad, and of the many
  false candidates to expect when only 3 remain.
- The terminal app, after a full review:
  - A candidate kept a designation by its tracklet's number, which another
    analysis of the set can change: two objects could share a designation.
    Designations and a person's decisions now follow the object's positions,
    numbers are never reused, and analyzing a set again asks first.
  - Images that don't name their telescope shared one profile, so another
    observatory's were reported under its code. They now ask every time.
  - A profile edited after the analysis reaches the report.
  - An analysis that failed on a set analyzed before showed the old results
    without a word.
  - A frame all of one value stopped [`analyze_set`](@ref): it is now left
    out as a bad frame. A set with a frame mostly of one value, or with an orbit without H in
    its field, failed at the very end of its analysis.
  - A .zip is no longer unpacked beside itself, which also failed in a
    read-only folder, and a set's copies of its frames are removed once analyzed.
  - The report check read UTF-8 as Latin-1, damaging accents in its corrected report.
  - Saving a report no longer replaces a different one beside it.
  - A damaged labels file, missing saved images, a wrong Julia path, an
    apostrophe in the install's path and brackets in a path no longer close
    the app or stop its set-up; `wl-copy` failing falls back on the next
    clipboard tool; `.ZIP` in capitals is found.
- MPC reports: accented letters in the header are written in TeX notation
  (`Olgu\\'{i}n`), as the MPC asks of names, and
  [`check_mpc_report`](@ref) warns of them and corrects them.
- [`extract_set`](@ref) and [`analyze_set`](@ref) took the hidden files macOS
  writes beside each file it zips (`._frame.fits`) for frames, and a set
  zipped on a Mac couldn't be analyzed.
- [`neocp_check`](@ref) skipped objects whose temporary designation has
  fewer than 7 characters, which the NEO Confirmation Page lists too: a
  candidate matching one was not recognised.
- While a copy of Gaia was down, every query waited out its timeouts
  before trying the next one: minutes for every frame of a sequence
  observed over many fields. A service that failed is tried last for the
  next 15 minutes.
- The terminal app treated a damaged labels or designations file as empty,
  and overwrote it at the next save: every label, or the designations
  already reported, lost. It now says the file is damaged and saves
  nothing; the designations file is written whole or not at all.
- The terminal app gave two sets of the same name (two `set.zip`, in different
  folders) one saved analysis, one set of labels and the same temporary
  designations: one of them was listed, and could be shown the other's
  candidates. Each now has its own.
- Known objects from the MPC's orbit file were lost on a set of more than
  about 50 frames: the observatory's positions were asked of JPL Horizons
  in one request, which it refuses (HTTP 502) past that, and the failure
  became a warning. They are asked 50 at a time.
- Two archives holding a folder of the same name (`images/`), unpacked in
  the same place, gave the second set the first one's frames:
  [`extract_set`](@ref) reused the folder. A folder holding other frames is
  left alone, and the set goes by the archive's name. Cleaned copies are
  kept by the set's frames as well as its name.
- One frame under cloud, in which no comparison star reached
  signal-to-noise 30, made [`moving_light_curve`](@ref) give up the
  comparison stars for the whole night, doubling its scatter; that frame
  alone now gets its own zero point.
- [`gaia_reference_stars`](@ref) gave a `NaN` position to a star whose
  proper motion a service writes as `NaN` instead of leaving it empty.
- Long sequences at short cadence, checked on 40 frames 30 s apart,
  which gave the asteroid in
  overlapping tracklets of at most 26 of its 40 frames, and false ones of
  3–4 frames: the hot-pixel filter took a slow mover, within 2 px of itself in the frames
  next to it, for a hot pixel (a hot pixel's recurrences now have to cover
  half the sequence's span), and a tracklet sharing two or more detections
  with a longer one is dropped as the same object. It is now one tracklet
  of 40 frames, and nothing else.
- [`load_wcs`](@ref) took the pixel coordinates IRAF writes into frames it
  processed (`WCSDIM`, `LTM1_1`…) for an astrometric solution, with no RA
  or Dec; it now requires a celestial one.
- [`refine_wcs`](@ref) no longer accepts a solution made of random matches
  on a damaged frame (smeared, or full of detector defects): thousands of
  spurious detections made a false offset win the vote, and the frame was
  reported refined with hundreds of matches. Such a frame now keeps its
  header WCS (`min_significance`). Found on 6 smeared ZTF frames of Hektor,
  where the magnitudes were wrong; the two damaged Pan-STARRS1 sets
  (`XY14_p10`, `XY73_p00`) were affected the same way. Every other frame
  checked, crowded Galactic-plane fields included, is refined as before.

## 1.2.0 — 2026-10-25

**Added**

- While SkyBoT is unavailable, known objects come from JPL's Small-Body
  Identification service: one query per field (a minute or two), its objects
  moved along their rates to each frame. Identification, the nearby-object
  check and the objects expected in the field keep working through a SkyBoT
  outage, with one short warning instead of one per search.

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
