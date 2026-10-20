# The pipeline, step by step

![The pipeline's stages](../assets/pipeline-overview.png)

An asteroid shows itself by moving: in a few images of the same field taken
minutes apart, stars stay put and the asteroid is somewhere else in each. The
pipeline turns that into a procedure — remove what stays put, find what is
left, keep what moves the way a real object moves, and drop what is already
known.

## 1. Clean and calibrate

| Problem in real images | What is done |
|:--|:--|
| invalid pixels flagged by a `BLANK` value, or gaps between detector cells filled with one constant | replaced, or masked ([`clean_blank_pixels`](@ref), [`fill_value_mask`](@ref)) |
| the header's astrometric solution is off by several arcseconds (6.8″ on one real Pan-STARRS1 set, 32″ on others) | the frames' WCS refitted to Gaia DR3 stars detected in them, all frames together ([`refine_wcs_set`](@ref)); residuals ~0.1″ |
| a frame full of detector defects | its 10σ sources counted per Gaia star; above 30, the frame is left out, and with fewer than 3 good frames the set is flagged ([`set_quality`](@ref)) |
| exposure times in different conventions | the time of mid-exposure, in UTC ([`frame_epoch`](@ref)) |

![Sources per reference star, on ten real image sets](../assets/set-quality.png)

## 2. Detect

Each frame is aligned to the others and the **median of the set** is
subtracted: stars and galaxies cancel, a moving object doesn't (it is in a
different place in each frame, so the median doesn't contain it). Sources are
detected at 6σ on the differences ([`detect_sources`](@ref),
[`stack_difference`](@ref)): with the stars gone, the threshold can be lower
than on a raw frame, where 8σ was needed not to drown in stars. Alternatively, a deep
reference built from other nights can be subtracted with ZOGY
([Deep references and ZOGY](../guides/difference-imaging.md)).

## 3. Filter

Differences leave things that aren't asteroids. Each kind is removed by a
test that was measured to separate it from real objects
([Telling real objects from artifacts](artifacts.md)):

| What | Removed by |
|:--|:--|
| single hot pixels, cosmic rays | `min_sharpness`: a real source is spread over the PSF |
| pieces of one trailed object | merged into one detection |
| detector defects (hot pixels, bad columns) | `detector_artifact_radius`: they stay put on the detector, while the sky moves with the telescope's dithers |
| halos, ghosts and diffraction spikes of bright stars | `max_background_excess`: they sit on extra light |

## 4. Link

Detections in different frames are joined into **tracklets** — the same
object seen several times — by [`link_candidates`](@ref): pairs of frames
give a trial motion; the other frames are searched where that motion
predicts, within `match_radius` (2″). A tracklet must then pass the tests
that real objects pass and coincidences don't:

| Test | Real objects | Why |
|:--|:--|:--|
| moves at 3–210″/h (`min_speed`, `max_speed`) | main belt ~20–40″/h, Trojans ~10–15″/h | slower is a star matched to itself |
| straight line, within 1 px RMS (`max_residual`) — or 3 times the astrometric residual, if larger (at most 2 px) | 0.1–0.3 px measured | chance alignments wander (median 2.6 px) |
| brightness constant within ×2.5 (`max_flux_ratio`) | | unrelated detections differ |
| seen in 3 or more frames (`min_frames`) | | a real object can fall in a gap in one frame |

## 5. Calibrate the tracklets

Pixel positions become right ascension and declination through each frame's
refitted WCS ([`astrometric_calibrate`](@ref)), and brightness becomes
calibrated magnitudes through each frame's own zero point against Gaia
([`candidate_magnitudes`](@ref), [`photometric_zeropoint`](@ref)).

## 6. Identify

Each tracklet is compared with the known asteroids predicted at its own
times ([`crossmatch_catalog`](@ref)): a match in three frames makes it known.
Every one left is then compared with all known objects within 10′, for one
moving exactly like it, and — for recent images — with the objects awaiting
confirmation on the MPC's NEO Confirmation Page
([Is a candidate new?](known-objects.md)). What remains are the
**candidates**: for a person to check, and then to report.

## 7. Report

[`mpc_report`](@ref) writes the candidates' observations in the MPC's
80-column format, or [`ades_psv`](@ref) in ADES; with no candidates,
[`blank_mpc_report`](@ref) writes the report for "no moving objects".
[Write and check MPC reports](../guides/reports.md).

## The settings

All of the above, as validated on real data, is
[`pipeline_config`](@ref)`(pixel_scale)`:

| Setting | Value | Meaning |
|:--|:--|:--|
| `threshold` | 6 | detection threshold, σ |
| `difference_stack`, `mask_fill_values`, `refine_astrometry` | on | steps 1–2 |
| `min_sharpness` | 0.3 | step 3 |
| `detector_artifact_radius` | 2 px | step 3 |
| `max_background_excess` | 1σ | step 3 |
| `match_radius` | 2″ | step 4 |
| `min_speed`, `max_speed` | 3″/h, ~210″/h | step 4 |
| `max_residual` | 1 px | step 4 |
| `max_flux_ratio` | 2.5 | step 4 |
| `min_frames` | 3 | step 4 |
