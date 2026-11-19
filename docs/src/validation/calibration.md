# Raw frames, reduced here

## Summary

[`plain_frame`](@ref) turns a camera's file into one plain image (a
tile-compressed `.fz` file, several extensions, amplifiers with their
overscan, binned readouts) and [`calibrate_frames`](@ref) subtracts bias and
dark and divides by a flat. Raw frames of Las Cumbres Observatory's public
archive were reduced here and compared with the same frames as the
archive's own pipeline (BANZAI) reduced them, each searched by
[`analyze_set`](@ref). On the 13 sets of the second and third runs, the
frames reduced here gave the same astrometric fits as BANZAI's pixels from
the same headers (within 0.01″ on every frame BANZAI solved itself), zero
points as uniform (scatter a median 0.000 mag from BANZAI's, at most 0.08
above), sky noise within 3% where measured, limiting magnitudes within 0.1, positions as
close to JPL's (median 0.070″ against 0.102″, and 0.210″ against 0.215″),
and the target's magnitudes within 0.01 mag. Of the 18 catalogued asteroids
found with BANZAI's frames, 14 were found here, and 1 that BANZAI's frames
missed. Of the 4 missing, 2 had the same signal to noise in both and were
at the detection threshold, 1 was lost with a frame whose astrometric fit
failed from its raw header (since fixed), and 1 was found and left
unidentified once, then identified when the same frames were searched
again. Two expectations failed: the solutions were not within 0.3″ of
BANZAI's (wrongly framed: BANZAI's own pixels were not either), and on two
of the four MuSCAT cameras the sky's signal to noise was 0.67–0.88 of
BANZAI's, from a different sky level, not more noise.

## 1. What is done

**`plain_frame`.** The first extension with a 2-D image, or a camera's
amplifiers assembled into one image by their `DATASEC` and `DETSEC`, each
with its overscan (`BIASSEC`) subtracted and its `GAIN` applied. `DETSEC`
counts unbinned detector pixels and a subarray's does not start at 1, so
`CCDSUM`'s binning and the subarray's origin are taken into account.

**`calibrate_frames`.** Master bias, dark and flat by median, the dark
scaled by exposure time, the flat normalised; then subtracted and divided.
Calibration frames read out in another mode than the frames (another size)
are refused with an error saying so.

## 2. Design

Declared before running, in `examples/calibration_test.jl`, which keeps
the full record.

**Data.** Sets of the archive whose target is a minor planet, with at
least 5 exposures of 15 s or more in one filter on one camera. Of each, 5
raw exposures, the same 5 as BANZAI reduced them, and raw calibration
frames of the same camera from the same night or the nearest of the 3
before it, read out in the exposures' mode: 5 bias frames, 3 darks and 5
sky flats in the filter. A set without flats in those nights is left out.

**Three reductions** of each set, each searched by `analyze_set` with its
defaults: "raw" (`plain_frame` only), "here" (`calibrate_frames`) and
"BANZAI".

**Measured.** Whether the set is reliable; the catalogued asteroids
identified; their positions against JPL Horizons' from the headers' site;
the target's magnitudes, here against BANZAI's; the limiting magnitude.
From the second run, also frame by frame, with the raw headers' WCS for all
three reductions: the astrometric fit, and its difference from BANZAI's
solution; the photometric zero point against Gaia and its scatter over the
frame's stars, which a flat field's errors raise; the sky's signal to noise.

## 3. First run: 14 sets, little tested

The nights of 2024-10-01, -04 and -09, the first 5 exposures of each set:
(19) Fortuna, (109) Felicitas and (187) Lamberta on 0.4 m telescopes
(QHY600 cameras), (203) Pompeja on 1 m telescopes (Sinistro, four
amplifiers per file). Three sets were left out (two without flats in the
filter, one by a refused Gaia query). All 11 were reliable in each
reduction, but the search tested little: on the 0.4 m telescopes the first
5 exposures span 1–2 minutes, so any main-belt asteroid is subtracted with
the stars, and Pompeja moves near the slowest rate 14 minutes can show.

The run found two faults. Calibration frames had been chosen by camera,
kind, filter and night only, and some were read out full-frame while the
exposures were a central subarray; they are now chosen in the exposures'
readout mode. And the astrometric fit failed, or only corrected an offset,
on every frame of 6 sets raw and 3 here, against 1 frame with BANZAI's: a
raw header holds the telescope's pointing, 9–106″ off and, where measured,
turned 0.37–0.67° and 0.15–0.6% out in scale, so one offset matches only the
stars near the centre. BANZAI's pixels failed the same way from the raw
headers, so the fault was the astrometry, not the reduction.
[`refine_wcs_set`](@ref) now also solves such a frame from its stars
([`solve_wcs`](@ref), near the header's centre and scale).

## 4. Second run: 4 sets

Other nights (2024-10-15 to 11-01), 5 exposures spread over the first hour
so that asteroids move. Of 6 sets, 2 were left out by the rules ((450649)
without flats in its filter, (3460) Shawna with fewer than 5 exposures in
an hour): (19) Fortuna, (201) Penelope, (617) Patroclus and 2024 TX13.

| Expected | Measured | |
|:--|:--|:--|
| frames refined wherever BANZAI's pixels are | 19 of 20 in each reduction, the same frame failing | met |
| solutions within 0.3″ of BANZAI's | a median 0.74″, at most 1.25″; BANZAI's own pixels 0.75″ and 1.24″ | not met, wrongly framed |
| zero-point scatter at most 0.01 mag above BANZAI's | 0.041 against 0.033 (median); the difference a median −0.000, at most +0.082 on one frame | met but one frame |
| raw frames' scatter above both | 0.059 | met |
| sky signal to noise at least 0.9 of BANZAI's | 0.98–0.99 | met |
| limiting magnitude within 0.1 of BANZAI's | +0.03, −0.02, −0.01 (Fortuna: no tracklets in any) | met |
| the same catalogued asteroids | 7 of BANZAI's 9 | not met |

BANZAI's solution is not the truth; against JPL the positions here were
the closer (35 points, median error 0.070″, errors over `pos_err` deviation
1.05; BANZAI's frames 45 points, 0.102″, 1.98). The target, Patroclus, was
identified in one set in each reduction, its magnitudes here a median
0.001 mag from BANZAI's.

Of the 2 asteroids missed here: **2010 KW118** has the same signal to noise
at its positions in both reductions (4–8 in a 4 px aperture); frame by frame
it is detected in 1 frame here and 2 in BANZAI's, and BANZAI's frames found
it only by stacking: it was at the threshold. **2001 DQ96** went with
Penelope's first frame, whose fit failed from the raw header in all three
reductions; that frame gave 0 detections here against 51 with BANZAI's,
which carry BANZAI's solution.

## 5. Changed after the second run

Both in `refine_wcs_set`, and checked on new data in the third run:

- A frame left with few stars starts again from the best frame's solution,
  carried over by the headers, but was searched only within 2″. Raw
  headers' pointing differs between exposures by more (Penelope's first
  frame 12″ off, the others 18–23″); an offset is now searched when that
  fails. Penelope's first frame is then refined in all three reductions,
  1.25″ from BANZAI's like the others.
- The frames' joint fit gave a frame without stars of its own the others'
  stars through the headers, unchecked: on a synthetic pair with headers
  15″ apart that frame came out "refined" and 15″ off. The joint fit now
  needs 3 stars of each frame's own.

On the 13 Pan-STARRS1 sets of the other validations, `refine_wcs_set` gives
the same result before and after every change made here (same frames
refined, same stars, the same solutions to the last digit), so those
validations stand.

## 6. Third run: 9 sets

Declared before it was run: the nights of 2024-11-05 to 11-19, a target a
minor planet if JPL's Small-Body Database calls it an asteroid, the
expectations as in the second run but with the astrometric one reframed (no
further from BANZAI's solution than BANZAI's pixels from the same headers,
plus 0.1″) and positions against JPL no worse than with BANZAI's frames. Of
13 sets, 4 were left out by the rules (three without flats in the filter,
one with fewer than 5 exposures in an hour). The 9: (156) Xanthippe on
0.4 m telescopes on four nights, (447755) on a 1 m, and the near-Earth
asteroid 2024 VC on the four cameras of MuSCAT (2 m, the four bands at
once).

| Expected | Measured | |
|:--|:--|:--|
| frames refined wherever BANZAI's pixels are | 45 of 45 in each reduction | met |
| solutions no further from BANZAI's than BANZAI's pixels' + 0.1″ | within ±0.01″ on the 24 frames BANZAI solved; see below for MuSCAT | met |
| positions against JPL no worse than with BANZAI's frames | median 0.210″ (32 points) against 0.215″ (38) | met |
| zero-point scatter at most 0.01 mag above BANZAI's | the difference a median −0.001, at most +0.011 | met |
| raw frames' scatter above both | 0.120 against 0.123 and 0.125 | not met |
| sky signal to noise at least 0.9 of BANZAI's | 0.91–1.12, but 0.88–0.90 and 0.67–0.77 on two MuSCAT cameras | not met on those |
| limiting magnitude within 0.1 of BANZAI's | −0.10 to +0.02 on the 6 sets with one in both | met |
| the same catalogued asteroids | 7 of BANZAI's 9, and 1 BANZAI's frames missed | not met |

**Astrometry on MuSCAT.** BANZAI marks its own solution failed
(`WCSERR = 4`) on all 20 MuSCAT frames and on one of Xanthippe's; its
solution there is the telescope's pointing, 10–132″ from the fits. On
MuSCAT the fits from the frames here and from BANZAI's pixels differed by
up to 14″, uniformly across a frame, both matching about 40 stars in all
four quadrants. The telescope tracked 2024 VC, which moved 1530″/h, 25″ in
each 60 s exposure: the stars are trails, which detection breaks into
pieces, and a fit can lock onto any point along them. `analyze_set` is not
made for frames tracked on a moving target.

**Zero points.** Through BANZAI's failed solutions no Gaia star could be
measured on MuSCAT's frames, so the zero points come from the other 5 sets.
There the raw frames' scatter was no higher than the calibrated ones', so
this measure could not show a flat field's work on these frames.

**Sky signal to noise on MuSCAT.** Split into its parts on the worst camera
(ep04), the noise here is 3% above BANZAI's, and the sky level is lower:
92–108 e⁻ against 134–136 e⁻. The level here is the raw frame's less the
night's bias; MuSCAT has no overscan, and BANZAI subtracts another level. The
measure declared mixes the two.

**Catalogued asteroids.** **2006 AE50** has the same signal to noise in both
reductions (3.8–6.1), is detected frame by frame in 1 frame in each, and
BANZAI's frames found it by stacking: at the threshold, as KW118 was.
**(10631) Koppes** was linked here and left unidentified, the one
unidentified tracklet of that set; searched again, the same frames give
Koppes identified, 4 points, as BANZAI's frames did. What differed is not
established: SkyBoT's answer is the only outside input, and no warning
showed it failing. **1999 VC5** was found here and raw, not with BANZAI's
frames. The target (447755) was identified in all three reductions, its
magnitudes here within 0.010 mag of BANZAI's (each within its error) and
its positions within 0.07″.

## 7. What it shows

Reduced here, raw frames of three kinds of camera (a single-image CMOS
camera, a four-amplifier CCD, a multi-band imager) give images that search
as BANZAI's do: the same astrometry from the same headers, zero points as
uniform, sky noise within a few percent, limiting magnitudes within 0.1,
and the same asteroids apart from those at the detection threshold, where
which reduction finds them is chance. Not tested: cameras with other
formats, darks of very long exposures, fringing in red bands, and frames
tracked on a moving object, on which the astrometry is ambiguous whatever
the reduction.
