# Other telescopes

## Summary

[`analyze_set`](@ref) was run with its defaults, nothing set for the
instrument, on 15 public image sets of three telescopes other than
Pan-STARRS1. The three differ in aperture, pixel scale and filter: ZTF
(1.2 m Schmidt, 1.01″/px), the Las Cumbres Observatory's 0.4 m and 1 m
telescopes (0.74 and 0.39″/px), and DECam on the Blanco 4 m (0.263″/px). On
each, 100 synthetic asteroids were added per set, and the catalogued
asteroids in the frames were checked against JPL.

The results:

- **Depth.** On all three, half of the synthetic asteroids are found frame
  by frame 0.8–1.0 mag brighter than the frames' 5σ limit, as on
  Pan-STARRS1 (0.89).
- **Catalogued asteroids** are found as often as synthetic ones of the same
  brightness.
- **Positions and their uncertainties** hold on synthetic asteroids, on
  every instrument.
- **Against JPL**, the catalogued asteroids' positions hold on DECam and
  Las Cumbres but not on ZTF: a deviation of 1.67 on 23 points.
- **False candidates.** Expectation 3, at most 2 unidentified tracklets per
  set, fails on every instrument. The failure is concentrated in crowded
  fields. One DECam set in z near the galactic plane left 2,998 tracklets,
  took 8 hours, and recovered none of its 100 synthetic asteroids.

The test led to six fixes in the package. The first run's magnitude
uncertainties were 3–5 times too small. The candidate photometry was
replaced, validated on Pan-STARRS1, and every set was run again. Both runs
are reported.

## 1. Data

The sets were chosen and downloaded by `examples/telescope_data.jl`, by
rules fixed before any of their images was looked at:

- **ZTF**: October 2023, 3 sequences at galactic latitude |b| ≥ 20° and 3 at
  |b| < 10°. Each is one quadrant of 4–5 exposures of 30 s. IRSA's survey
  limit `MAGLIM` is given for comparison.
- **Las Cumbres Observatory**: frames reduced by the archive's own pipeline
  (BANZAI), 1–14 December 2024, 3 sets from 0.4 m and 3 from 1 m telescopes.
  Each set is one target, camera, filter and night, 5 exposures spanning
  30 min or more, with the telescope tracking the stars.
- **DECam**: frames calibrated by the community pipeline (InstCal) in March
  2023. Each set is 5 exposures and one of DECam's 59 CCDs, both drawn at
  random.

All 15 sets are within 15° of the ecliptic.

## 2. Design

Declared before running, in `examples/telescope_test.jl`.

**Search.** [`analyze_set`](@ref) with its defaults (stacking on), the
instrument's MPC code (for Las Cumbres, the one nearest each set's site),
and the MPC's orbits.

**Measured on the frames as they are:**

- the 5σ limiting magnitude L;
- the catalogued asteroids predicted on the frames, found or not;
- the unidentified tracklets;
- the catalogued asteroids' points against JPL Horizons.

**Measured with 100 synthetic asteroids added per set**
([`inject_movers`](@ref), V from L − 2 to L + 0.5, 5–100″/h):

- the fraction found frame by frame and with stacking, by V − L;
- V50, the V − L at which half are found, with sets resampled for its
  interval;
- their errors over the uncertainties given.

**Expected,** if the search works as on Pan-STARRS1:

1. half found frame by frame between L − 1.2 and L − 0.6 (the median over
   the sets);
2. catalogued asteroids found as often as synthetic ones of the same
   brightness, with stacking (their 95% intervals overlap);
3. at most 2 unidentified tracklets per set, on average;
4. synthetic asteroids' position errors over `pos_err` with a deviation of
   at most 1.5, and 90% or more within 2;
5. catalogued asteroids' errors against JPL with a deviation of at most 1.5.

## 3. Changes made during the test

Each change was declared in the script's header when made, with what had
been seen. Five of them are fixes in the package:

- **The injector read the wrong extension.** [`inject_movers`](@ref) read
  only a file's first extension. DECam's and Las Cumbres' files hold the
  image in the second, tile-compressed. It now reads such files as
  [`analyze_set`](@ref) does. This was found on the first DECam set, before
  any result.
- **The recovery count matched by frame number.**
  [`injection_recovery`](@ref) matched points by frame number. The injector
  numbers frames in the order given, the search in time order. DECam's files
  are named by checksum and so do not sort by time, so 0 of 100 were counted
  recovered although the search's tracklets lay 0.07–0.46″ from them. It now
  matches by epoch.
- **A request waited forever.** During DECam's second set a SkyBoT request
  waited 3 h for an answer that never came: HTTP.jl's idle timeout did not
  end the TLS read. SkyBoT and JPL requests now end at a deadline of their
  own and are retried.
- **A failing Gaia copy was not passed over.** One night all three copies
  of Gaia DR3 failed. ESA's archive answered 401, then its maintenance page
  behind status 200, which was cached as a catalogue; VizieR timed out and
  NOIRLab answered 502. A copy that refuses, or answers with a page, is now
  passed over for the next. The 5 sets that failed that night were run
  again.
- **No limit without candidates.** The limiting magnitude was measured
  only on frames with candidates, so a set with none got no limit, and its
  injection was skipped (ZTF f0382_c15). It is now measured on every frame.

One change is in the test itself:

- **Points without an uncertainty.** Points whose aperture measured no flux
  had no uncertainty and made expectations 4 and 5 NaN. They are now left
  out of both, and counted.

After the first run, the candidate photometry was replaced (section 5) and
every set was run again with the same seeds. The expectations were not
changed.

## 4. Results

### By set, second run

V50 is relative to L and computed frame by frame; it is shown only where
the sets are counted.

| Set | Frames | L | ZTF's limit | Unidentified | V50 frame by frame | V50 with stacking | Time |
|:--|--:|--:|--:|--:|--:|--:|--:|
| **ZTF** f0382_c06 (\|b\| < 10°) | 5 | 19.00 | 19.14 | 34 | −0.71 | −0.71 | 590 s |
| f0382_c13 (\|b\| < 10°) | 4 | 18.85 | 18.79 | 17 | −1.30 | −1.27 | 932 s |
| f0382_c15 (\|b\| < 10°) | 4 | 18.91 | 19.20 | 0 | −1.54 | −1.54 | 288 s |
| f0468_c05 | 5 | 19.64 | 20.25 | 0 | −0.51 | −0.51 | 138 s |
| f0574_c06 | 5 | 19.48 | 19.54 | 0 | −0.56 | −0.56 | 138 s |
| f0621_c10 | 4 | 20.38 | 20.84 | 0 | −1.42 | −1.42 | 117 s |
| **Las Cumbres** AD Leonis, H-alpha | 5 | 16.73 | | 5 | — | — | 230 s |
| ASAS J064437, V | 5 | 18.98 | | 2 | −0.56 | −0.56 | 367 s |
| K2-25, r′ | 5 | 20.11 | | 0 | −0.76 | −0.27 | 84 s |
| TIC 118555069, i′ | 5 | 18.58 | | 1 | −0.79 | −0.62 | 100 s |
| TOI-3884, 12 Dec, r′ | 5 | 22.08 | | 6 | — | −1.25 | 381 s |
| TOI-3884, 14 Dec, r′ | 5 | 20.70 | | 13 | −1.26 | −1.26 | 391 s |
| **DECam** 15 Mar, r, CCD 37 | 5 | 23.77 | | 0 | −0.85 | −0.23 | 184 s |
| 18 Mar, z, CCD 9 | 5 | 20.03 | | 2,998 | — | — | 28,489 s |
| 21 Mar, r, CCD 32 | 5 | 23.77 | | 1 | −0.95 | −0.24 | 141 s |

A dash is a V50 that the fractions found do not define: half is never
reached, or is reached more than once.

### Expectations, both runs

The first run's value is given first, the second's after the arrow.

| | ZTF | Las Cumbres | DECam |
|:--|:--|:--|:--|
| 1. V50 frame by frame, median | −0.71 (5 sets) → −1.00 (6 sets): met | −0.78 → −0.78: met | −0.90 → −0.90: met |
| 2. catalogued against synthetic | met → met | met → met | met → met |
| 3. unidentified per set | 10.2 → 8.5: **not met** | 4.5 → 4.5: **not met** | 1,000 → 1,000: **not met** |
| 4. synthetic, positions: deviation, within 2 | 1.38, 91% → 1.36, 92%: met | 1.49, 90% → 1.47, 90%: **not met** | 1.16, 92% → 1.19, 92%: met |
| 5. catalogued against JPL: deviation (points) | 1.67 (23) → 1.67 (23): **not met** | 1.44 (27) → 1.20 (27): met | 1.40 (8) → 1.19 (12): met |
| synthetic, magnitudes: deviation (robust) | 4.04 (1.11) → 3.06 (0.73) | 5.51 (1.17) → 1.18 (0.85) | 3.05 (1.17) → 0.88 (0.84) |

On Las Cumbres, expectation 4's deviation is within 1.5 but the fraction
within 2 falls just under 90% (it rounds to 90%).

In the first run, ZTF's V50 is the median over 5 sets because f0382_c15 had
no limit (section 3). DECam's catalogued points rose from 8 to 12 between
the runs: the asteroids were the same, and the cause was not traced.

### Fraction found by V − L, second run, all sets together

| V − L | ZTF: injected | frame by frame | with stacking | Las Cumbres: injected | frame by frame | with stacking | DECam: injected | frame by frame | with stacking |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| −2.00 to −1.75 | 68 | 54% | 54% | 68 | 59% | 79% | 30 | 60% | 63% |
| −1.75 to −1.50 | 52 | 67% | 67% | 57 | 67% | 79% | 30 | 73% | 73% |
| −1.50 to −1.25 | 63 | 62% | 62% | 53 | 60% | 74% | 34 | 74% | 74% |
| −1.25 to −1.00 | 54 | 67% | 70% | 60 | 55% | 60% | 26 | 38% | 42% |
| −1.00 to −0.75 | 63 | 59% | 59% | 70 | 37% | 43% | 29 | 31% | 48% |
| −0.75 to −0.50 | 60 | 42% | 42% | 62 | 24% | 42% | 35 | 11% | 46% |
| −0.50 to −0.25 | 66 | 17% | 17% | 59 | 7% | 25% | 20 | 5% | 45% |
| −0.25 to 0.00 | 56 | 4% | 4% | 50 | 0% | 12% | 31 | 0% | 23% |
| 0.00 to +0.25 | 51 | 0% | 0% | 63 | 0% | 5% | 28 | 0% | 11% |
| +0.25 to +0.50 | 67 | 0% | 0% | 58 | 0% | 0% | 37 | 0% | 8% |

DECam's columns include the z set, where none was found.

## 5. The magnitudes, and the photometry replaced

In the first run, the synthetic asteroids' magnitude errors over their
uncertainties had a deviation of 3–5.5, while the robust deviation was
1.1–1.2. A few points were far off. On DECam these were points up to 1.3 mag
too bright: a star in the aperture.

The candidate photometry was then changed in two ways:

- **What lies still under a candidate is taken out.** The same aperture is
  measured where the candidate was in the frames where its tracklet has
  moved clear of it (by twice the aperture plus the FWHM). The median of
  those, scaled by the frames' zero points, is subtracted, with its error.
- **The noise is measured, not computed.** The aperture's noise is measured
  in 400 apertures placed at random in empty places of the frame. The
  frame's pixel scatter missed the sky's correlations.

The change was declared, then checked on the 13 Pan-STARRS1 sets of
[Uncertainties](uncertainties.md) with the same synthetic asteroids
(section 2 there, fourth run), before this test was run again. There it met
all four expectations:

- magnitude deviation from 2.35 to 1.06, robust from 1.24 to 0.97;
- points beyond 5σ from 0.9% to 0.1%;
- positions 0.95 → 0.91.

On the other telescopes, the second run brought DECam to 0.88 (robust 0.84)
and Las Cumbres to 1.18 (0.85). ZTF fell only to 3.06, with a robust 0.73:
its bulk is right, but a few points remain far off. The cause has not been
found.

## 6. What it means

**What holds.** On instruments of 1.01 to 0.263″/px and apertures of 0.4 to
4 m, with nothing set for them, the search reaches the same depth relative
to the frames' limit as on Pan-STARRS1. It finds catalogued asteroids as
often as synthetic ones. Its position uncertainties hold on synthetic
asteroids. Against JPL they hold on DECam and Las Cumbres, on small samples
(12 and 27 points).

**False candidates in crowded fields.** The search leaves too many there:

- ZTF's f0382_c06 and f0382_c13 left 34 and 17;
- Las Cumbres' TOI-3884 field on 14 December left 13;
- DECam's z set left 2,998 and found none of its synthetic asteroids.

Elsewhere the counts are 0–6 per set. Expectation 3 fails on every
instrument, and this is the next thing to measure and fix.

**ZTF.** Even bright synthetic asteroids are found only 54–70% of the time.
Its sets span 13–22 min, so objects slower than 38–69″/h move less than four
seeings across the frames and are subtracted with the stars: `analyze_set`
warns of this on every ZTF set. Its catalogued positions are off by more
than their uncertainties (robust deviation 1.75) on 23 points of 3 sets. Its
magnitudes have tails the new photometry did not remove.

**Not tested here:** fields far from the ecliptic, and instruments other
than these three.
