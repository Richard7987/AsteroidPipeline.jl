# Fainter objects by stacking

## Summary

[`stack_search`](@ref) adds each frame's light along trial motions, to find
moving objects too faint to be detected in any one frame. On four
Pan-STARRS1 fields with 80 injected asteroids each, it raised the recovered
ones from 82 to 101, all of the 19 added between V 21.4 and 23.0, with no
false candidate among them. On the 13 local sets without injections it added
two catalogued asteroids the frame-by-frame search had missed and four
candidates, of which two look like real faint objects and one is a diffuse
moving patch it can't tell from one. On a field whose noise is far from
Gaussian, it added nothing.

## 1. Method

The frames' differences from their median ([`stack_difference`](@ref)) are
filtered with a Gaussian of the measured FWHM, a matched filter for a point
source, and renormalized to unit noise. Peaks at 3σ in the first and last
frame are paired into motions, from twice the FWHM over the sequence (slower
objects subtract themselves in the median) to 100″/h, and every frame is
read where the motion puts the object. The score is the faintest frame's
value times √N, so every frame must show it.

Three things keep false candidates out, all measured on the set itself or
taken from the frame-by-frame search, none tuned to these images:

- **The threshold** is what noise reaches in the same frames: the search is
  repeated with the frames out of time order, where nothing real lines up,
  and the threshold is 1.1 times its best score, at least 7.
- **The artifacts** the frame-by-frame search drops: peaks on one detector
  pixel in several frames (a hot pixel, dithered) are masked, and a track
  near a masked pixel or over a bright star's light is left out.
- **One track per spot**: several tracks through one place are a patch of
  extended light, not a point source, and all are dropped.
- **Round light**: a track whose light within 4 FWHM is drawn out into a line
  (axis ratio above 1.7 in most frames) is a bright star's diffraction spike.
  Injected point sources measured 1.0–1.44, three spikes 1.85–2.16: the cut
  is between them, from those few cases.

[`analyze_set`](@ref) runs it after the frame-by-frame search and adds the
tracklets that search didn't find (`stacked` in its result; `stack=false`
turns it off).

## 2. Results

Four sets, 80 asteroids injected into each (V 21–24, 5–100″/h, seeds
101, 201, 202, 203), with [`inject_movers`](@ref):

| Field | Frame by frame | With stacking | Added (V) | False added | Positions, median / max |
|:--|:--|:--|:--|:--|:--|
| XY54_p10 | 20 | 26 | 22.1–23.0 | 0 | 0.13″ / 0.36″ |
| XY72_p01 | 22 | 27 | 21.6–22.3 | 0 | 0.31″ / 0.98″ |
| XY45_p01 | 14 | 22 | 21.4–22.7 | 0 | 0.14″ / 0.41″ |
| XY26_p01 | 26 | 26 | — | 0 | — |
| **Total** | **82** | **101** | | **0** | |

One more candidate on the XY45_p01 set is not injected: the one found on the
same field without injections, below.

The same fields without injections, and XY42_p11, the set with the most
detector defects:

| Set | Threshold | Tracklets added |
|:--|:--|:--|
| XY54_p10 | 7.4 | 0 |
| XY72_p01 | 8.8 | 0 |
| XY45_p01 | 7.1 | 1: 33″/h, faint in all four frames |
| XY26_p01 | 19.9 | 0 |
| XY42_p11 | 11.4 | 0 |

Every local set (IASC, Pan-STARRS1), before and after:

| Set | Known | Unknown | Added by stacking |
|:--|:--|:--|:--|
| XY15_p01 | 0 → 2 | 0 → 1 | the 2 known, and a diffuse patch at 5″/h, V 20.3 |
| XY25_p10 | 4 → 4 | 48 → 49 | 1: V 23.3, 41″/h, barely visible |
| XY45_p01 | 2 → 2 | 0 → 1 | 1: V 21.8, 33″/h, a faint point in each frame |
| XY45_p10 | 2 → 2 | 0 → 1 | 1: V 21.7, 34″/h, a faint point in each frame |
| the other 9 | unchanged | unchanged | 0 |

No catalogued asteroid was lost, and the damaged set is still flagged.
Before the round-light test, three diffraction spikes on XY72_p11 and
XY73_p01 were added too; it removed them and none of the injected asteroids.

A set's analysis took about 10–50 s longer, against runs on other days, so
with the catalogue services' own variation in it.

## 3. Discussion

**Where it gains.** On fields with Gaussian noise: the asteroids it adds are
at V 21.4–23.0, around and past the sets' own limits (V 22.4–23.0), where the
frame-by-frame search finds few.

**What it still lets through.** The diffuse patch on XY15_p01 moves with the
track, is round and leaves no light behind: it may be an optical ghost, which
moves as the telescope does, or a diffuse object. Nothing measured here
tells them apart; a person looking at it can.

**Where it doesn't.** XY26_p01's noise reaches a score of 18 with the frames
out of order, so the threshold rises to 19.9 and only bright objects pass,
which the frame-by-frame search already finds. That is the intended
behaviour, but it means the gain depends on the field.

**What it rests on.** The threshold is one number from the scrambled search,
its single best score: one slow artifact can raise it, and on a set where
the scrambled frames happen to miss an artifact, it can be too low. Before
the slowest motion searched was matched to the median's limit, XY26_p01's
threshold fell to 12.2 and seven false tracks passed.

**What wasn't tested.** Sequences of other lengths than four frames, and
telescopes other than Pan-STARRS1. With more frames, requiring every frame
loses objects that cross a gap; it doesn't add false ones.

## 4. Reproducing

`examples/stacking_test.jl` runs one field, with or without injected
asteroids, and prints what each search found:

```sh
julia --project=. examples/stacking_test.jl data/real/iasc/XY54_p10 101   # seed 101
julia --project=. examples/stacking_test.jl data/real/iasc/XY54_p10       # no injections
```
