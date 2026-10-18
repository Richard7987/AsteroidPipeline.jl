# A known variable star

## Summary

[`search_field`](@ref) finds a catalogued eclipsing binary in real ZTF images
as a variable source, 1.7″ from its catalogued position. The allowance for
systematic photometric errors, set on three fields with confirmed variables,
keeps all three detectable while flagging 0.6–2% of constant stars.

## 1. Introduction

Asteroid searches see variable stars too: a star that changes brightness
between exposures leaves a residual in a difference image. The pipeline's
variable-source search had only been checked against constant stars; this
tests it on a star known to vary.

## 2. Data

| | |
|:--|:--|
| Target | ASASSN-V J183620.31, contact eclipsing binary (EW), period 0.322 d ([VSX](https://www.aavso.org/vsx/)) |
| Images | ZTF field 487, CCD 12, quadrant 1, r band, 2019-06-10: 144 exposures over 2.45 h, 29 used |
| More fields | V1012 Mon and ASASSN-V J072906.85-090518.2, each with its own field |

The 2.45 hours cover about a third of one period, so recovering the period
from this night was not expected — only the variability.

## 3. Methods

`examples/variable_star_demo.jl`: [`search_field`](@ref), then cross-match of
the variable candidates with VSX; [`light_curve`](@ref) and
[`recover_rotation_period`](@ref) on the recovered star. The field is near the
Galactic plane (~12,900 sources per frame at 5σ), so detection used 60σ: the
target is far brighter than that.

## 4. Results

| | |
|:--|:--|
| detections | 638 |
| variable candidates | 22 |
| matched to VSX variables | 2, including the target, at 1.7″ |
| period from this night | 0.5 d, false-alarm probability ≈ 0 (catalogued: 0.322 d) |

![False-positive rate and the variables' significance, against the systematic error allowed](../assets/systematic-error-floor-sweep.png)

| Systematic error allowed | 0% | 1% | 2% | **3%** | 5% |
|:--|--:|--:|--:|--:|--:|
| constant stars flagged, field 451 | 13.2% | 3.3% | 2.0% | **1.3%** | 0.0% |
| constant stars flagged, V1012 Mon | 62.2% | 6.1% | 2.6% | **2.0%** | 1.0% |
| constant stars flagged, ASASSN-V J072906 | 49.7% | 2.8% | 0.6% | **0.6%** | 0.6% |

At 3%, all three variables stay above the detection threshold (reduced
χ² > 10); at 5%, two of them fall below it.

## 5. Discussion

The search finds a real variable among real constant stars, and the error
allowance controls its false positives. The 0.5-day period is what a light
curve covering a third of the true period gives — an alias, not a failure of
the search; recovering a period needs more than one full cycle. With 0.6–2%
of constant stars still flagged, a variable candidate is a lead to confirm
(a catalogue, a period), not a result.
