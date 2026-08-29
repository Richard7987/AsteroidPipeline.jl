# MPC / ADES Export

Formats an [`astrometric_calibrate`](@ref) candidate table as either an
ADES PSV observation table (`ades_psv`) — the format the Minor Planet
Center currently requires for astrometric submissions — or the legacy
fixed-width 80-column format (`mpc80_report`), still required by some
programs (e.g. IASC, as of 2026) even though the MPC's own submissions
now prefer ADES.

## ADES example

```julia
ades_psv(candidates, "I41")
```

produces:

| trkSub | mode | stn | obsTime | ra | dec |
|:--|:--|:--|:--|--:|--:|
| 1 | CCD | I41 | 2000-01-01T12:00:00.000Z | 150.1234568 | 20.9876543 |
| 1 | CCD | I41 | 2000-01-01T12:00:00.864Z | 150.1235568 | 20.9877543 |
| 2 | CCD | I41 | 2000-01-01T12:00:00.000Z | 200.5000000 | -10.2500000 |

(shown as a table here; the real output is pipe-separated text, one line
per row, ready to write to a `.psv` file.) Both rows sharing `id=1` in
`candidates` share the same `trkSub`, which is how the Minor Planet
Center correlates them back into one tracklet.

## 80-column example

```julia
mpc80_report(candidates, "I41")
```

produces (one fixed-width 80-character line per row, `·` standing in
here for a literal space so the column boundaries stay visible):

```
·····1·······C2000·01·01.50000010·00·29.630+20·59·15.56·····························I41
·····1·······C2000·01·01.51000010·00·29.654+20·59·15.92·····························I41
·····2·······C2000·01·01.50000013·22·00.000-10·15·00.00·····························I41
```

Columns 1-5 (the permanent-number field) are always blank — this
pipeline's candidates are locally-numbered tracklets, never
MPC-numbered objects — and columns 6-12 carry the same `id`-derived
tracking label as ADES's `trkSub` (`1` for both rows above), just
capped one character shorter (7, not 8) and left-justified in a fixed
field instead of pipe-delimited. See `mpc80_report`'s own docstring for
the full column layout and the same real gaps `ades_psv` has (no
magnitude, no discovery asterisk).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["mpc_export.jl"]
```
