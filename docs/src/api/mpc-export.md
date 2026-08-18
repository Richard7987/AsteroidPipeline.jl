# MPC / ADES Export

Formats an [`astrometric_calibrate`](@ref) candidate table as an ADES
PSV observation table (`ades_psv`) — the format the Minor Planet Center
currently requires for astrometric submissions. Not the legacy 80-column
format; see `ades_psv`'s own docstring for why.

## Example

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

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["mpc_export.jl"]
```
