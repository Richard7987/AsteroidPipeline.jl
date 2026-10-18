# Write and check MPC reports

The Minor Planet Center (MPC) takes observations in two formats:

| Format | Function | When |
|:--|:--|:--|
| **ADES** (PSV) | [`ades_psv`](@ref) | the MPC's current format |
| **80-column** ("MPC1992") | [`mpc80_report`](@ref), [`mpc_report`](@ref) | still required by some programs and campaigns |

## A whole report

[`mpc_report`](@ref) writes the header, one line per observation, and the end
line; [`blank_mpc_report`](@ref), the report for "no moving objects":

```julia
header = (observers="J. Bulger, T. Lowe", measurers="A. Name, Team",
          telescope="1.8-m f/4.4 Ritchey-Chretien + CCD", net="Gaia DR3")
print(mpc_report(candidates, "F51"; trksub_prefix="P", header...))
```

```
COD F51
OBS J. Bulger, T. Lowe
MEA A. Name, Team
TEL 1.8-m f/4.4 Ritchey-Chretien + CCD
ACK MPCReport file updated 2026.10.18 15:03:46
NET Gaia DR3
     P2       C2025 09 16.41824600 12 25.816-10 43 50.18         21.4 G      F51
     P2       C2025 09 16.42926600 12 25.279-10 43 54.80         21.4 G      F51
     P2       C2025 09 16.44030300 12 24.750-10 43 59.11         21.7 G      F51
     P2       C2025 09 16.45135300 12 24.206-10 44 03.48         21.8 G      F51
----- end -----
```

| Header line | Holds | From |
|:--|:--|:--|
| `COD` | the observatory's MPC code — must come first | your observatory; see [Use a different telescope](telescopes.md) |
| `OBS` | who observed | the telescope's team |
| `MEA` | who measured | you |
| `TEL` | the telescope | its description |
| `NET` | the reference catalogue | Gaia DR3 (what positions are calibrated against) |

Each observation line carries a temporary designation (here `P2`, from the
tracklet id), the time of mid-exposure, the position and the magnitude.

### Your own designations

Observers usually number their finds with a prefix of their own. Give them
with `designations` (tracklet id ⇒ designation, up to 7 letters and digits):

```julia
mpc_report(candidates, "F51"; designations=Dict(2 => "NHU0002"), header...)
```

```
     NHU0002  C2025 09 16.41824600 12 25.816-10 43 50.18         21.4 G      F51
```

A candidate that is already on the MPC's NEO Confirmation Page goes under the
NEOCP's designation ([Identify known objects](known-objects.md)).
[`ades_psv`](@ref) takes the same keyword.
Known objects are reported by their packed designation:
[`pack_designation`](@ref)`("2016 SB115") == "K16SB5B"`.

## ADES

```julia
ades_psv(candidates, "F51")
```

gives a pipe-separated table, one row per observation, ready to write to a
`.psv` file; rows of one tracklet share its `trkSub`.

## Check a report before submitting it

[`check_mpc_report`](@ref) catches the mistakes that make a submission
bounce, or say something false — whatever program wrote the report:

```julia
c = check_mpc_report(read("MPCReport.txt", String); observatory="F51",
                     epochs=r.epochs, field=frame_footprint(r.wcs, r.image_size),
                     known=r.expected)
c.ok          # false if there is an error
c.issues      # (code, line, detail) for each problem
c.cleaned     # the corrected report
```

| Code | Problem | Error? |
|:--|:--|:--|
| `:missing_cod`, `:cod_not_first` | no `COD` line, or not first | yes |
| `:wrong_code`, `:station_mismatch` | another observatory's code | yes |
| `:other_set` | a line from another image set (other time, or other part of the sky) | yes |
| `:malformed_line`, `:empty` | not the 80-column format | yes |
| `:known_object` | a known asteroid measured by mistake | warning |
| `:single_observation` | an object observed once (the MPC needs two or more) | warning |
| `:missing_end` | no `----- end -----` | warning |

`epochs`, `field` and `known` are optional; each enables its checks. All
three come from an [`analyze_set`](@ref) result.
