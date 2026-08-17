# Validating `search_field` against a real, independently-confirmed variable star

Part of this project's [Investigation Log](investigation-log.md) —
split onto its own page since it's a self-contained validation story,
distinct from the ZTF field 451 bug-hunting narrative there.

Every real-data check so far had confirmed known *moving* objects (via
SkyBoT) but never a known *variable* — the one real dataset checked had
only one catalogued VSX variable in its footprint, too faint to serve as
a useful positive control. Fixed by choosing a real target with an actual
VSX catalog entry: ASASSN-V J183620.31 (type EW, a contact eclipsing
binary, period 0.322427 d), in ZTF field 487/CCD 12/quadrant 1/zr, night
2019-06-10 — a real high-cadence campaign with 144 exposures over 2.45 h,
thinned to 29 for `examples/variable_star_demo.jl`. Declared in advance:
2.45 h covers only ~32% of one period, so full period recovery from this
single night was not expected — the actual test was whether
`find_variable_sources` flags the star as variable at all.

First attempt did not run: `search_field`'s default `threshold=5.0` (used
at `8.0` in `real_data_demo.jl`) produced ~12,900 detections per frame
here, against field 451's ~130 — this field sits near the galactic plane.
`link_candidates`'s tracklet search is pairwise in detections/frame, and
did not finish in 35 minutes at that density; killed and confirmed via a
standalone check (`detect_sources` alone, one frame) that the counts were
real, not a hang. The target star is extremely bright at this field's
noise level (flux ≈ 87,000, ~1.7 px from its WCS-predicted position at
every threshold tested from 10 to 100), so raising the demo's threshold
to 60 — cutting detections/frame to ~1,260 — loses none of the signal
this run actually needs; a deliberate, field-specific tradeoff, not the
pipeline's own default. A second bug surfaced once linking finished:
`light_curve`'s default `timestamp_key="MJD-OBS"` doesn't exist in this
survey's headers (`search_field` was already correctly called with
`"OBSMJD"`) — a copy-paste omission in the demo script, not a pipeline
bug, fixed by passing the same key.

With both fixed, the run found 22 variable candidates from 638
detections; 2 matched a known VSX variable, including the target itself —
**ASASSN-V J183620.31, recovered at a 1.7" offset from its catalogued
position** — the positive-control criterion this whole exercise was
built to test, and it passed. A Lomb-Scargle fit to the recovered
candidate's own forced-photometry light curve (`light_curve` +
`recover_rotation_period`) found a real, highly significant periodic
signal (false-alarm probability ≈ 0) at 0.5 d, not the catalogued
0.322 d — the declared-in-advance outcome of fitting a period search to
a light curve covering less than a third of that period (almost
certainly an alias, not evidence against the true period), and not a
failure of `find_variable_sources` itself, which is what the crossmatch
recovery above actually validates.
