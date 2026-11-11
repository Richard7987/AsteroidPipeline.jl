# Limitations

What the pipeline can't do, or doesn't do well — measured, with what would
fix it.

| Limitation | Measured | Why | What would help |
|:--|:--|:--|:--|
| Faint objects | frame by frame, 82% found to V≈21.5; 21% at 22–22.5; none fainter (Pan-STARRS1, 45 s). [`stack_search`](@ref) adds objects at V 21.4–23.0 where the noise is Gaussian (82 → 101 of 320 injected), nothing where it isn't ([Fainter objects by stacking](../validation/stacking.md)) | single-frame noise | a deep reference ([ZOGY](../guides/difference-imaging.md)) |
| Slow objects | 50% found at 5–10″/h; in general, lost below about 4 × seeing ÷ the sequence's span ([`analyze_set`](@ref) warns when that passes 10″/h: a 5-minute sequence at 1.9″ seeing lost an asteroid at 83″/h) | they overlap their own positions in the set's median and partly subtract themselves | a longer sequence, or a reference from other nights |
| Fast objects | 57% at 100–400″/h, 20% at 400–1000″/h ([Pan-STARRS1](../validation/pan-starrs1.md)) | they leave the frame, or trail beyond recognition, within one or two frames | detection of trails |
| Objects on bright stars | lost in that frame (2 of 180 in a test) | masked as star light | — (3 frames out of 4 still give a tracklet) |
| Only 3 frames | on 6 Pan-STARRS1 sets with one frame removed: 100 unidentified candidates left, 52 of them on one set, after 3,940 of 4,053 are left out; the same 13 known objects and 92 injected ones found | no further frame confirms a tracklet, and their median removes less | more frames; [`analyze_set`](@ref) leaves out what is fainter than the set's limit, so a real object that faint is lost |
| Known objects with poor orbits | flagged `:possibly_known` only if displaced along the track | an orbit's error can be in any direction | — |
| Objects awaiting confirmation | NEOCP checked for images up to 30 days old, within 30′ | NEOCP predictions come from arcs of hours | — |
| Orbits | from 3 nights, a median 4.8% off in a; from 30 days, 0.03%; a sparse arc's error can exceed its linear 1σ many times ([Preliminary orbits](../validation/orbits.md)) | two-body, no planets; linear covariance | more nights; ranging or sampling of the uncertainty |
| Rotation periods from surveys | from one season of ZTF data, the period or its half within 1% for 13 of 29 asteroids chosen by a rule, against 3 of 29 with the magnitudes shuffled ([Rotation periods of thirty asteroids](../validation/period-sample.md)) | one visit a night: daily aliases, and a period can't be told from its half | several visits a night, or a dense light curve |
| Linking across nights | 3+ nights: all 47 real objects linked, none wrong, among 969 false tracklets; 2 nights: 110 wrong ([Linking across nights](../validation/linking.md)) | two tracklets often fit some orbit; a tracklet with one point beyond 3σ is left out | a third night; tested on catalogued objects only, not a campaign's candidates |
| Variable stars | 0.6–2% of constant stars flagged | photometric systematics | confirm with a catalogue or a period |
| Network | needs Gaia, SkyBoT, JPL Horizons and the MPC | calibration and identification are against them | a service that fails is reported in `warnings`, and the rest goes on |

## Things that are not limitations

- **Any telescope** works whose images have, or can be given, a WCS: the
  settings scale with the pixel scale ([Use a different telescope](../guides/telescopes.md)).
- **Sets of any length** ≥ 3 frames.
- **Positions** agree with an independent reduction to 0.1–0.2″
  ([Validation](../validation/pan-starrs1.md)).
