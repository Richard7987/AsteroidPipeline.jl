# Limitations

What the pipeline can't do, or doesn't do well — measured, with what would
fix it.

| Limitation | Measured | Why | What would help |
|:--|:--|:--|:--|
| Faint objects | ~80% found to V≈21.5; 16% at 22–22.5; none fainter (Pan-STARRS1, 45 s) | single-frame noise | a deep reference ([ZOGY](../guides/difference-imaging.md)), or stacking along trial motions |
| Slow objects | 50% found at 5–10″/h | they overlap their own positions in the set's median and partly subtract themselves | a reference from other nights |
| Fast objects | 44% at 100–400″/h, 7% above 700″/h | they leave the frame, or trail beyond recognition, within one or two frames | detection of trails |
| Objects on bright stars | lost in that frame (2 of 180 in a test) | masked as star light | — (3 frames out of 4 still give a tracklet) |
| Known objects with poor orbits | flagged `:possibly_known` only if displaced along the track | an orbit's error can be in any direction | — |
| Objects awaiting confirmation | NEOCP checked for images up to 30 days old, within 30′ | NEOCP predictions come from arcs of hours | — |
| Linking across nights | synthetic tests only | no real multi-night object in the test data | a real campaign's data |
| Variable stars | 0.6–2% of constant stars flagged | photometric systematics | confirm with a catalogue or a period |
| Network | needs Gaia, SkyBoT, JPL Horizons and the MPC | calibration and identification are against them | a service that fails is reported in `warnings`, and the rest goes on |

## Things that are not limitations

- **Any telescope** works whose images have, or can be given, a WCS: the
  settings scale with the pixel scale ([Use a different telescope](../guides/telescopes.md)).
- **Sets of any length** ≥ 3 frames.
- **Positions** agree with an independent reduction to 0.1–0.2″
  ([Validation](../validation/pan-starrs1.md)).
