# Detection & Linking

Per-frame point-source detection (`detect_sources`, with
`fill_value_mask` for detector-gap fill values), difference imaging
against the median of a short sequence (`stack_difference` — the
numerical equivalent of a human's blink, for sets with no deep
reference), and cross-frame linear-motion matching into asteroid
candidate tracklets (`link_candidates`).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["detection.jl", "stack_difference.jl", "linking.jl"]
```
