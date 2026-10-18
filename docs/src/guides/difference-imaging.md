# Deep references and ZOGY

Two ways to take the static sky out of the images before searching them:

| | Median of the set itself (default) | Deep reference + ZOGY |
|:--|:--|:--|
| Needs | only the set | ~30 earlier images of the same field |
| How | each frame minus the median of all frames (`difference_stack=true`) | each frame minus a deep stack, PSF-matched ([Zackay, Ofek & Gal-Yam 2016](https://doi.org/10.3847/0004-637X/830/1/27)) |
| Good for | short sets (3–5 frames), any telescope | going fainter than one frame's noise; long series |
| Weakness | an object slower than ~10″/h partly subtracts itself | building the reference takes tens of minutes |

## Median differencing

On by default in [`pipeline_config`](@ref). Each frame is aligned to the
others through its WCS (whole pixels; `subpixel_alignment=true` interpolates,
at a measured cost in completeness), the median is subtracted, and sources
are detected on what is left.

## ZOGY with a deep reference

```julia
ref = build_reference(reference_paths)              # 30 frames: tens of minutes
candidates = run_pipeline(science_paths; reference=ref, threshold=5.0)
```

[`build_reference`](@ref) reprojects and stacks the reference frames;
[`estimate_psf`](@ref) measures the PSF from bright isolated stars (falling
back to a [`fit_moffat_psf`](@ref) model if there are too few);
[`zogy_subtract`](@ref) returns the difference and its significance map,
`S_corr`, which is what sources are detected on.

`quality_max_std` (default 1.5) leaves out a frame whose `S_corr` noise is
far from 1 — on real ZTF data, the one frame taken through thin cloud.
[Validation: ZTF](../validation/ztf.md) has the numbers.
