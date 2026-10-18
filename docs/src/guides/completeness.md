# Measure completeness

*Completeness* is the fraction of real moving objects a search finds. It
depends on the images (depth, seeing, crowding), so it is measured on them:
synthetic asteroids are added to the real frames, the pipeline is run, and
the recovered ones are counted.

```julia
injected, truth = inject_movers(paths, "injected/"; n=60, mag_range=(18.0, 23.0),
                                speed_range=(5.0, 100.0))        # "/h
found = run_pipeline(injected; pipeline_config(0.2563)...)
rec = injection_recovery(found, truth)
count(r -> r.recovered, rec) / length(rec)
```

[`inject_movers`](@ref) places each object through the frame's own refined
WCS, at a brightness set by the frame's own zero point, with the frame's
measured PSF, Poisson noise and the trail the object leaves during the
exposure. [`injection_recovery`](@ref) counts an object as recovered when one
tracklet follows it in at least three frames.

`examples/injection_test.jl` runs the full measurement and prints it by
magnitude and by rate. On a real Pan-STARRS1 set:

![Fraction of injected asteroids recovered, by magnitude and by rate](../assets/completeness.png)

About 80% are recovered down to V ≈ 21.5, falling to 16% at 22–22.5 and none
fainter. Fast objects are harder: above ~400″/h most cross the chip, or
trail beyond recognition, within one or two frames. Run it on your own
images: the numbers are theirs, not these.
