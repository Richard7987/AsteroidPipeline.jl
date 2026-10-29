# Variable stars and light curves

The same detections that find moving objects also find sources that stay put
and change brightness.

## Find variables

```julia
r = search_field(paths; timestamp_key="OBSMJD", threshold=60.0)
r.movers        # what run_pipeline would return
r.variables     # stationary sources whose flux isn't constant
```

[`find_variable_sources`](@ref) matches sources by position across frames
and keeps those whose fluxes are inconsistent with a constant one,
by [`variability_chi2`](@ref) (reduced χ² against the mean, threshold 10).

Every flux measurement has a systematic error the photon noise doesn't
cover (flat-fielding, PSF changes); without allowing for it, bright stars
look variable. `systematic_error_fraction` (default 3%) adds it. The default
was chosen on three real fields, each with a confirmed variable star:

![False-positive rate and the confirmed variables' signal, against the systematic error allowed](../assets/systematic-error-floor-sweep.png)

At 3%, 0.6–2% of constant stars are flagged, and all three real variables
stay above the threshold; at 5%, none are flagged, but two of the three real
variables are lost too. Confirm a variable candidate (a catalogue match, a
period) before trusting it.

## Light curves and periods

For a source at a known position, over a dedicated series of exposures:

```julia
times, flux, flux_err = light_curve(paths, ra, dec)
p = recover_rotation_period(times, flux; minimum_period=0.02, maximum_period=1.0)
p.period, p.false_alarm_probability
```

[`recover_rotation_period`](@ref) runs a Lomb–Scargle periodogram
([LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl)). Bracket the
periods that make physical sense for the object, and cover more than one
full period: on real data, a light curve covering a third of a 0.32-day
period gave a highly significant — and wrong — 0.5-day period
([Validation: a variable star](../validation/variable-star.md)).

## A moving asteroid's light curve

[`moving_light_curve`](@ref) measures an asteroid where it is in each frame.
Give it the frames of a night and a few of the object's positions, such as
its tracklet from a set the pipeline analyzed:

```julia
r = analyze_set("set.zip"; workdir="work", observatory="F51")
tracklet = [c for c in r.candidates if c.id == 1]
lc = moving_light_curve(night_paths, tracklet)
good = .!lc.star_nearby .& isfinite.(lc.mag)
p = recover_rotation_period(lc.epoch[good], lc.mag[good]; minimum_period=0.05, maximum_period=1.0)
```

Its motion is fitted as a straight line in RA and Dec, which holds within a
night, and the aperture follows it into every frame. Where the object is
bright enough the aperture moves onto its centroid; where it isn't, the flux
is measured at the predicted position, so frames where the pipeline didn't
detect it still count. Each frame has its own zero point, so thin cloud
doesn't show up as a dip in the light curve, measured on the same comparison
stars in every frame: on unfiltered frames each star's colour shifts it by
its own amount, which a changing set of stars would turn into noise. On 40
unfiltered frames of 4775 Hansen from the 0.84 m telescope at San Pedro
Mártir, a fixed set of stars halved the scatter, from 0.016 to 0.008 mag.
`mag_err` includes the zero point's error.

Drop the frames marked `star_nearby`. An asteroid crossing a star looks
brighter: on ZTF frames, 2002 UY45 measured 16.8 on a star and 18.0–18.4
in the other frames.

An IASC set's four exposures in 45 minutes are too few for a period, which
needs a series covering one or more full rotations: hours, or several
nights.

## Several nights

From one night to the next the asteroid's distances from the Sun and from
Earth change, and so does its phase angle: its magnitude drifts by more than
its rotation changes it. [`reduced_magnitudes`](@ref) removes that drift,
with the geometry of each epoch from JPL Horizons
([`object_geometry`](@ref)), and moves each epoch back by the light time:

```julia
epochs, mags = lc.epoch[good], lc.mag[good]      # every night's, together
g = object_geometry(624, epochs; observatory="I41")
red = reduced_magnitudes(epochs, mags, g)
p = recover_rotation_period(red.epoch, red.mag; minimum_period=2/24, maximum_period=20/24)
```

The phase angle's effect follows the H-G model with `G=0.15` unless another
slope is given. Most asteroids' light curves have two maxima per rotation,
so the periodogram's highest peak is usually half the period: on 107 ZTF
frames of (624) Hektor over six months, it gave 3.460 h for a 6.9205 h
period.
