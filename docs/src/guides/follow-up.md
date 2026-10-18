# Follow-up tools

## Plate solving

For an image with no WCS in its header, [`plate_solve`](@ref) gets one from
[nova.astrometry.net](https://nova.astrometry.net/) (a free account gives
the API key):

```julia
plate_solve("frame.fits"; api_key=key)                 # one frame
run_pipeline(paths; plate_solve_api_key=key)           # every frame without a WCS
```

It uploads the image and waits for the solution: allow a minute or more per
frame.

## Orbit class: digest2

[`digest2_score`](@ref) runs the MPC's [digest2](https://github.com/Smithsonian/digest2)
on each tracklet: the likelihood that an object with this motion is a
near-Earth object, a main-belt asteroid, and so on. It is a separate C
program:

```sh
git clone https://github.com/Smithsonian/digest2.git && cd digest2/digest2 && make
```

```julia
digest2_score(candidates, "I41"; digest2_path="digest2/digest2/digest2")
```

| id | rms | int_score | neo_score | n22_score | n18_score |
|--:|--:|--:|--:|--:|--:|
| 1 | 0.02 | 100.0 | 100.0 | 23.0 | 0.0 |

!!! warning "It scores motion, not reality"
    digest2 says what kind of object a tracklet would be *if it is real*.
    On real ZTF data, star "tracklets" — a star's centroid jittering by less
    than 1″ — scored 100 as NEOs, while the two real main-belt asteroids
    scored 5 and 7. Score tracklets only after the false ones are gone.

## Linking across nights

[`link_across_nights`](@ref) groups tracklets from different nights that fit
one object moving in a straight line, allowing for curvature with
`max_accel_arcsec_per_day2`:

```julia
nights = [search_field(p; timestamp_key="OBSMJD").movers for p in (night1, night2, night3)]
groups = link_across_nights(nights)            # each group: (night, id) pairs
```

!!! note "Not yet tested on real data"
    It recovers objects exactly in synthetic tests, but no real object seen
    on several nights in the same field has been available to test it on.
    The defaults (5″, no curvature allowance) are starting points.
