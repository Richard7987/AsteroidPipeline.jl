# Follow-up tools

## Plate solving

An image with no WCS in its header is solved from where the telescope
pointed, which most headers record ([`frame_pointing`](@ref)): within 10′ of
it, [`solve_wcs`](@ref) finds the pixel scale (0.2–3″/pixel), the
orientation and a mirrored axis by matching triangles of stars against Gaia
DR3. [`run_pipeline`](@ref) and [`analyze_set`](@ref) do this by themselves
for every frame without a WCS:

```julia
img = permutedims(Float64.(read(FITS("frame.fits")[1])))
ra, dec = frame_pointing(read_header(FITS("frame.fits")[1]))
stars = gaia_reference_stars(ra, dec, 0.3; epoch_jd=2457618.8, mag_range=(6.0, 19.0))
r = solve_wcs(img, ra, dec, stars)       # r.wcs, r.scale, r.rms_arcsec
```

Frames from the 0.84 m telescope at San Pedro Mártir, whose headers hold
only IRAF's pixel coordinates, solved at 0.478″/pixel with 0.08–0.18″
residuals; ZTF frames with their WCS removed and the pointing 4′ off, at
1.014″/pixel, within 0.13″ of ZTF's own solution.

For a frame without its pointing in the header, [`plate_solve`](@ref) gets
a solution from [nova.astrometry.net](https://nova.astrometry.net/) (a free
account gives the API key):

```julia
plate_solve("frame.fits"; api_key=key)                 # one frame
run_pipeline(paths; plate_solve_api_key=key)           # where solve_wcs fails
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
