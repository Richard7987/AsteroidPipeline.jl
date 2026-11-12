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

Frames whose headers hold only IRAF's pixel coordinates have been solved
this way; ZTF frames with their WCS removed and the pointing 4′ off, at
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

[`link_nights`](@ref) groups tracklets from different nights that one orbit
fits: the candidates of a field's sets from several nights, each set's
table as one night.

```julia
sets = [analyze_set(dir; observatory="F51", workdir=joinpath("work", dir)) for dir in ("night1", "night2", "night3")]
linkages = link_nights([s.candidates for s in sets], "F51")
for l in linkages
    println(l.tracklets, "  a = ", l.orbit.elements.a, " AU, rms ", l.orbit.rms, "″")
end
```

Each linkage lists its tracklets as `(night, id)` and carries the
[`fit_orbit`](@ref) result for their observations together. Pass the
astrometry's uncertainty as `sigma` (0.5″ by default): every observation must
lie within 3σ of the orbit, except one point of a tracklet of 3 or more, which
is then left out of the orbit and listed in `dropped`. Tested on real astrometry of 73 asteroids over
12 nights, with and without false tracklets: [Linking across nights](../validation/linking.md).

Linkages need 3 nights by default (`min_nights`): two tracklets moving like
asteroids often fit some orbit, so a 2-night linkage is a guess to
check, not an identification. An object seen on one night only can't be
linked at all: its field has to be observed on two more nights, each within
`max_days` (15) of another.

[`link_across_nights`](@ref), the earlier linker, extrapolates each tracklet
in a straight line instead; it needs no orbit and no network, but its
tolerance (5″ by default) has not been tested on real data.
