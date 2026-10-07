# Astrometry & Plate-Solving

Pixel-to-sky calibration via WCS (`load_wcs`, `pix_to_sky`,
`astrometric_calibrate`), refinement of a header's WCS against Gaia DR3
stars detected in the frame itself (`gaia_reference_stars`, `refine_wcs`
— needed for real IASC/Pan-STARRS1 headers, off by ~7" on a real 2025
set), and `plate_solve` as a fallback for frames with no WCS in their
header at all (via the nova.astrometry.net API).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["astrometry.jl", "astrometric_refinement.jl", "platesolve.jl"]
```
