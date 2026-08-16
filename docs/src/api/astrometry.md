# Astrometry & Plate-Solving

Pixel-to-sky calibration via WCS (`load_wcs`, `pix_to_sky`,
`astrometric_calibrate`), and `plate_solve` as a fallback for frames with
no WCS in their header at all (via the nova.astrometry.net API).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["astrometry.jl", "platesolve.jl"]
```
