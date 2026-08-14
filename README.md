# AsteroidPipeline.jl

An open-source Julia pipeline for asteroid search campaigns, developed for
use with the [International Astronomical Search Collaboration
(IASC)](https://iasc.cosmosearch.org/).

## Purpose

The pipeline processes sequences of FITS frames from a survey field to:

1. **Detect** point sources in each frame (`detect_sources`).
2. **Link** detections across frames by consistent linear motion to form
   asteroid candidate tracklets (`link_candidates`).
3. **Calibrate** tracklet pixel positions to sky coordinates using each
   frame's WCS astrometric solution (`load_wcs`, `pix_to_sky`,
   `astrometric_calibrate`).
4. **Cross-match** candidates against known-object catalogs — SkyBoT, VSX,
   SIMBAD — to separate previously cataloged objects from candidates that
   warrant human verification (`crossmatch_catalog`).

`run_pipeline` runs steps 1-3 end to end on a sequence of FITS file paths,
returning a candidate table ready for `crossmatch_catalog`.

Planned extensions: variable-star and transient detection, and rotation
period recovery via Lomb-Scargle periodograms for confirmed discoveries.

## Status

Early development. `detect_sources`, `link_candidates`, the WCS
calibration step, and `crossmatch_catalog` are implemented and wired
together end to end in `run_pipeline`, validated against synthetic FITS
frames with a known injected source track. Not yet run on a real IASC
dataset.

## Known limitations

- **No plate-solving fallback.** `load_wcs` only parses a WCS solution
  already present in the FITS header; it cannot derive one from an
  unsolved frame. IASC campaign frames are generally pre-solved, so this
  is not currently blocking, but it will need addressing before the
  pipeline can be used on frames from other sources (e.g. own
  blazar/exoplanet-timing imaging). See the `TODO` on `load_wcs` in
  `src/astrometry.jl` for candidate approaches.

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

A reproducible development environment is also provided via Nix
(`flake.nix`):

```
nix develop
```

## Dependencies

- [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) — FITS I/O
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) — source detection and photometry
- [WCS.jl](https://github.com/JuliaAstro/WCS.jl) — astrometric (pixel-to-sky) calibration
- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) — periodogram analysis
- [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl) — catalog cross-match queries

## License

MIT. See `LICENSE`.
