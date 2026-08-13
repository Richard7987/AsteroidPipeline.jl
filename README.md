# AsteroidPipeline.jl

An open-source Julia pipeline for asteroid search campaigns, developed for
use with the [International Astronomical Search Collaboration
(IASC)](https://iasc.cosmosearch.org/).

## Purpose

The pipeline processes sequences of FITS frames from a survey field to:

1. **Detect** point sources in each frame (`detect_sources`).
2. **Link** detections across frames by consistent linear motion to form
   asteroid candidate tracklets (`link_candidates`).
3. **Cross-match** candidates against known-object catalogs — SkyBoT, VSX,
   SIMBAD — to separate previously cataloged objects from candidates that
   warrant human verification (`crossmatch_catalog`).

Planned extensions: variable-star and transient detection, and rotation
period recovery via Lomb-Scargle periodograms for confirmed discoveries.

## Status

Early development. The public API (`detect_sources`, `link_candidates`,
`crossmatch_catalog`) is defined but not yet implemented.

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
- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) — periodogram analysis

## License

MIT. See `LICENSE`.
