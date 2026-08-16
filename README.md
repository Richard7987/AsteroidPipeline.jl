# AsteroidPipeline.jl

[![Documenter](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/Documenter.yml)
[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://richard7987.github.io/AsteroidPipeline.jl/dev/)

An open-source Julia pipeline for asteroid search campaigns, developed for
use with the [International Astronomical Search Collaboration
(IASC)](https://iasc.cosmosearch.org/). Detects moving and variable
sources across a sequence of FITS frames, optionally via ZOGY difference
imaging, and calibrates and cross-matches the result against known-object
catalogs (SkyBoT, VSX, SIMBAD).

**Full documentation, including project status, real-data validation
results, known limitations, and the complete API reference, is at
<https://richard7987.github.io/AsteroidPipeline.jl/dev/>.**

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

## Contributing

`docs/make.jl` builds the documentation site
([Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) +
[DocumenterVitepress.jl](https://github.com/LuxDL/DocumenterVitepress.jl));
build it locally with:

```
julia --project=docs docs/make.jl
```

`.github/workflows/Documenter.yml` builds and deploys it to GitHub Pages
on every push to `main` and on tags.

## License

MIT. See `LICENSE`.
