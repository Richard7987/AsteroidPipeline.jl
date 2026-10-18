# AsteroidPipeline.jl

[![CI](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml)
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://richard7987.github.io/AsteroidPipeline.jl/stable/)
[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://richard7987.github.io/AsteroidPipeline.jl/dev/)

Find asteroids in astronomical images. Give it a few exposures of one field
taken minutes apart — a `.zip`, a folder or FITS files, from any telescope —
and it calibrates them against Gaia, finds what moves, tells known asteroids
from new candidates, and writes the report for the Minor Planet Center.

It comes two ways:

- **`asteroidpipeline`**, a terminal app: open an image set, watch the images
  animate with every object marked, confirm or reject each candidate, copy
  the report. Every decision is saved, with the candidate's pixels, as
  training data for a future classifier.
- **AsteroidPipeline.jl**, the Julia library that does all the analysis —
  for scripts, batches and anything the app doesn't do.

**Documentation: <https://richard7987.github.io/AsteroidPipeline.jl/stable/>** (for the latest release; [dev](https://richard7987.github.io/AsteroidPipeline.jl/dev/) for `main`)

## Install the app

On any Linux (or macOS) with Python 3.11 or newer:

```sh
curl -fsSL https://raw.githubusercontent.com/Richard7987/AsteroidPipeline.jl/main/install.sh | sh
asteroidpipeline
```

or with [pipx](https://pipx.pypa.io/), from a release's wheel:

```sh
pipx install https://github.com/Richard7987/AsteroidPipeline.jl/releases/download/v0.2.0/asteroidpipeline-0.2.0-py3-none-any.whl
asteroidpipeline --setup      # optional: install the analysis engine now
```

The first run installs the analysis engine: Julia (downloaded if you don't
have it) and its packages, several minutes once. Images are shown at full
resolution in terminals with a graphics protocol (Kitty, WezTerm, Ghostty,
foot…), as coloured blocks elsewhere. With Nix: `nix run github:Richard7987/AsteroidPipeline.jl`.

## Use the library

```julia
using Pkg; Pkg.add(url="https://github.com/Richard7987/AsteroidPipeline.jl")

using AsteroidPipeline
result = analyze_set("images.zip"; observatory="F51", workdir="work")
print(result.report)
```

Many sets at once: `julia --project=. examples/analyze_sets.jl <folder> <observatory code>`.
Start with [Getting started](https://richard7987.github.io/AsteroidPipeline.jl/dev/getting-started).

## Repository

| Path | What |
|:--|:--|
| `src/` | the Julia library |
| `test/` | its tests (`julia --project=. -e 'using Pkg; Pkg.test()'`) |
| `docs/` | the documentation site (`julia --project=docs docs/make.jl`) |
| `examples/` | scripts for real data, batches and completeness tests |
| `ui/` | the terminal app (Python, [Textual](https://textual.textualize.io/)); see [ui/README.md](ui/README.md) |
| `pyproject.toml`, `install.sh` | packaging of the app, with the library inside |
| `flake.nix` | Nix: `nix develop` for a development shell, `nix build` for the app |

## License

MIT. See `LICENSE`.
