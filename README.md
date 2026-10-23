# AsteroidPipeline.jl

[![CI](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://richard7987.github.io/AsteroidPipeline.jl/stable/)
[![Coverage](https://codecov.io/gh/Richard7987/AsteroidPipeline.jl/graph/badge.svg)](https://codecov.io/gh/Richard7987/AsteroidPipeline.jl)
[![PyPI](https://img.shields.io/pypi/v/asteroidpipeline.svg)](https://pypi.org/project/asteroidpipeline/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

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

**Documentation: <https://richard7987.github.io/AsteroidPipeline.jl/stable/>**
— for the latest release; the version picker at the top right has each
release's, and [dev](https://richard7987.github.io/AsteroidPipeline.jl/dev/)
for `main`. What changed: [release notes](https://richard7987.github.io/AsteroidPipeline.jl/stable/release-notes).

## What it does

- **Any image set**: a `.zip`, a folder or FITS files, from any telescope
  whose images carry a WCS; the observatory's MPC code is found from the
  FITS headers.
- **Calibrated**: every frame's astrometry refitted to Gaia DR3 (positions
  agree with an independent reduction to 0.12–0.17″), magnitudes against
  Gaia's.
- **The search**: each frame minus the set's median, artifacts removed (hot
  pixels, bright-star halos, ghosts and diffraction spikes), detections
  linked along straight lines at constant rates; damaged frames and sets
  recognised.
- **Is it new?** Known asteroids identified against SkyBoT and the MPC's own
  orbit file (downloaded daily); each candidate compared with every known
  object within 10′ and, for recent images, with the MPC's NEO Confirmation
  Page.
- **The report**: MPC 80-column or ADES, under your own designations
  (`NHU0002`) or the NEOCP's; any report checked before it is submitted.

How well it works, measured on real data, is in the
[validation reports](https://richard7987.github.io/AsteroidPipeline.jl/stable/validation/summary).

## Install

| | Command | Needs |
|:--|:--|:--|
| the app | `pipx install asteroidpipeline` | Python 3.11+ |
| the library | `julia -e 'using Pkg; Pkg.add("AsteroidPipeline")'` | Julia 1.11+ |

The app carries the library inside it and sets up Julia by itself; install
the library only to use it from your own code. Details:
[Installation](https://richard7987.github.io/AsteroidPipeline.jl/stable/installation).

### The app

On any Linux (or macOS) with Python 3.11 or newer, from PyPI with
[pipx](https://pipx.pypa.io/) (or pip):

```sh
pipx install asteroidpipeline
asteroidpipeline
```

or with the installer, which makes its own Python environment:

```sh
curl -fsSL https://raw.githubusercontent.com/Richard7987/AsteroidPipeline.jl/main/install.sh | sh
```

The first run sets up the analysis engine — Julia (downloaded if you don't
have it) and its packages — which takes several minutes, once;
`asteroidpipeline --setup` does it ahead of time. Images are shown at full
resolution in terminals with a graphics protocol (Kitty, WezTerm, Ghostty,
foot…), as coloured blocks elsewhere. The app's own guide is
[ui/README.md](ui/README.md).

### The library

Julia 1.11 or newer, from Julia's General registry:

```julia
using Pkg
Pkg.add("AsteroidPipeline")

using AsteroidPipeline
result = analyze_set("images.zip"; observatory="F51", workdir="work")
print(result.report)
```

Many sets at once: `julia --project=. examples/analyze_sets.jl <folder> <observatory code>`.
Start with [Getting started](https://richard7987.github.io/AsteroidPipeline.jl/stable/getting-started).

## Repository

| Path | What |
|:--|:--|
| `src/` | the Julia library |
| `test/` | its tests (`julia --project=. -e 'using Pkg; Pkg.test()'`), run on Julia 1.11 and the newest release |
| `docs/` | the documentation site (`julia --project=docs docs/make.jl`; every version: `docs/build_all_versions.sh`) |
| `examples/` | scripts for real data, batches and completeness tests |
| `ui/` | the terminal app (Python, [Textual](https://textual.textualize.io/)) |
| `pyproject.toml`, `install.sh` | packaging of the app, with the library inside |
| `flake.nix` | optional, for Nix users: a development shell and a package of the app |
| `.github/workflows/` | tests, the documentation, and releases (PyPI, GitHub) |

## Development

On any Linux (or macOS), with git, Python 3.11+ and Julia 1.11+ (the
easiest way to get Julia is [juliaup](https://julialang.org/install):
`curl -fsSL https://install.julialang.org | sh`):

```sh
git clone https://github.com/Richard7987/AsteroidPipeline.jl.git
cd AsteroidPipeline.jl

# the Julia library
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# the documentation
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl                  # site in docs/build/1

# the terminal app, editable: it uses this checkout's Julia code
python3 -m venv .venv && . .venv/bin/activate
pip install -e .
asteroidpipeline
```

With [Nix](https://nixos.org), `nix develop` gives a shell with Julia and
Python instead, and `nix run .` runs the app — optional, the steps above
need nothing of it.

Releasing a version — the tag, PyPI, the documentation, Julia's registry —
is described in [Releasing a version](https://richard7987.github.io/AsteroidPipeline.jl/stable/development/releasing).

## License

MIT. See `LICENSE`.
