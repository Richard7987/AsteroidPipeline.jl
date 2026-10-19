# AsteroidPipeline.jl

[![CI](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Richard7987/AsteroidPipeline.jl/actions/workflows/CI.yml)
[![PyPI](https://img.shields.io/pypi/v/asteroidpipeline.svg)](https://pypi.org/project/asteroidpipeline/)
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

**Documentation: <https://richard7987.github.io/AsteroidPipeline.jl/stable/>**
— for the latest release; the version picker at the top right has each
release's, and [dev](https://richard7987.github.io/AsteroidPipeline.jl/dev/)
for `main`. What changed: [release notes](https://richard7987.github.io/AsteroidPipeline.jl/stable/release-notes).

## What it does

- **Any image set**: a `.zip`, a folder or FITS files, from any telescope
  whose images carry a WCS; the observatory's MPC code is found from the
  FITS headers.
- **Calibrated**: every frame's astrometry refitted to Gaia DR3 (positions
  agree with an independent reduction to 0.13–0.19″), magnitudes against
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

## Install the app

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

With Nix: `nix run github:Richard7987/AsteroidPipeline.jl`.

The first run sets up the analysis engine — Julia (downloaded if you don't
have it) and its packages — which takes several minutes, once;
`asteroidpipeline --setup` does it ahead of time. Images are shown at full
resolution in terminals with a graphics protocol (Kitty, WezTerm, Ghostty,
foot…), as coloured blocks elsewhere. The app's own guide is
[ui/README.md](ui/README.md).

## Use the library

Julia 1.11 or newer:

```julia
using Pkg
Pkg.add(url="https://github.com/Richard7987/AsteroidPipeline.jl")
# Pkg.add("AsteroidPipeline") once its registration in Julia's General registry completes

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
| `flake.nix` | Nix: `nix develop` for a development shell, `nix build` for the app |
| `.github/workflows/` | tests, the documentation, and releases (PyPI, GitHub) |

Releasing a version — the tag, PyPI, the documentation, Julia's registry —
is described in [Releasing a version](https://richard7987.github.io/AsteroidPipeline.jl/dev/development/releasing).

## License

MIT. See `LICENSE`.
