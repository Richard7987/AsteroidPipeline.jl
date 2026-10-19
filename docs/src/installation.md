# Installation

AsteroidPipeline comes two ways. Install the one you need — or both.

| You want to… | Install | Needs |
|:--|:--|:--|
| open image sets, look at the candidates and copy the report, without writing code | the **app**, from PyPI (Python) | Python 3.11+ |
| analyze from scripts, process many sets, change the settings | the **library**, from Julia's registry | Julia 1.11+ |

The app carries the library inside it: installing the app is enough to use
the program, and it sets up Julia by itself.

## The app, with Python

```sh
pipx install asteroidpipeline
asteroidpipeline
```

[pipx](https://pipx.pypa.io/) installs it in an environment of its own, so it
doesn't mix with your other Python packages; `pip install asteroidpipeline`
works too. The first run sets up the analysis engine — Julia (downloaded if
you don't have it) and its packages — which takes several minutes, once;
`asteroidpipeline --setup` does it ahead of time.

Without pipx, the installer makes its own environment:

```sh
curl -fsSL https://raw.githubusercontent.com/Richard7987/AsteroidPipeline.jl/main/install.sh | sh
```

To upgrade: `pipx upgrade asteroidpipeline`.

## The library, with Julia

```julia
using Pkg
Pkg.add("AsteroidPipeline")
```

!!! note "Until 22 October 2026"
    The package's registration in Julia's General registry completes after
    its first-time waiting period, around 22 October 2026. Until then,
    install it from its repository:
    `Pkg.add(url="https://github.com/Richard7987/AsteroidPipeline.jl")`.

Then:

```julia
using AsteroidPipeline
result = analyze_set("images.zip"; observatory="F51", workdir="work")
```

[Getting started](getting-started.md) goes on from there. To upgrade:
`Pkg.update("AsteroidPipeline")`.

## Getting Python, pipx and Julia

Most Linux distributions have Python 3.11 or newer already. pipx and Julia:

| Distribution | pipx |
|:--|:--|
| Debian, Ubuntu, Mint | `sudo apt install pipx` |
| Fedora | `sudo dnf install pipx` |
| Arch, Manjaro | `sudo pacman -S python-pipx` |
| macOS (Homebrew) | `brew install pipx` |

Julia, on any of them, with [juliaup](https://julialang.org/install), its
official installer (which also keeps it up to date: `juliaup update`):

```sh
curl -fsSL https://install.julialang.org | sh
```

The app doesn't need it: it downloads Julia itself if there is none.

## Checking it works

| | |
|:--|:--|
| app | `asteroidpipeline --version` |
| library | `julia -e 'using AsteroidPipeline; println(pkgversion(AsteroidPipeline))'` |

Both need an internet connection when they run: positions are calibrated
against Gaia, and known asteroids come from the MPC and SkyBoT.

!!! note "Nix"
    With [Nix](https://nixos.org) the app is also
    `nix run github:Richard7987/AsteroidPipeline.jl` — an option, not a
    requirement: everything above works on any Linux.
