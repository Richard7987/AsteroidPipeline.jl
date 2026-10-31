# Setting up for development

To change the code, run the tests or build this documentation. Any Linux
(or macOS) works; you need git and Julia 1.11 or newer —
[juliaup](https://julialang.org/install) installs it on any distribution:

```sh
curl -fsSL https://install.julialang.org | sh
```

## The library

```sh
git clone https://github.com/Richard7987/AsteroidPipeline.jl.git
cd AsteroidPipeline.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then, in the repository:

| To | Run |
|:--|:--|
| use your working copy | `julia --project=.`, then `using AsteroidPipeline` |
| run the tests | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| try it on real data | the scripts in `examples/` (`examples/fetch_data.sh` downloads public ZTF images) |

Some tests use the network (SkyBoT, JPL Horizons, the MPC); without it they
are skipped, not failed. Gaia, SIMBAD and VizieR, the slowest services,
answer from recorded answers in `test/fixtures`:

| To | Set |
|:--|:--|
| query them instead | `ASTEROIDPIPELINE_LIVE_TESTS=true` |
| record their answers again (after changing a query) | `ASTEROIDPIPELINE_RECORD=true` |

A query with no recorded answer fails with a message saying so. The tests
run on GitHub on Julia 1.11, with the recorded answers, and on the newest
release, live — keep the code working on both.

## The documentation

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The site lands in `docs/build/1`; `docs/build_all_versions.sh` builds every
version, as it is published. Pages are in `docs/src`; figures are drawn by
`docs/make_figures.jl` from the measured numbers, on every build.

## With Nix (optional)

The repository's `flake.nix` gives the same tools without installing them:
`nix develop` opens a shell with Julia and Python. Nothing above needs it.
