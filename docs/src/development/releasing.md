# Releasing a version

How a version goes out: one signed tag sets off everything else.

## The repositories

The project's own repository is a self-hosted Forgejo instance; GitHub
(`Richard7987/AsteroidPipeline.jl`) is a mirror of it, kept identical on
every push. Work is pushed to Forgejo only; GitHub runs the workflows.
Anything created only on GitHub's side as a git ref (a branch, a tag) is
erased by the next mirror sync — which is why nothing here pushes to GitHub
and why JuliaRegistries' TagBot is not used.

## Steps

1. **Version numbers.** The same version in `Project.toml` (`version`),
   `ui/asteroid_ui/__init__.py` (`__version__`), `install.sh` (`VERSION`)
   and the Status line of `docs/src/index.md`; an entry in
   `docs/src/release-notes.md`.
2. **Tests.** `julia --project=. -e 'using Pkg; Pkg.test()'`; the docs build
   (`julia --project=docs docs/make.jl`).
3. **Tag and push**, to Forgejo:

   ```sh
   git tag -s v0.2.2 -m "AsteroidPipeline 0.2.2 …"
   git push origin main v0.2.2
   ```

4. **On GitHub, automatically**, once the mirror has the tag:

   | Workflow | Does |
   |:--|:--|
   | `Release.yml` | builds the app's wheel (with the library inside) and source distribution, checks the wheel installs, publishes both to PyPI, and attaches them with their checksums to a GitHub release |
   | `Documenter.yml` | rebuilds the documentation, every version (`docs/build_all_versions.sh`): `stable` becomes this release |
   | `CI.yml` | the tests, on Julia 1.11 and the newest release |

5. **Julia's General registry.** On GitHub, comment
   `@JuliaRegistrator register` on the tagged commit. The JuliaRegistrator
   app opens a pull request in JuliaRegistries/General, whose automatic
   checks install and load the package on the oldest Julia it supports;
   once they pass, it is merged (after three days, the first time).

## Set up once

- **PyPI**: a trusted publisher for the project `asteroidpipeline`
  (repository `Richard7987/AsteroidPipeline.jl`, workflow `Release.yml`,
  environment `pypi`) — no token is stored anywhere.
- **GitHub**: the environments `pypi`, and `github-pages` allowed to deploy
  from `main` and from `v*` tags; the JuliaRegistrator app installed on the
  repository.
