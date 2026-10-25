# Releasing a version

How a version goes out: registered in Julia's General registry first, then one
signed tag sets off everything else.

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
3. **Push the commit**, untagged, to Forgejo: `git push origin main`.
4. **Julia's General registry, first.** On GitHub, once the mirror has the
   commit, comment `@JuliaRegistrator register` on it. The JuliaRegistrator
   app opens a pull request in JuliaRegistries/General, whose automatic
   checks install and load the package on the oldest Julia it supports;
   once they pass, it is merged (after three days, the first time). If they
   fail, fix, push and comment again on the new commit: nothing has been
   tagged yet.
5. **Tag, after the registration is merged**, the registered commit, on
   Forgejo — tags must follow the registration, never precede it: a tag for
   a version the registry turned down names a release that doesn't exist.
   (TagBot, which tags automatically, can't be used: it pushes to GitHub,
   and the mirror would erase the tag.)

   ```sh
   git tag -s v1.0.0 <registered commit> -m "AsteroidPipeline 1.0.0 …"
   git push origin v1.0.0
   ```

6. **On GitHub, automatically**, once the mirror has the tag:

   | Workflow | Does |
   |:--|:--|
   | `Release.yml` | builds the app's wheel (with the library inside) and source distribution, checks the wheel installs, publishes both to PyPI, and attaches them with their checksums to a GitHub release |
   | `Documenter.yml` | rebuilds the documentation, every version (`docs/build_all_versions.sh`): `stable` becomes this release |

   `CI.yml` runs the tests on every push to `main`, on Julia 1.11 and the
   newest release, and uploads the newest's coverage to Codecov.

7. **Publish `stable`** with the next push to `main`. GitHub Pages names a
   deployment after its commit, and the tagged commit was already deployed
   from `main`, so the tag's own deployment reports success but leaves the
   site as it was: `stable` changes once a later commit is deployed.

## Set up once

- **PyPI**: a trusted publisher for the project `asteroidpipeline`
  (repository `Richard7987/AsteroidPipeline.jl`, workflow `Release.yml`,
  environment `pypi`) — no token is stored anywhere.
- **GitHub**: the environments `pypi`, and `github-pages` allowed to deploy
  from `main` and from `v*` tags; the JuliaRegistrator app installed on the
  repository.
- **Codecov**: the repository enabled at <https://app.codecov.io> (signed in
  with GitHub); uploads authenticate with GitHub's OIDC, no token stored.
