"""The analysis engine: AsteroidPipeline.jl, run by Julia.

The Julia package travels inside this Python package (julia/AsteroidPipeline,
put there when the wheel is built; in a checkout of the repository, the
repository itself is used). The first time the app runs — and again after an
upgrade — it gets an environment ready for it:

1. find Julia: $ASTEROIDPIPELINE_JULIA, `julia` on the PATH, or, if neither,
   a private copy that `juliapkg` downloads;
2. make a Julia environment in the user's data folder with AsteroidPipeline
   and the few packages server.jl needs, and precompile it (several minutes,
   once).

After that, server.jl starts in seconds.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from collections.abc import Callable
from pathlib import Path

from asteroid_ui import __version__
from asteroid_ui.config import BACKEND_LOG, DATA_DIR

ENGINE_DIR = Path(__file__).resolve().parent
SERVER = ENGINE_DIR / "server.jl"
_BUNDLED = ENGINE_DIR.parent / "julia" / "AsteroidPipeline"
_CHECKOUT = ENGINE_DIR.parents[2]          # ui/asteroid_ui/engine -> repository root
_SERVER_PACKAGES = ["FITSIO", "JSON", "WCS", "TypedTables", "Dates", "Statistics"]
MIN_JULIA = (1, 11)


class EngineError(Exception):
    """The engine can't be found or set up; the message says what to do."""


def julia_package() -> Path:
    """Where AsteroidPipeline.jl's source is: bundled in the wheel, or this checkout."""
    for candidate in (_BUNDLED, _CHECKOUT):
        if (candidate / "Project.toml").is_file() and (candidate / "src" / "AsteroidPipeline.jl").is_file():
            return candidate
    raise EngineError("AsteroidPipeline.jl's source is missing from this installation; reinstall it.")


def env_dir() -> Path:
    return DATA_DIR / "engine"


def julia_executable() -> str:
    explicit = os.environ.get("ASTEROIDPIPELINE_JULIA")
    if explicit:
        return explicit
    found = shutil.which("julia")
    if found:
        return found
    try:
        import juliapkg  # downloads a private Julia if there is none
    except ImportError:
        raise EngineError("Julia isn't installed. Install it from https://julialang.org/install "
                          "(or `pip install juliapkg` to have it downloaded) and start the app again.") from None
    return juliapkg.executable()


def julia_version(julia: str) -> tuple[int, ...]:
    out = subprocess.run([julia, "--startup-file=no", "-e", "print(VERSION.major, ' ', VERSION.minor)"],
                         capture_output=True, text=True, timeout=120)
    try:
        return tuple(int(x) for x in out.stdout.split())
    except ValueError:
        raise EngineError(f"Could not run Julia ({julia}): {out.stderr.strip()[:300]}") from None


def _fingerprint(package: Path) -> str:
    """Changes whenever the engine must be set up again: app version, package sources' manifest, Julia."""
    h = hashlib.sha256(__version__.encode())
    h.update((package / "Project.toml").read_bytes())
    h.update(str(package).encode())
    for f in sorted((package / "src").glob("*.jl")):
        h.update(f.name.encode())
        h.update(str(f.stat().st_size).encode())
    return h.hexdigest()


def needs_setup() -> bool:
    try:
        stamp = (env_dir() / "ready").read_text().strip()
    except OSError:
        return True
    return stamp != _fingerprint(julia_package())


def setup(julia: str, on_line: Callable[[str], None] = lambda line: None) -> None:
    """Make the engine's Julia environment; lines of Julia's output go to `on_line` and the log."""
    version = julia_version(julia)
    if version < MIN_JULIA:
        raise EngineError(f"Julia {'.'.join(map(str, version))} is too old; AsteroidPipeline needs "
                          f"{'.'.join(map(str, MIN_JULIA))} or newer.")
    package = julia_package()
    env = env_dir()
    if env.exists():
        shutil.rmtree(env)
    env.mkdir(parents=True)
    script = (
        "using Pkg; "
        f"Pkg.develop(path={str(package)!r}); "
        f"Pkg.add({_SERVER_PACKAGES!r}); "
        "Pkg.instantiate(); Pkg.precompile()"
    ).replace("'", '"')
    BACKEND_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(BACKEND_LOG, "a", buffering=1) as log:
        log.write(f"--- setting up the engine in {env} with {julia}\n")
        proc = subprocess.Popen([julia, f"--project={env}", "--startup-file=no", "-e", script],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            log.write(line)
            if line.strip():
                on_line(line.rstrip())
        if proc.wait() != 0:
            raise EngineError(f"Setting up the analysis engine failed; details in {BACKEND_LOG}")
    (env / "ready").write_text(_fingerprint(package))


def server_command(julia: str) -> list[str]:
    return [julia, f"--project={env_dir()}", "--startup-file=no", "--threads=auto", str(SERVER)]
