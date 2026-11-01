"""Paths and preferences.

Everything the user can change lives in ~/.config/asteroid-ui/config.json and
in the Settings screen; the rest follows from where this package sits in the
repository.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

_XDG_CONFIG = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
_XDG_CACHE = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
_XDG_DATA = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
CONFIG_FILE = _XDG_CONFIG / "asteroid-ui" / "config.json"
CACHE_DIR = _XDG_CACHE / "asteroid-ui"
SETS_CACHE = CACHE_DIR / "sets"
BACKEND_LOG = CACHE_DIR / "engine.log"
# Every candidate a person has labelled, with its pixels: training data for a
# future classifier. Kept apart from the cache, which can be thrown away.
DATA_DIR = _XDG_DATA / "asteroid-ui"
DATASET_DIR = DATA_DIR / "dataset"


@dataclass
class Config:
    # Shown big on the start screen: your team's name, or the program's.
    title: str = "NEO Hunters"
    # Folder searched (with its subfolders) for image-set .zip files.
    downloads: str = str(Path.home() / "Downloads")
    # The MPC report's MEA and CON lines: who measures (you, your team) and
    # how to reach them. The telescope's own lines (COD, OBS, TEL) come with
    # its profile (profiles.py).
    measurers: str = ""
    contact: str = ""
    # Your own temporary designations for the candidates you report: a prefix
    # ("NHU") and the next number (NHU0002). Empty: the program's own (P2).
    designation_prefix: str = ""
    next_designation: int = 1
    # A report file to check before submitting it (from this program or any other).
    report_to_check: str = ""
    # Animation speed, images per second (at most 3: no fast flashing).
    fps: float = 1.5

    def save(self) -> None:
        write_atomic(CONFIG_FILE, json.dumps(asdict(self), indent=2, ensure_ascii=False))


def write_atomic(path: Path, text: str) -> None:
    """Write the whole file or none of it: a program stopped halfway leaves the old one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def load_config() -> Config:
    try:
        data = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return Config()
    if not isinstance(data, dict):
        return Config()
    defaults = Config()
    # a value of the wrong type, edited by hand, is left at its default
    return Config(**{k: v for k, v in data.items() if hasattr(defaults, k)
                     and isinstance(v, (int, float) if isinstance(getattr(defaults, k), float) else type(getattr(defaults, k)))})
