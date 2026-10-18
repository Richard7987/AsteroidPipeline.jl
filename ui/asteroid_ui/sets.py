"""Image sets: .zip files waiting to be analyzed, and saved analyses."""

from __future__ import annotations

import json
import re
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from asteroid_ui.config import SETS_CACHE

FITS_SUFFIXES = (".fits", ".fit", ".fts")
_zip_cache: dict[tuple[str, float], bool] = {}


@dataclass
class SetEntry:
    name: str
    # what to analyze: a .zip, a folder, or FITS files (None: only the saved result is left)
    source: Path | list[Path] | None
    date: str
    result: dict[str, Any] | None

    @property
    def analyzed(self) -> bool:
        return self.result is not None

    @property
    def cache_dir(self) -> Path:
        return SETS_CACHE / self.name

    @property
    def source_label(self) -> str:
        if isinstance(self.source, list):
            return f"{len(self.source)} FITS files in {self.source[0].parent}"
        return str(self.source) if self.source else "(source no longer available)"


def is_fits(path: Path) -> bool:
    return path.suffix.lower() in FITS_SUFFIXES


def slug(text: str) -> str:
    """A safe folder name for the cache."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "set"


def zip_has_fits(path: Path) -> bool:
    key = (str(path), path.stat().st_mtime)
    if key not in _zip_cache:
        try:
            with zipfile.ZipFile(path) as z:
                _zip_cache[key] = sum(is_fits(Path(n)) for n in z.namelist()) >= 3
        except (OSError, zipfile.BadZipFile):
            _zip_cache[key] = False
    return _zip_cache[key]


def entry_for(source: Path | list[Path]) -> SetEntry:
    """The entry for a .zip, a folder of FITS images, or a list of FITS files."""
    if isinstance(source, list):
        first = source[0]
        name = slug(f"{first.parent.name}_{len(source)}files") if len(source) else "set"
        mtime = first.stat().st_mtime
    else:
        name = slug(source.stem if source.is_file() else source.name)
        mtime = source.stat().st_mtime
    date = datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")
    result = load_result(name)
    if result is not None:
        date = result_date(result) or date
    return SetEntry(name, source, date, result)


def load_result(name: str) -> dict[str, Any] | None:
    try:
        return json.loads((SETS_CACHE / name / "result.json").read_text())
    except (OSError, ValueError):
        return None


def result_date(result: dict[str, Any]) -> str:
    frames = result.get("frames") or []
    return frames[0]["time_utc"][:10] if frames else ""


def discover(folder: str, max_depth: int = 3) -> list[SetEntry]:
    """.zip files of FITS images under `folder`, and every saved analysis, newest first."""
    root = Path(folder).expanduser()
    found: dict[str, SetEntry] = {}
    zips: list[Path] = []
    if root.is_dir():
        for path in root.rglob("*.zip"):
            try:
                if len(path.relative_to(root).parts) <= max_depth and zip_has_fits(path):
                    zips.append(path)
            except OSError:
                continue
    for path in sorted(zips, key=lambda p: p.stat().st_mtime, reverse=True):
        entry = entry_for(path)
        found.setdefault(entry.name, entry)
    if SETS_CACHE.is_dir():
        for d in sorted(SETS_CACHE.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
            if d.name in found or (result := load_result(d.name)) is None:
                continue
            source = result.get("source")
            if isinstance(source, list):
                src: Path | list[Path] | None = [Path(p) for p in source] if all(Path(p).exists() for p in source) else None
            else:
                src = Path(source) if source and Path(source).exists() else None
            found[d.name] = SetEntry(d.name, src, result_date(result), result)
    return list(found.values())


def summary(result: dict[str, Any], decided: set[int] = frozenset()) -> str:
    """One line for lists: what the analysis found, and what is left to check."""
    if result["report_kind"] == "unreliable":
        return "damaged image set"
    new = [t["id"] for t in result["tracklets"] if not t.get("known")]
    n_known = sum(1 for k in result["known"] if k["detected"])
    if not new:
        return f"no new objects · {n_known} known"
    left = sum(1 for i in new if i not in decided)
    return (f"⚠ {left} of {len(new)} candidate(s) to check" if left else f"{len(new)} candidate(s) checked") + \
        f" · {n_known} known"
