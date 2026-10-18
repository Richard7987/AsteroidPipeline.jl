"""People's decisions on candidates, kept as training data.

Every candidate someone confirms as real or rejects is written to the dataset
folder (config.DATASET_DIR) with its pixels: the tracklet's cut-out cube
(63 x 63 pixels x frames, FITS, from the analysis engine) and one row in
labels.csv. Known asteroids the program detected are added too, labelled
real by the catalogs — free positive examples. The goal is a classifier
trained on these that can one day stand in for the human check.

labels.csv columns:
  key           image set / tracklet id, unique
  set, tracklet
  label         real | bogus
  reason        why it was rejected (bogus): hot_pixel, bright_star, noise,
                other; for a real one the person recognised as a known
                object, known:<its name>
  source        human | catalog
  speed_arcsec_h, mag, n_points
  observatory, epoch_first
  stamp         path of the cut-out cube, relative to the dataset folder
  labelled_at   ISO time
"""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

from asteroid_ui.config import DATASET_DIR

REAL, BOGUS = "real", "bogus"
REASONS = {
    "hot_pixel": "Hot pixel or detector defect (stays put on the detector)",
    "bright_star": "Light of a bright star (halo, ghost, diffraction spike)",
    "noise": "Noise: nothing is really there",
    "other": "Something else",
}
_COLUMNS = ["key", "set", "tracklet", "label", "reason", "source", "speed_arcsec_h", "mag", "n_points",
            "observatory", "epoch_first", "stamp", "labelled_at"]
_INDEX = DATASET_DIR / "labels.json"


def _load_index() -> dict[str, dict[str, Any]]:
    try:
        return json.loads(_INDEX.read_text())
    except (OSError, ValueError):
        return {}


def _save_index(index: dict[str, dict[str, Any]]) -> None:
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    tmp = _INDEX.with_suffix(".tmp")
    tmp.write_text(json.dumps(index, indent=1))
    tmp.replace(_INDEX)
    with open(DATASET_DIR / "labels.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_COLUMNS)
        writer.writeheader()
        for row in index.values():
            writer.writerow({c: row.get(c, "") for c in _COLUMNS})


def labels_for(set_name: str) -> dict[int, dict[str, Any]]:
    """Tracklet id -> its label row, for one image set."""
    return {int(r["tracklet"]): r for r in _load_index().values() if r["set"] == set_name}


def _row(set_name: str, result: dict[str, Any], tracklet: dict[str, Any], cache_dir: Path,
         label: str, reason: str, source: str) -> dict[str, Any]:
    stamp_rel = f"stamps/{set_name}/tracklet_{tracklet['id']}.fits"
    stamp_src = cache_dir / tracklet.get("stamp", "")
    if stamp_src.is_file():
        (DATASET_DIR / stamp_rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(stamp_src, DATASET_DIR / stamp_rel)
    else:
        stamp_rel = ""
    points = tracklet["points"]
    return {
        "key": f"{set_name}/{tracklet['id']}", "set": set_name, "tracklet": tracklet["id"], "label": label,
        "reason": reason, "source": source,
        "speed_arcsec_h": round(tracklet["speed_arcsec_h"], 2),
        "mag": round(tracklet["mag"], 2) if tracklet.get("mag") is not None else "",
        "n_points": len(points), "observatory": result.get("observatory", ""),
        "epoch_first": points[0]["epoch"] if points else "", "stamp": stamp_rel,
        "labelled_at": datetime.now().isoformat(timespec="seconds"),
    }


def set_label(set_name: str, result: dict[str, Any], cache_dir: Path, tracklet_id: int,
              label: str | None, reason: str = "") -> None:
    """Record (or, with `label=None`, remove) a person's decision on a candidate."""
    index = _load_index()
    key = f"{set_name}/{tracklet_id}"
    if label is None:
        index.pop(key, None)
    else:
        tracklet = next(t for t in result["tracklets"] if t["id"] == tracklet_id)
        index[key] = _row(set_name, result, tracklet, cache_dir, label, reason, "human")
    _save_index(index)


def add_catalog_labels(set_name: str, result: dict[str, Any], cache_dir: Path) -> None:
    """After an analysis: drop the set's old labels (its tracklet ids changed) and
    add its known asteroids as real, confirmed by the catalogs."""
    index = {k: v for k, v in _load_index().items() if v["set"] != set_name}
    for t in result["tracklets"]:
        if t.get("known"):
            row = _row(set_name, result, t, cache_dir, REAL, "", "catalog")
            index[row["key"]] = row
    _save_index(index)


def dataset_counts() -> dict[str, int]:
    rows = _load_index().values()
    return {
        "real_human": sum(1 for r in rows if r["label"] == REAL and r["source"] == "human"),
        "real_catalog": sum(1 for r in rows if r["label"] == REAL and r["source"] == "catalog"),
        "bogus": sum(1 for r in rows if r["label"] == BOGUS),
    }
