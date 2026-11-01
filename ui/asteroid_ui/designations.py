"""Temporary designations for the candidates you report.

The MPC's 80-column format has 7 characters for an observer's own
designation of an object it doesn't know yet. Teams number their finds with
a prefix of their own — "NHU0001", "NHU0002"… — set in Settings. Each
candidate keeps the designation it was first given, whenever its report is
written again, even if another analysis of the set numbers its tracklets
differently (candidates are told apart by their positions); numbers are
never reused.

A candidate on the MPC's NEO Confirmation Page is reported under the NEOCP's
designation instead: those observations help confirm someone's find.
Without a prefix, reports use the program's own: "P" and the tracklet's
number.
"""

from __future__ import annotations

import json
from typing import Any

from asteroid_ui.config import CONFIG_FILE, write_atomic
from asteroid_ui.sets import same_track, track_of

DESIGNATIONS_FILE = CONFIG_FILE.parent / "designations.json"
MAX_LENGTH = 7


def _load() -> dict[str, Any]:
    # only a missing file is empty: an unreadable one would be overwritten at
    # the next save, and reported candidates given new designations
    try:
        return json.loads(DESIGNATIONS_FILE.read_text())
    except FileNotFoundError:
        return {"assigned": {}}
    except ValueError as e:
        raise ValueError(f"{DESIGNATIONS_FILE} is damaged ({e}): repair it before reporting again") from e


def _save(data: dict[str, Any]) -> None:
    write_atomic(DESIGNATIONS_FILE, json.dumps(data, indent=1))


def valid_prefix(prefix: str) -> bool:
    return prefix == "" or (prefix.isalnum() and prefix.isascii() and len(prefix) <= MAX_LENGTH - 2)


def neocp_designation(tracklet: dict[str, Any]) -> str | None:
    check = tracklet.get("neocp") or {}
    if check.get("verdict") != "on_neocp":
        return None
    agreeing = [o for o in check.get("objects", []) if o.get("agreement") == "motion"]
    return agreeing[0]["designation"] if agreeing else None


def _given_before(set_name: str, tracklet: dict[str, Any], data: dict[str, Any], taken: set[str]) -> str | None:
    """The designation this object already has: by its positions, or by its id if given before they were kept."""
    track = track_of(tracklet)
    for designation, entry in data["tracks"].items():
        if entry["set"] == set_name and designation not in taken and same_track(entry["track"], track):
            return designation
    old = data["assigned"].get(f"{set_name}/{tracklet['id']}")
    return old if old is not None and old not in data["tracks"] and old not in taken else None


def assign(config, set_name: str, tracklets: list[dict[str, Any]]) -> dict[str, str]:
    """Designation for each of `tracklets` (tracklet id, as a string, => designation),
    giving new numbers to those that have none yet."""
    data = _load()
    data.setdefault("assigned", {})
    data.setdefault("tracks", {})
    out: dict[str, str] = {}
    prefix = config.designation_prefix.strip().upper()
    number = max(int(config.next_designation), 1)
    used = set(data["assigned"].values()) | set(data["tracks"])
    taken: set[str] = set()
    for t in tracklets:
        neocp = neocp_designation(t)
        if neocp:
            out[str(t["id"])] = neocp
            continue
        designation = _given_before(set_name, t, data, taken)
        if designation is None:
            if not prefix:
                continue
            digits = MAX_LENGTH - len(prefix)
            while True:
                if len(str(number)) > digits:
                    raise ValueError(f"{prefix} has room for {digits} digits; {number} doesn't fit")
                designation = f"{prefix}{number:0{digits}d}"
                number += 1
                if designation not in used:
                    break
            used.add(designation)
        taken.add(designation)
        data["assigned"][f"{set_name}/{t['id']}"] = designation
        data["tracks"][designation] = {"set": set_name, "track": track_of(t)}
        out[str(t["id"])] = designation
    _save(data)
    if number != max(int(config.next_designation), 1):
        config.next_designation = number
        config.save()
    return out
