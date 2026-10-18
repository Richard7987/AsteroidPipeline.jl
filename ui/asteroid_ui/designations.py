"""Temporary designations for the candidates you report.

The MPC's 80-column format has 7 characters for an observer's own
designation of an object it doesn't know yet. Teams number their finds with
a prefix of their own — "NHU0001", "NHU0002"… — set in Settings. Each
candidate keeps the designation it was first given, whenever its report is
written again; numbers are never reused.

A candidate on the MPC's NEO Confirmation Page is reported under the NEOCP's
designation instead: those observations help confirm someone's find.
Without a prefix, reports use the program's own: "P" and the tracklet's
number.
"""

from __future__ import annotations

import json
from typing import Any

from asteroid_ui.config import CONFIG_FILE

DESIGNATIONS_FILE = CONFIG_FILE.parent / "designations.json"
MAX_LENGTH = 7


def _load() -> dict[str, Any]:
    try:
        return json.loads(DESIGNATIONS_FILE.read_text())
    except (OSError, ValueError):
        return {"assigned": {}}


def _save(data: dict[str, Any]) -> None:
    DESIGNATIONS_FILE.parent.mkdir(parents=True, exist_ok=True)
    DESIGNATIONS_FILE.write_text(json.dumps(data, indent=1))


def valid_prefix(prefix: str) -> bool:
    return prefix == "" or (prefix.isalnum() and prefix.isascii() and len(prefix) <= MAX_LENGTH - 2)


def neocp_designation(tracklet: dict[str, Any]) -> str | None:
    check = tracklet.get("neocp") or {}
    if check.get("verdict") != "on_neocp":
        return None
    agreeing = [o for o in check.get("objects", []) if o.get("agreement") == "motion"]
    return agreeing[0]["designation"] if agreeing else None


def assign(config, set_name: str, tracklets: list[dict[str, Any]]) -> dict[str, str]:
    """Designation for each of `tracklets` (tracklet id, as a string, => designation),
    giving new numbers to those that have none yet."""
    data = _load()
    assigned: dict[str, str] = data.setdefault("assigned", {})
    out: dict[str, str] = {}
    prefix = config.designation_prefix.strip().upper()
    for t in tracklets:
        neocp = neocp_designation(t)
        if neocp:
            out[str(t["id"])] = neocp
            continue
        if not prefix:
            continue
        key = f"{set_name}/{t['id']}"
        if key not in assigned:
            digits = MAX_LENGTH - len(prefix)
            number = max(int(config.next_designation), 1)
            if len(str(number)) > digits:
                raise ValueError(f"{prefix} has room for {digits} digits; {number} doesn't fit")
            assigned[key] = f"{prefix}{number:0{digits}d}"
            config.next_designation = number + 1
            config.save()
        out[str(t["id"])] = assigned[key]
    _save(data)
    return out
