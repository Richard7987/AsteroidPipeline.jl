"""Telescope profiles: the report header's lines that come with the images.

Every image set's frames say which telescope took them (TELESCOP, INSTRUME...
in their FITS headers) but not its MPC code, its observers or its full
description. A profile holds those for one telescope: the first time a
telescope's images arrive, the engine suggests the code (from the site in the
headers and the MPC's list of observatories) and the person confirms it; from
then on every set from that telescope gets its profile on its own. The
measurers (MEA) are whoever does the analysis, so they are in Settings, not
here.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any

from asteroid_ui.config import CONFIG_FILE

PROFILES_FILE = CONFIG_FILE.parent / "profiles.json"


@dataclass
class Profile:
    key: str             # which telescope, as its images name it
    observatory: str     # MPC code (COD)
    observers: str = ""  # OBS
    telescope: str = ""  # TEL
    builtin: bool = False

    def describe(self) -> str:
        return f"{self.observatory} — {self.key}" + (" (built in)" if self.builtin else "")


# Pan-STARRS 1, as IASC's image sets come: the header IASC's own Astrometrica
# settings (PS1.cfg) write.
_BUILTIN = {
    "PS1 · gpc1": Profile("PS1 · gpc1", "F51", "J. Bulger, T. Lowe, A. Schultz, M. Willman",
                          "1.8-m f/4.4 Ritchey-Chretien + CCD", builtin=True),
}


def profile_key(site: dict[str, Any]) -> str:
    """The telescope's name as the images give it: telescope · instrument, else origin or observatory."""
    names = [site.get("telescope") or "", site.get("instrument") or ""]
    key = " · ".join(n for n in names if n)
    return key or site.get("origin") or site.get("observatory") or "(no telescope name in the images)"


def _load() -> dict[str, Profile]:
    try:
        data = json.loads(PROFILES_FILE.read_text())
    except (OSError, ValueError):
        return {}
    return {k: Profile(**{f: v for f, v in p.items() if f in Profile.__dataclass_fields__}) for k, p in data.items()}


def all_profiles() -> list[Profile]:
    """Saved profiles, then built-in ones not overridden."""
    saved = _load()
    return list(saved.values()) + [p for k, p in _BUILTIN.items() if k not in saved]


def find_profile(site: dict[str, Any]) -> Profile | None:
    key = profile_key(site)
    return _load().get(key) or _BUILTIN.get(key)


def save_profile(profile: Profile) -> None:
    profiles = _load()
    profile.builtin = False
    profiles[profile.key] = profile
    PROFILES_FILE.parent.mkdir(parents=True, exist_ok=True)
    PROFILES_FILE.write_text(json.dumps({k: asdict(p) for k, p in profiles.items()}, indent=2, ensure_ascii=False))


def delete_profile(key: str) -> None:
    profiles = _load()
    if profiles.pop(key, None) is not None:
        PROFILES_FILE.write_text(json.dumps({k: asdict(p) for k, p in profiles.items()}, indent=2, ensure_ascii=False))


def report_header(config, profile: Profile) -> dict[str, str]:
    """The report header's lines: the telescope's from its profile, the measurers' from Settings."""
    return {"observers": profile.observers, "telescope": profile.telescope,
            "measurers": config.measurers, "contact": config.contact}


def save_profile_used(entry, profile: Profile) -> None:
    """Record in the set's saved result which profile it was analyzed with (the report needs it later)."""
    path = entry.cache_dir / "result.json"
    entry.result["profile"] = asdict(profile)
    path.write_text(json.dumps(entry.result))


def profile_of(result: dict[str, Any]) -> Profile | None:
    data = result.get("profile")
    return Profile(**{f: v for f, v in data.items() if f in Profile.__dataclass_fields__}) if data else None
