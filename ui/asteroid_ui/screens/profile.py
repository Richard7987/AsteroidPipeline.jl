"""Confirm or edit a telescope's profile: its MPC code and its report header lines."""

from __future__ import annotations

from typing import Any

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import Button, Input, Label, Select, Static

from asteroid_ui.profiles import Profile, profile_key, save_profile

_OTHER = "__other__"


def describe_site(site: dict[str, Any]) -> str:
    parts = []
    for key, label in (("telescope", "telescope"), ("instrument", "instrument"), ("origin", "origin"),
                       ("observatory", "observatory")):
        if site.get(key):
            parts.append(f"{label} «{site[key]}»")
    if site.get("lat") is not None and site.get("lon") is not None:
        lat, lon = site["lat"], site["lon"]
        where = f"{abs(lat):.3f}° {'N' if lat >= 0 else 'S'}, {abs(lon):.3f}° {'E' if lon >= 0 else 'W'}"
        if site.get("elev") is not None:
            where += f", {site['elev']:.0f} m"
        parts.append(f"site {where}")
    return " · ".join(parts) if parts else "nothing about where they were taken"


class ProfileScreen(ModalScreen[Profile | None]):
    BINDINGS = [Binding("ctrl+s", "save", "Save"), Binding("escape", "cancel", "Cancel")]

    def __init__(self, site: dict[str, Any], current: Profile | None = None) -> None:
        super().__init__()
        self.site = site
        self.key = current.key if current else profile_key(site)
        self.current = current

    def compose(self) -> ComposeResult:
        candidates = self.site.get("candidates") or []
        options = []
        for c in candidates:
            where = f", {c['distance_km']:.1f} km from the site in the images" if c.get("distance_km") is not None else ""
            options.append((f"{c['code']} — {c['name']}{where}", c["code"]))
        options.append(("Another code (type it below)", _OTHER))
        codes = [c["code"] for c in candidates]
        chosen = self.current.observatory if self.current else (codes[0] if codes else _OTHER)
        if chosen not in codes:
            chosen = _OTHER
        manual = self.current.observatory if self.current and self.current.observatory not in codes else ""

        with VerticalScroll(id="profile-dialog", classes="dialog"):
            yield Static(f"Telescope profile: {self.key}", classes="dialog-title")
            yield Static(Text.assemble("The images say: ", (describe_site(self.site), "bold")), classes="dialog-help")
            yield Static("Confirm it once: every image set from this telescope will use it from now on. "
                         "You can change it later in Settings → Telescope profiles.", classes="dialog-help")
            yield Label("Observatory (MPC code) — the best match first", classes="form-label")
            yield Select(options, value=chosen, allow_blank=False, id="profile-code")
            yield Input(manual, placeholder="three characters, e.g. F51", id="profile-code-manual",
                        classes="" if chosen == _OTHER else "hidden")
            yield Label("Report header, OBS: who observed with this telescope", classes="form-label")
            yield Input(self.current.observers if self.current else "", placeholder="e.g. J. Bulger, T. Lowe",
                        id="profile-observers")
            yield Label("Report header, TEL: the telescope", classes="form-label")
            yield Input(self.current.telescope if self.current else self.site.get("telescope", ""),
                        placeholder="e.g. 1.8-m f/4.4 Ritchey-Chretien + CCD", id="profile-telescope")
            with Horizontal(classes="buttons"):
                yield Button("Save (Ctrl+S)", id="save", variant="primary")
                yield Button("Cancel (Esc)", id="cancel")

    def on_mount(self) -> None:
        self.query_one("#profile-code", Select).focus()

    def on_select_changed(self, event: Select.Changed) -> None:
        manual = self.query_one("#profile-code-manual", Input)
        manual.set_class(event.value != _OTHER, "hidden")
        if event.value == _OTHER:
            manual.focus()

    def action_save(self) -> None:
        value = self.query_one("#profile-code", Select).value
        code = (self.query_one("#profile-code-manual", Input).value if value == _OTHER else str(value)).strip().upper()
        if len(code) != 3:
            self.app.notify("The observatory code has exactly 3 characters (e.g. F51).", severity="error")
            return
        profile = Profile(self.key, code, self.query_one("#profile-observers", Input).value.strip(),
                          self.query_one("#profile-telescope", Input).value.strip())
        save_profile(profile)
        self.dismiss(profile)

    def action_cancel(self) -> None:
        self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        {"save": self.action_save, "cancel": self.action_cancel}[event.button.id]()
