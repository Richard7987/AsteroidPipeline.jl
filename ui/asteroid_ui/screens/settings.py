"""Settings: where image sets are, who measures, and the telescope profiles."""

from __future__ import annotations

from pathlib import Path

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import Screen
from textual.widgets import Button, Footer, Input, Label, OptionList, Static
from textual.widgets.option_list import Option

from asteroid_ui.config import DATASET_DIR
from asteroid_ui.labels import dataset_counts
from asteroid_ui.profiles import all_profiles, delete_profile
from asteroid_ui.widgets import StatusBar

FIELDS = [
    ("title", "Title shown on the start screen — your team's name, for example", "NEO Hunters"),
    ("downloads", "Folder with your image-set .zip files (its subfolders are searched too)", ""),
    ("measurers", "Report header, MEA: who measures — you and your team (required)", "e.g. A. Name, Team Name"),
    ("contact", "Report header, CON: contact name and e-mail (optional)", "e.g. A. Name [name@example.org]"),
    ("fps", "Animation speed, images per second (0.5 to 3)", ""),
]


class SettingsScreen(Screen):
    BINDINGS = [Binding("ctrl+s", "save", "Save"), Binding("delete", "delete_profile", "Delete profile"),
                Binding("escape", "back", "Back")]

    def compose(self) -> ComposeResult:
        config = self.app.config
        with VerticalScroll(id="settings"):
            yield Static("Settings", classes="screen-title")
            for key, label, placeholder in FIELDS:
                with Vertical(classes="form-row"):
                    yield Label(label, classes="form-label")
                    yield Input(str(getattr(config, key)), placeholder=placeholder, id=f"set-{key}")
            yield Label("Telescope profiles — the observatory code and the OBS and TEL lines of every telescope "
                        "whose images you have opened. Enter edits one, Delete removes it.", classes="form-label")
            yield OptionList(id="profiles-list")
            counts = dataset_counts()
            yield Static(f"Training data so far ({DATASET_DIR}): {counts['real_human']} confirmed real, "
                         f"{counts['real_catalog']} known asteroids, {counts['bogus']} not real.",
                         id="settings-note")
            with Horizontal(classes="buttons"):
                yield Button("Save (Ctrl+S)", id="save", variant="primary")
                yield Button("Back (Esc)", id="back")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self._fill_profiles()
        self.query_one("#set-downloads", Input).focus()

    def _fill_profiles(self) -> None:
        options = self.query_one("#profiles-list", OptionList)
        options.clear_options()
        self.profiles = all_profiles()
        for p in self.profiles:
            options.add_option(Option(Text.assemble((f"{p.observatory}  ", "bold"), p.key,
                                                    (" (built in)" if p.builtin else "", "dim"),
                                                    (f"   TEL {p.telescope}" if p.telescope else "", "dim")), id=p.key))
        if not self.profiles:
            options.add_option(Option("No telescope profiles yet: they are made the first time you open "
                                      "a telescope's images.", disabled=True))

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        profile = next((p for p in self.profiles if p.key == event.option.id), None)
        if profile is None:
            return
        from asteroid_ui.screens.profile import ProfileScreen

        site = {"telescope": profile.key, "candidates": [{"code": profile.observatory, "name": "(current)"}]}
        self.app.push_screen(ProfileScreen(site, profile), lambda _: self._fill_profiles())

    def action_delete_profile(self) -> None:
        options = self.query_one("#profiles-list", OptionList)
        if not options.has_focus or options.highlighted is None:
            return
        key = options.get_option_at_index(options.highlighted).id
        profile = next((p for p in self.profiles if p.key == key), None)
        if profile is None or profile.builtin:
            self.app.notify("Built-in profiles can be edited, not removed.")
            return
        delete_profile(key)
        self.app.notify(f"Profile «{key}» removed: its images will ask again next time.")
        self._fill_profiles()

    def action_save(self) -> None:
        config = self.app.config
        values = {key: self.query_one(f"#set-{key}", Input).value.strip() for key, _, _ in FIELDS}
        try:
            fps = min(3.0, max(0.5, float(values["fps"].replace(",", "."))))
        except ValueError:
            self.app.notify("The animation speed must be a number between 0.5 and 3.", severity="error")
            return
        if not Path(values["downloads"]).expanduser().is_dir():
            self.app.notify(f"The folder {values['downloads']} doesn't exist.", severity="error")
            return
        if not values["measurers"]:
            self.app.notify("Say who measures: it goes in every report.", severity="error")
            return
        config.downloads, config.measurers, config.contact = values["downloads"], values["measurers"], values["contact"]
        config.title = values["title"] or "NEO Hunters"
        config.fps = fps
        config.save()
        self.app.notify("Settings saved.")
        self.app.pop_screen()

    def action_back(self) -> None:
        self.app.pop_screen()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        {"save": self.action_save, "back": self.action_back}[event.button.id]()
