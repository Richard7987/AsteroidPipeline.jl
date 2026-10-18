"""Analysis in progress: one or more image sets, step by step."""

from __future__ import annotations

import time

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import Screen
from textual.widgets import Footer, ProgressBar, Static

from asteroid_ui.backend import BackendError
from asteroid_ui.labels import add_catalog_labels
from asteroid_ui.profiles import Profile, find_profile, report_header, save_profile_used
from asteroid_ui.screens.profile import ProfileScreen
from asteroid_ui.sets import SetEntry, load_result
from asteroid_ui.texts import ICON_DONE, ICON_ERROR, ICON_PENDING, ICON_RUNNING, STAGES
from asteroid_ui.widgets import StatusBar

_ORDER = [code for code, _ in STAGES]


class AnalysisScreen(Screen):
    BINDINGS = [Binding("escape", "back", "Back")]

    def __init__(self, entries: list[SetEntry]) -> None:
        super().__init__()
        self.entries = entries
        self.errors: dict[str, str] = {}
        self.stage: str | None = None
        self.started = time.monotonic()
        self.finished = False

    def compose(self) -> ComposeResult:
        with Vertical(id="analysis"):
            yield Static("", id="analysis-title", classes="screen-title")
            yield Static("", id="analysis-source")
            yield Static("", id="analysis-profile")
            yield Static("", id="steps")
            yield ProgressBar(total=100, show_eta=False, id="progress")
            yield Static("", id="elapsed")
            yield Static("", id="analysis-error")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self.set_interval(1, self._tick)
        self.run_worker(self._run_all(), exclusive=True)

    def _tick(self) -> None:
        if not self.finished:
            s = int(time.monotonic() - self.started)
            self.query_one("#elapsed", Static).update(
                f"Elapsed: {s // 60}:{s % 60:02d}   (a set usually takes 1–3 minutes; the first one, longer)")

    def _show_steps(self, failed: bool = False) -> None:
        current = _ORDER.index(self.stage) if self.stage in _ORDER else -1
        lines = Text()
        for i, (_, label) in enumerate(STAGES):
            if i < current:
                lines.append(f" {ICON_DONE}  {label}\n", style="green")
            elif i == current:
                icon, style = (ICON_ERROR, "bold red") if failed else (ICON_RUNNING, "bold yellow")
                lines.append(f" {icon}  {label}\n", style=style)
            else:
                lines.append(f" {ICON_PENDING}  {label}\n", style="dim")
        self.query_one("#steps", Static).update(lines)

    async def _on_event(self, event: dict) -> None:
        if event.get("event") == "progress":
            stage = event.get("stage")
            if stage == "done":
                stage = "images"
            if stage in _ORDER:
                self.stage = stage
                self._show_steps()
            self.query_one("#progress", ProgressBar).update(progress=100 * float(event.get("fraction", 0)))

    async def _run_all(self) -> None:
        app = self.app
        for i, entry in enumerate(self.entries, 1):
            counter = f" — set {i} of {len(self.entries)}" if len(self.entries) > 1 else ""
            self.query_one("#analysis-title", Static).update(f"Analyzing {entry.name}{counter}")
            self.query_one("#analysis-source", Static).update(f"From: {entry.source_label}")
            self.stage = "extract"
            self._show_steps()
            self.query_one("#progress", ProgressBar).update(progress=0)
            app.engine_status = f"Analyzing {entry.name}…"
            source = ([str(p) for p in entry.source] if isinstance(entry.source, list) else str(entry.source))
            try:
                profile = await self._profile_for(source)
                if profile is None:
                    raise BackendError("no observatory was chosen for these images")
                await app.backend.request(
                    "analyze", on_event=self._on_event, source=source, out=str(entry.cache_dir),
                    name=entry.name, observatory=profile.observatory, header=report_header(app.config, profile))
                entry.result = load_result(entry.name)
                if entry.result is not None:
                    save_profile_used(entry, profile)
                    add_catalog_labels(entry.name, entry.result, entry.cache_dir)
            except BackendError as e:
                self.errors[entry.name] = str(e)
                self._show_steps(failed=True)
        self.finished = True
        app.engine_status = "Analysis engine ready."
        done = [e for e in self.entries if e.result is not None]
        if len(self.entries) == 1 and done:
            from asteroid_ui.screens.results import ResultsScreen

            app.switch_screen(ResultsScreen(done[0]))
        elif len(self.entries) > 1:
            from asteroid_ui.screens.summary import SummaryScreen

            app.switch_screen(SummaryScreen(self.entries, self.errors))
        else:
            message = next(iter(self.errors.values()), "unknown error")
            self.query_one("#analysis-error", Static).update(
                Text.assemble((f"{ICON_ERROR} The analysis failed: ", "bold red"), message,
                              "\n\nPress Esc to go back. The engine's log has the details."))

    async def _profile_for(self, source: str | list[str]) -> Profile | None:
        """The telescope's profile, from the images' headers — asking the first time."""
        self.app.engine_status = "Reading the images' headers…"
        site = await self.app.backend.request("inspect", source=source)
        profile = find_profile(site)
        if profile is None:
            profile = await self.app.push_screen_wait(ProfileScreen(site))
        if profile is not None:
            self.query_one("#analysis-profile", Static).update(f"Observatory: {profile.describe()}")
        return profile

    def action_back(self) -> None:
        if self.finished:
            self.app.pop_screen()
        else:
            self.app.notify("The analysis is still running; it can't be stopped halfway.", severity="warning")
