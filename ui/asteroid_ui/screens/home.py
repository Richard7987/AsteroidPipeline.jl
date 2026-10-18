"""Start screen: the image sets found, and how to open others."""

from __future__ import annotations

from pathlib import Path

import pyfiglet
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import Footer, Label, OptionList, Static
from textual.widgets.option_list import Option

from asteroid_ui.labels import labels_for
from asteroid_ui.sets import SetEntry, discover, entry_for, summary
from asteroid_ui.texts import HELP_HOME
from asteroid_ui.widgets import StatusBar



def banner(title: str) -> str:
    """The title in large letters, if it fits a terminal; as it is, otherwise."""
    text = pyfiglet.Figlet(font="small", width=200).renderText(title).rstrip("\n")
    return text if max(map(len, text.splitlines()), default=0) <= 120 else title


def start_analysis(screen: Screen, entries: list[SetEntry]) -> None:
    """Open the analysis screen for `entries`, if everything it needs is set."""
    app = screen.app
    if not entries:
        app.notify("There is nothing new to analyze.")
        return
    if not app.config.measurers.strip():
        app.notify("First, say who measures (Settings): it goes in every report.", severity="error", timeout=8)
        from asteroid_ui.screens.settings import SettingsScreen

        app.push_screen(SettingsScreen())
        return
    if app.backend.busy:
        app.notify("An analysis is already running; wait for it to finish.", severity="warning")
        return
    from asteroid_ui.screens.analysis import AnalysisScreen

    app.push_screen(AnalysisScreen(entries))


class HomeScreen(Screen):
    BINDINGS = [
        Binding("o", "open", "Open files"),
        Binding("t", "analyze_all", "Analyze all new"),
        Binding("a", "reanalyze", "Analyze again"),
        Binding("v", "check_report", "Check a report"),
        Binding("s", "settings", "Settings"),
        Binding("u", "refresh", "Refresh"),
        Binding("question_mark", "help", "Help"),
        Binding("q", "app.quit", "Quit"),
    ]

    def __init__(self) -> None:
        super().__init__()
        self.entries: list[SetEntry] = []

    def compose(self) -> ComposeResult:
        with Vertical(id="home"):
            yield Static(banner(self.app.config.title), id="banner")
            yield Static("find moving objects in your images · write the MPC report", id="subtitle")
            with Horizontal(id="home-body"):
                with Vertical(id="sets-pane"):
                    yield Label("", id="sets-title", classes="section-title")
                    yield OptionList(id="sets")
                yield Static(HELP_HOME, id="help-pane")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self.refresh_sets()

    def on_screen_resume(self) -> None:
        self.query_one("#banner", Static).update(banner(self.app.config.title))
        self.refresh_sets()

    def refresh_sets(self) -> None:
        folder = self.app.config.downloads
        self.entries = discover(folder)
        self.query_one("#sets-title", Label).update(f"Image sets — .zip files in {folder}, and earlier analyses")
        options = self.query_one("#sets", OptionList)
        highlighted = options.highlighted
        options.clear_options()
        if not self.entries:
            options.add_option(Option(Text(f"No image sets found in {folder}. Press o to open files.",
                                           style="italic"), disabled=True))
            return
        for i, e in enumerate(self.entries):
            if e.analyzed:
                state = Text("✓ ", style="bold green")
                decided = {t for t, r in labels_for(e.name).items() if r["source"] == "human"}
                note = Text(summary(e.result, decided))
            else:
                state = Text("● ", style="bold yellow")
                note = Text("not analyzed yet — Enter to analyze", style="italic")
            label = Text.assemble(state, (f"{e.name:<24}", "bold"), (f" {e.date:<11} ", "dim"), note)
            options.add_option(Option(label, id=str(i)))
        options.highlighted = min(highlighted or 0, len(self.entries) - 1)
        options.focus()

    def _highlighted(self) -> SetEntry | None:
        options = self.query_one("#sets", OptionList)
        if options.highlighted is None or not self.entries:
            return None
        option = options.get_option_at_index(options.highlighted)
        return self.entries[int(option.id)] if option.id is not None else None

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        if event.option.id is None:
            return
        entry = self.entries[int(event.option.id)]
        if entry.analyzed:
            from asteroid_ui.screens.results import ResultsScreen

            self.app.push_screen(ResultsScreen(entry))
        else:
            start_analysis(self, [entry])

    def action_reanalyze(self) -> None:
        entry = self._highlighted()
        if entry is None:
            return
        if entry.source is None:
            self.app.notify("The images of this set are no longer where they were.", severity="error")
            return
        start_analysis(self, [entry])

    def action_analyze_all(self) -> None:
        start_analysis(self, [e for e in self.entries if not e.analyzed and e.source is not None])

    def action_open(self) -> None:
        from asteroid_ui.screens.open import OpenScreen

        def opened(source: Path | list[Path] | None) -> None:
            if source is not None:
                start_analysis(self, [entry_for(source)])

        self.app.push_screen(OpenScreen(self.app.config.downloads), opened)

    def action_check_report(self) -> None:
        from asteroid_ui.screens.check import CheckScreen

        entry = self._highlighted()
        self.app.push_screen(CheckScreen(entry if entry and entry.analyzed else None))

    def action_settings(self) -> None:
        from asteroid_ui.screens.settings import SettingsScreen

        self.app.push_screen(SettingsScreen())

    def action_refresh(self) -> None:
        self.refresh_sets()

    def action_help(self) -> None:
        from asteroid_ui.screens.help import HelpScreen

        self.app.push_screen(HelpScreen())
