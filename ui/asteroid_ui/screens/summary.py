"""After analyzing several image sets at once: one row per set."""

from __future__ import annotations

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import Screen
from textual.widgets import DataTable, Footer, Static

from asteroid_ui.sets import SetEntry
from asteroid_ui.texts import ICON_ERROR, VERDICTS
from asteroid_ui.widgets import StatusBar


class SummaryScreen(Screen):
    BINDINGS = [Binding("escape", "back", "Back")]

    def __init__(self, entries: list[SetEntry], errors: dict[str, str]) -> None:
        super().__init__()
        self.entries = entries
        self.errors = errors

    def compose(self) -> ComposeResult:
        with Vertical(id="summary"):
            yield Static(f"{len(self.entries)} image sets analyzed", classes="screen-title")
            yield Static("Enter opens a set: its objects, its images and its report.")
            yield DataTable(cursor_type="row", zebra_stripes=True, id="summary-table")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        table = self.query_one(DataTable)
        table.add_columns("Set", "Date", "Result", "New", "Known", "Images used")
        for e in self.entries:
            if e.result is None:
                table.add_row(e.name, e.date, Text(f"{ICON_ERROR} failed: {self.errors.get(e.name, '')[:60]}",
                                                   style="red"), "", "", "", key=e.name)
                continue
            r = e.result
            icon, title, style, _ = VERDICTS[r["report_kind"]]
            colour = {"success": "green", "warning": "yellow", "error": "red"}[style]
            new = sum(1 for t in r["tracklets"] if not t.get("known"))
            known = sum(1 for k in r["known"] if k["detected"])
            used = sum(1 for f in r["frames"] if f["used"])
            table.add_row(e.name, e.date, Text(f"{icon} {title}", style=colour), str(new), str(known),
                          f"{used} of {len(r['frames'])}", key=e.name)
        table.focus()

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        entry = next(e for e in self.entries if e.name == event.row_key.value)
        if entry.result is not None:
            from asteroid_ui.screens.results import ResultsScreen

            self.app.push_screen(ResultsScreen(entry))

    def action_back(self) -> None:
        self.app.pop_screen()
