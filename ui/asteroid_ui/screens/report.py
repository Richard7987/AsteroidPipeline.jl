"""The MPC report of an image set — the candidates confirmed as real, or blank."""

from __future__ import annotations

from pathlib import Path

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, Footer, Static, TextArea

from asteroid_ui.backend import BackendError
from asteroid_ui.clipboard import copy_text
from asteroid_ui.labels import REAL
from asteroid_ui.profiles import Profile, profile_of, report_header
from asteroid_ui.sets import SetEntry
from asteroid_ui.texts import ICON_ERROR, ICON_OK, ICON_WARNING, REPORT_HELP, REPORT_TITLES


def report_path(entry: SetEntry) -> Path:
    """Where a report is saved: next to the .zip or the images it came from."""
    source = entry.source
    if isinstance(source, list):
        folder = source[0].parent
    elif source is not None:
        folder = source.parent if source.is_file() else source
    else:
        folder = Path((entry.result or {}).get("set_dir", Path.home()))
    return folder / f"{entry.name}_MPCReport.txt"


class ReportScreen(ModalScreen[None]):
    BINDINGS = [
        Binding("c", "copy", "Copy"),
        Binding("s", "save", "Save to a file"),
        Binding("escape", "close", "Close"),
    ]

    def __init__(self, entry: SetEntry, decisions: dict[int, str], reasons: dict[int, str] | None = None) -> None:
        super().__init__()
        self.entry = entry
        self.result = entry.result or {}
        candidates = [t["id"] for t in self.result.get("tracklets", []) if not t.get("known")]
        reasons = reasons or {}
        # real and new: a candidate recognised as a known object isn't reported
        self.confirmed = [i for i in candidates if decisions.get(i) == REAL
                          and not reasons.get(i, "").startswith("known:")]
        self.undecided = [i for i in candidates if i not in decisions]
        if self.confirmed:
            self.kind = "confirmed"
        else:
            self.kind = "unreliable" if self.result.get("report_kind") == "unreliable" else "blank"
        self.report = ""
        self.profile = profile_of(self.result) or Profile("", self.result.get("observatory", ""))

    def compose(self) -> ComposeResult:
        icon = ICON_ERROR if self.kind == "unreliable" else ICON_OK
        with Vertical(id="report-dialog", classes="dialog"):
            yield Static(f"{icon} {REPORT_TITLES[self.kind]} — {self.entry.name}", classes="dialog-title")
            yield Static(REPORT_HELP[self.kind], classes="dialog-help")
            if self.undecided:
                yield Static(Text(f"{ICON_WARNING} {len(self.undecided)} candidate(s) not checked yet are left out. "
                                  "Check them first (y / n) if they might be real.", style="bold yellow"),
                             classes="dialog-help")
            yield TextArea("Writing the report…", read_only=True, soft_wrap=False, id="report-text")
            with Horizontal(classes="buttons"):
                yield Button("Copy (c)", id="copy", variant="primary")
                yield Button("Save to a file (s)", id="save")
                yield Button("Close (Esc)", id="close")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one("#copy", Button).focus()
        self.run_worker(self._build())

    async def _build(self) -> None:
        from asteroid_ui.designations import assign

        try:
            confirmed = [t for t in self.result.get("tracklets", []) if t["id"] in self.confirmed]
            designations = assign(self.app.config, self.entry.name, confirmed)
        except ValueError as e:
            self.query_one("#report-text", TextArea).text = f"The designations could not be given: {e}"
            return
        try:
            answer = await self.app.backend.request(
                "build_report", result=str(self.entry.cache_dir / "result.json"),
                tracklets=self.confirmed, header=report_header(self.app.config, self.profile),
                designations=designations)
        except BackendError as e:
            self.query_one("#report-text", TextArea).text = f"The report could not be written: {e}"
            return
        self.report = answer["text"]
        self.query_one("#report-text", TextArea).text = self.report

    def action_copy(self) -> None:
        if self.report:
            copy_text(self.app, self.report)
            self.app.notify("Report copied. Paste it where you submit your reports.")

    def action_save(self) -> None:
        if not self.report:
            return
        path = report_path(self.entry)
        try:
            path.write_text(self.report)
        except OSError as e:
            self.app.notify(f"Could not save it: {e}", severity="error")
            return
        self.app.notify(f"Saved to {path}", timeout=8)

    def action_close(self) -> None:
        self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        {"copy": self.action_copy, "save": self.action_save, "close": self.action_close}[event.button.id]()
