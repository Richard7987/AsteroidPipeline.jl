"""Check an MPC report file before submitting it — from this program or any other."""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import Button, Footer, Input, Label, Select, Static, TextArea

from asteroid_ui.backend import BackendError
from asteroid_ui.clipboard import copy_text
from asteroid_ui.config import SETS_CACHE
from asteroid_ui.sets import SetEntry, discover, load_result
from asteroid_ui.texts import ICON_ERROR, ICON_OK, ICON_WARNING, ISSUES
from asteroid_ui.widgets import StatusBar


class CheckScreen(Screen):
    BINDINGS = [
        Binding("ctrl+l", "reload", "Check again"),
        Binding("ctrl+y", "copy", "Copy corrected report"),
        Binding("escape", "back", "Back"),
    ]

    def __init__(self, entry: SetEntry | None) -> None:
        super().__init__()
        self.entry = entry
        self.cleaned = ""

    def compose(self) -> ComposeResult:
        with Vertical(id="check"):
            yield Static("Check a report before submitting it", classes="screen-title")
            yield Label("Report file (MPC 80-column format):", classes="form-label")
            yield Input(self.app.config.report_to_check, placeholder="/path/to/MPCReport.txt", id="check-path")
            yield Select([], prompt="Which image set is this report for? (enables more checks)", id="check-set")
            yield Static("", id="issues")
            with Horizontal(id="check-panes"):
                with Vertical():
                    yield Label("Your report", classes="section-title")
                    yield TextArea("", read_only=True, soft_wrap=False, id="original")
                with Vertical():
                    yield Label("Corrected report", classes="section-title")
                    yield TextArea("", read_only=True, soft_wrap=False, id="cleaned")
            with Horizontal(classes="buttons"):
                yield Button("Check (Enter / Ctrl+L)", id="reload", variant="primary")
                yield Button("Copy corrected report (Ctrl+Y)", id="copy")
                yield Button("Back (Esc)", id="back")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        analyzed = [e for e in discover(self.app.config.downloads) if e.analyzed]
        select = self.query_one("#check-set", Select)
        select.set_options([(f"{e.name}  ({e.date})", e.name) for e in analyzed])
        if self.entry is not None and any(e.name == self.entry.name for e in analyzed):
            select.value = self.entry.name
        path_input = self.query_one("#check-path", Input)
        path_input.focus()
        if path_input.value:
            self.action_reload()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.action_reload()

    def on_select_changed(self, event: Select.Changed) -> None:
        if self.query_one("#check-path", Input).value.strip():
            self.action_reload()

    def action_reload(self) -> None:
        self.run_worker(self._check(), exclusive=True)

    async def _check(self) -> None:
        raw = self.query_one("#check-path", Input).value.strip()
        issues = self.query_one("#issues", Static)
        path = Path(raw).expanduser()
        if not raw or not path.is_file():
            issues.update(Text(f"{ICON_ERROR} There is no file at «{raw}».", style="bold red"))
            return
        self.app.config.report_to_check = str(path)
        self.app.config.save()
        changed = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        text = path.read_text(encoding="latin-1")
        self.query_one("#original", TextArea).text = text

        value = self.query_one("#check-set", Select).value
        result = str(SETS_CACHE / value / "result.json") if isinstance(value, str) else None
        # the set's observatory; without a set, the report's own COD line
        observatory = (load_result(value) or {}).get("observatory") if isinstance(value, str) else None
        if not observatory:
            cod = re.search(r"^COD +(\S{3})", text, re.MULTILINE)
            observatory = cod[1] if cod else "XXX"
        issues.update("Checking…")
        try:
            answer = await self.app.backend.request(
                "check_report", text=text, observatory=observatory, result=result)
        except BackendError as e:
            issues.update(Text(f"{ICON_ERROR} The check failed: {e}", style="bold red"))
            return
        self.cleaned = answer["cleaned"]
        self.query_one("#cleaned", TextArea).text = self.cleaned
        lines = Text(f"File last changed {changed}.\n", style="dim")
        if not answer["issues"]:
            lines.append(f"{ICON_OK} The report is fine: you can submit it as it is.", style="bold green")
        else:
            verdict = ("The report has errors: submit the corrected one (on the right) instead."
                       if not answer["ok"] else "The report can be submitted, but check these points:")
            lines.append(f"{verdict}\n", style="bold")
            for issue in answer["issues"]:
                is_error, template = ISSUES.get(issue["code"], (True, issue["code"] + " {detail}"))
                icon, style = (ICON_ERROR, "red") if is_error else (ICON_WARNING, "yellow")
                where = f"line {issue['line']}: " if issue["line"] else ""
                lines.append(f" {icon} {where}", style=f"bold {style}")
                lines.append(template.format(detail=issue["detail"].strip()) + "\n")
            if result is None:
                lines.append("Choose the image set above to also catch other sets' lines and known asteroids.",
                             style="dim italic")
        issues.update(lines)

    def action_copy(self) -> None:
        if self.cleaned:
            copy_text(self.app, self.cleaned)
            self.app.notify("Corrected report copied.")

    def action_back(self) -> None:
        self.app.pop_screen()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        {"reload": self.action_reload, "copy": self.action_copy, "back": self.action_back}[event.button.id]()
