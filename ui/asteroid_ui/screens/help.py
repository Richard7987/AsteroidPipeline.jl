"""Keys and marks, at a glance."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, Static

from asteroid_ui.texts import HELP_HOME, HELP_KEYS


class HelpScreen(ModalScreen[None]):
    BINDINGS = [Binding("escape,question_mark,q", "close", "Close")]

    def compose(self) -> ComposeResult:
        with Vertical(id="help-dialog", classes="dialog"):
            yield Static("Help", classes="dialog-title")
            yield Static(HELP_HOME)
            yield Static(HELP_KEYS)
            yield Button("Close (Esc)", id="close", variant="primary")

    def action_close(self) -> None:
        self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.action_close()
