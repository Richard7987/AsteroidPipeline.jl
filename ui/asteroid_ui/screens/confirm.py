"""A yes or no before something that can't be undone."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, Static


class ConfirmScreen(ModalScreen[bool]):
    BINDINGS = [Binding("y", "yes", "Yes"), Binding("n,escape", "no", "No")]

    def __init__(self, title: str, question: str) -> None:
        super().__init__()
        self.title_text = title
        self.question = question

    def compose(self) -> ComposeResult:
        with Vertical(id="reason-dialog", classes="dialog"):
            yield Static(self.title_text, classes="dialog-title", markup=False)
            yield Static(self.question, classes="dialog-help", markup=False)
            with Horizontal(classes="buttons"):
                yield Button("Yes (y)", id="yes", variant="primary")
                yield Button("No (n)", id="no")

    def action_yes(self) -> None:
        self.dismiss(True)

    def action_no(self) -> None:
        self.dismiss(False)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.dismiss(event.button.id == "yes")
