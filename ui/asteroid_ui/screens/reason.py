"""Why a candidate isn't real — one key per reason, so labelling stays quick."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import OptionList, Static
from textual.widgets.option_list import Option

from asteroid_ui.labels import REASONS


class ReasonScreen(ModalScreen[str | None]):
    BINDINGS = [Binding("escape", "cancel", "Cancel")] + [
        Binding(str(i), f"pick('{key}')", show=False) for i, key in enumerate(REASONS, 1)]

    def __init__(self, candidate: str) -> None:
        super().__init__()
        self.candidate = candidate

    def compose(self) -> ComposeResult:
        with Vertical(id="reason-dialog", classes="dialog"):
            yield Static(f"Why is {self.candidate} not real?", classes="dialog-title")
            yield OptionList(*[Option(f"{i}. {text}", id=key) for i, (key, text) in enumerate(REASONS.items(), 1)],
                             id="reasons")
            yield Static("Enter or the number chooses · Esc cancels", classes="dialog-help")

    def on_mount(self) -> None:
        self.query_one(OptionList).focus()

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        self.dismiss(event.option.id)

    def action_pick(self, key: str) -> None:
        self.dismiss(key)

    def action_cancel(self) -> None:
        self.dismiss(None)
