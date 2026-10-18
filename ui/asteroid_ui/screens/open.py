"""Choose what to analyze: a .zip, a folder of FITS images, or FITS files."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from rich.style import Style
from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.message import Message
from textual.screen import ModalScreen
from textual.widgets import Button, DirectoryTree, Footer, Static
from textual.widgets._directory_tree import DirEntry
from textual.widgets._tree import TreeNode

from asteroid_ui.sets import FITS_SUFFIXES, is_fits

INSTRUCTIONS = (
    "[b]Enter[/b] on a [b].zip[/b] opens it.   "
    "[b]f[/b] on a folder opens the FITS images in it.\n"
    "[b]Space[/b] on FITS files marks them one by one; then [b]a[/b] analyzes the marked files (3 or more).   "
    "[b]Backspace[/b] goes up one folder."
)


class SetTree(DirectoryTree):
    """Folders, .zip files and FITS images only; FITS images can be marked."""

    BINDINGS = [
        Binding("space", "toggle_mark", "Mark / unmark", show=True),
        Binding("f", "use_folder", "Open this folder", show=True),
        Binding("backspace", "go_up", "Up one folder", show=True),
    ]

    def __init__(self, path: Path, **kwargs) -> None:
        super().__init__(path, **kwargs)
        self.marked: list[Path] = []

    def filter_paths(self, paths: Iterable[Path]) -> Iterable[Path]:
        return [p for p in paths if not p.name.startswith(".")
                and (p.is_dir() or p.suffix.lower() == ".zip" or p.suffix.lower() in FITS_SUFFIXES)]

    def render_label(self, node: TreeNode[DirEntry], base_style: Style, style: Style) -> Text:
        label = super().render_label(node, base_style, style)
        if node.data is not None and is_fits(node.data.path):
            mark = "☑ " if node.data.path in self.marked else "☐ "
            label = Text.assemble((mark, "bold"), label)
        return label

    def action_toggle_mark(self) -> None:
        node = self.cursor_node
        if node is None or node.data is None:
            return
        path = node.data.path
        if is_fits(path):
            if path in self.marked:
                self.marked.remove(path)
            else:
                self.marked.append(path)
            self.refresh()
            self.post_message(SetTree.MarksChanged())
            self.action_cursor_down()
        else:
            self.action_toggle_node()

    def action_use_folder(self) -> None:
        node = self.cursor_node
        if node is not None and node.data is not None:
            folder = node.data.path if node.data.path.is_dir() else node.data.path.parent
            self.screen.use_folder(folder)  # type: ignore[attr-defined]

    def action_go_up(self) -> None:
        parent = Path(self.path).resolve().parent
        if parent != Path(self.path).resolve():
            self.path = parent

    class MarksChanged(Message):
        pass


class OpenScreen(ModalScreen[Path | list[Path] | None]):
    BINDINGS = [
        Binding("a", "analyze_marked", "Analyze marked files"),
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self, start: str) -> None:
        super().__init__()
        start_path = Path(start).expanduser()
        self.start = start_path if start_path.is_dir() else Path.home()

    def compose(self) -> ComposeResult:
        with Vertical(id="open-dialog", classes="dialog"):
            yield Static("Open an image set", classes="dialog-title")
            yield Static(INSTRUCTIONS, id="open-help")
            yield SetTree(self.start, id="open-tree")
            yield Static("No files marked.", id="open-marked")
            with Horizontal(classes="buttons"):
                yield Button("Analyze marked files (a)", id="analyze-marked", variant="primary")
                yield Button("Cancel (Esc)", id="cancel")
        yield Footer()

    def on_mount(self) -> None:
        self.query_one(SetTree).focus()

    def on_set_tree_marks_changed(self) -> None:
        marked = self.query_one(SetTree).marked
        text = (f"{len(marked)} file(s) marked: " + ", ".join(p.name for p in marked[:4])
                + (" …" if len(marked) > 4 else "")) if marked else "No files marked."
        self.query_one("#open-marked", Static).update(text)

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        if event.path.suffix.lower() == ".zip":
            self.dismiss(event.path)
        elif is_fits(event.path):
            self.query_one(SetTree).action_toggle_mark()

    def use_folder(self, folder: Path) -> None:
        n = sum(1 for p in folder.iterdir() if p.is_file() and is_fits(p))
        if n < 3:
            self.app.notify(f"{folder.name} has {n} FITS image(s); a set needs at least 3.", severity="warning")
            return
        self.dismiss(folder)

    def action_analyze_marked(self) -> None:
        marked = self.query_one(SetTree).marked
        if len(marked) < 3:
            self.app.notify("Mark at least 3 FITS images of the same field (Space).", severity="warning")
            return
        self.dismiss(sorted(marked))

    def action_cancel(self) -> None:
        self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "analyze-marked":
            self.action_analyze_marked()
        else:
            self.action_cancel()
