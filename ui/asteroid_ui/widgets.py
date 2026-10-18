"""Custom widgets: the status bar and the image viewer."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Container
from textual.reactive import reactive
from textual.widgets import Static

try:
    from textual_image.widget import Image as TerminalImage
    from textual_image.widget import get_cell_size
except ImportError:  # no textual-image: the viewer says so
    TerminalImage = None  # type: ignore[assignment,misc]

from asteroid_ui.imaging import SetImages

# Canvas resolution cap: more doesn't look better in a terminal and slows the
# animation down.
_MAX_CANVAS = (1800, 1300)


class StatusBar(Static):
    """One line at the bottom with the analysis engine's state."""

    def on_mount(self) -> None:
        self.watch(self.app, "engine_status", self._update)

    def _update(self, status: str) -> None:
        self.update(f"⚙  {status}")


class SkyViewer(Container):
    """One of the set's images, with marks, at exactly the panel's size."""

    frame: reactive[int] = reactive(1)
    marks: reactive[bool] = reactive(True)
    zoom_on: reactive[str | None] = reactive(None)
    selected: reactive[str | None] = reactive(None)

    def __init__(self, images: SetImages, **kwargs) -> None:
        super().__init__(**kwargs)
        self.images = images

    def compose(self) -> ComposeResult:
        if TerminalImage is None:
            yield Static("Images can't be shown: the textual-image package is missing.", id="no-image")
        else:
            yield TerminalImage(id="sky")

    def on_resize(self) -> None:
        self.redraw()

    def watch_frame(self) -> None:
        self.redraw()

    def watch_marks(self) -> None:
        self.redraw()

    def watch_zoom_on(self) -> None:
        self.redraw()

    def watch_selected(self) -> None:
        self.redraw()

    def redraw(self) -> None:
        if TerminalImage is None or not self.is_mounted:
            return
        cols, rows = self.content_size.width, self.content_size.height
        if cols < 4 or rows < 2:
            return
        cell = get_cell_size()
        w, h = cols * cell.width, rows * cell.height
        shrink = min(1.0, _MAX_CANVAS[0] / w, _MAX_CANVAS[1] / h)
        size = (int(w * shrink), int(h * shrink))
        canvas = self.images.render(self.frame, size, marks=self.marks, zoom_on=self.zoom_on, selected=self.selected)
        self.query_one("#sky", TerminalImage).image = canvas
