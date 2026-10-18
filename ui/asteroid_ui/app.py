from __future__ import annotations

from textual.app import App
from textual.reactive import reactive

# textual-image finds out which graphics the terminal supports (Kitty, Sixel...)
# by asking it with escape sequences; that only works before Textual takes the
# terminal over, so it is imported here, at start-up.
try:
    import textual_image.widget  # noqa: F401
except ImportError:
    pass

from asteroid_ui.backend import Backend
from asteroid_ui.config import load_config
from asteroid_ui.texts import STATUS_MESSAGES
from asteroid_ui.theme import GRUVBOX


class AsteroidApp(App):
    CSS_PATH = "styles.tcss"
    TITLE = "asteroids"
    SUB_TITLE = "AsteroidPipeline"

    engine_status: reactive[str] = reactive("Starting the analysis engine (Julia)… this can take a minute.")

    def __init__(self) -> None:
        super().__init__()
        self.config = load_config()
        self.backend = Backend()
        self.backend.on_status = lambda code: setattr(self, "engine_status", STATUS_MESSAGES.get(code, code))

    def on_mount(self) -> None:
        self.register_theme(GRUVBOX)
        self.theme = "gruvbox-asteroides"
        from asteroid_ui.screens.home import HomeScreen

        self.push_screen(HomeScreen())
        self.run_worker(self._start_backend())

    async def _start_backend(self) -> None:
        await self.backend.start()
        await self.backend.ready.wait()
        self.engine_status = f"✗ {self.backend.failed}" if self.backend.failed else "Analysis engine ready."

    async def on_unmount(self) -> None:
        await self.backend.stop()


USAGE = """\
usage: asteroidpipeline [--setup | --version | --help]

Find moving objects in astronomical image sets (.zip files, folders or FITS
files), check the candidates, and write MPC reports.

  --setup     install the analysis engine (Julia packages) and exit;
              the app otherwise does it on its first run
  --version   print the version and exit
"""


def main() -> None:
    import sys

    args = sys.argv[1:]
    if "--help" in args or "-h" in args:
        print(USAGE)
        return
    if "--version" in args:
        from asteroid_ui import __version__

        print(f"asteroidpipeline {__version__}")
        return
    if "--setup" in args:
        from asteroid_ui import engine

        try:
            julia = engine.julia_executable()
            print(f"Setting up the analysis engine with {julia} (several minutes)…")
            engine.setup(julia, print)
        except engine.EngineError as e:
            sys.exit(str(e))
        print("Done.")
        return
    AsteroidApp().run()


if __name__ == "__main__":
    main()
