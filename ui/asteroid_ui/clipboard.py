"""Copy text to the clipboard.

Textual's own copy (an OSC 52 escape sequence) works in most terminals, but not
all of them allow it; the desktop's clipboard tool is used as well when there
is one, so a copied report is really there to paste.
"""

from __future__ import annotations

import shutil
import subprocess

from textual.app import App

_TOOLS = (["wl-copy"], ["xclip", "-selection", "clipboard"], ["xsel", "--clipboard", "--input"])


def copy_text(app: App, text: str) -> None:
    app.copy_to_clipboard(text)
    for cmd in _TOOLS:
        if shutil.which(cmd[0]):
            try:
                subprocess.run(cmd, input=text.encode(), timeout=3, check=False,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except (OSError, subprocess.TimeoutExpired):
                continue
            break
