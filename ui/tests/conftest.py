"""Every test gets its own config, cache and data folders: never the person's."""

import os
import tempfile
from pathlib import Path

_base = Path(tempfile.mkdtemp(prefix="asteroid-ui-tests"))
for _name in ("CONFIG", "CACHE", "DATA"):
    os.environ[f"XDG_{_name}_HOME"] = str(_base / _name.lower())
