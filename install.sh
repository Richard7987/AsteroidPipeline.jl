#!/bin/sh
# Installs the AsteroidPipeline terminal app on any Linux (or macOS):
#
#   curl -fsSL https://raw.githubusercontent.com/Richard7987/AsteroidPipeline.jl/main/install.sh | sh
#
# or, from a downloaded release:  sh install.sh asteroidpipeline-0.2.1-py3-none-any.whl
#
# What it does: a private Python environment in ~/.local/share/asteroidpipeline,
# the app installed into it, the `asteroidpipeline` command in ~/.local/bin,
# and the analysis engine set up (Julia — downloaded if you don't have it —
# and its packages; several minutes, once). Needs Python 3.11 or newer.
# Remove everything with:  rm -rf ~/.local/share/asteroidpipeline ~/.local/bin/asteroidpipeline
set -eu

VERSION="0.2.1"
SOURCE="${1:-https://github.com/Richard7987/AsteroidPipeline.jl/releases/download/v$VERSION/asteroidpipeline-$VERSION-py3-none-any.whl}"
PREFIX="${XDG_DATA_HOME:-$HOME/.local/share}/asteroidpipeline"
BIN="$HOME/.local/bin"

say() { printf '\033[1m%s\033[0m\n' "$*"; }
fail() { printf 'error: %s\n' "$*" >&2; exit 1; }

PYTHON=""
for candidate in python3.13 python3.12 python3.11 python3; do
    if command -v "$candidate" >/dev/null 2>&1 &&
       "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
        PYTHON="$(command -v "$candidate")"
        break
    fi
done
[ -n "$PYTHON" ] || fail "Python 3.11 or newer is needed. Install it with your package manager
       (Debian/Ubuntu: sudo apt install python3 python3-venv; Fedora: sudo dnf install python3)."
"$PYTHON" -c 'import venv, ensurepip' 2>/dev/null ||
    fail "Python's venv module is missing (Debian/Ubuntu: sudo apt install python3-venv)."

say "Installing AsteroidPipeline $VERSION into $PREFIX"
"$PYTHON" -m venv --clear "$PREFIX"
"$PREFIX/bin/python" -m pip install --quiet --upgrade pip
"$PREFIX/bin/python" -m pip install --quiet "$SOURCE"

mkdir -p "$BIN"
ln -sf "$PREFIX/bin/asteroidpipeline" "$BIN/asteroidpipeline"

say "Setting up the analysis engine (Julia; several minutes, once)"
"$BIN/asteroidpipeline" --setup

say "Done. Start it with:  asteroidpipeline"
case ":$PATH:" in
    *":$BIN:"*) ;;
    *) echo "($BIN is not on your PATH: add it, or run $BIN/asteroidpipeline)" ;;
esac
