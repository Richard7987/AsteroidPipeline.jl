# asteroidpipeline — the terminal app

A terminal interface for AsteroidPipeline.jl: give it an image set — a `.zip`,
a folder of FITS images, or the FITS files themselves — and it finds the
moving objects in it, tells known asteroids from new candidates, shows you the
images animated with every object marked, lets you decide whether each
candidate is real, and writes the MPC report.

Everything here is interface. The analysis is the AsteroidPipeline.jl library
(`../src`), run by a long-lived Julia process the app talks to
(`asteroid_ui/engine/server.jl`).

## Installing and running

See the repository's README: `pipx install asteroidpipeline` (PyPI),
`install.sh` on any Linux, or `nix run`. The package is built from `../pyproject.toml`; the wheel
carries the Julia library inside it (`asteroid_ui/julia/AsteroidPipeline`).

On first run (or `asteroidpipeline --setup`) the engine is set up
(`asteroid_ui/engine/__init__.py`): Julia is found — `$ASTEROIDPIPELINE_JULIA`,
`julia` on the PATH, or a copy downloaded by `juliapkg` — and a Julia
environment with the library is made in `~/.local/share/asteroid-ui/engine`
and precompiled. It is made again after an upgrade.

From a checkout, for development:

```bash
nix develop                 # or any Python 3.11+ with textual, textual-image, pillow, pyfiglet
python -m asteroid_ui       # from the repository root (the dev shell puts ui/ on PYTHONPATH)
```

Images are shown at full resolution in terminals with a graphics protocol
(Kitty, WezTerm, Ghostty, foot with Sixel…) and as coloured blocks elsewhere.

## Using it

1. **Settings** (`s`): who measures (the report's MEA line, required),
   optionally a contact (CON), and your designations for new candidates — a
   prefix and the next number (`NHU` and `2` give NHU0002, NHU0003…; empty,
   the program's own P2). Everything else in the report's header comes with
   the images — see *Telescope profiles* below.
2. **Choose an image set**: the start screen lists the `.zip` files of FITS
   images in your downloads folder; `o` opens any `.zip`, folder or FITS
   files. `t` analyzes every new `.zip` at once.
3. **Analysis**, 1–3 minutes per set: frames calibrated against Gaia DR3,
   damaged frames left out, moving objects searched for, known asteroids
   identified against the MPC's orbit file (downloaded daily) and SkyBoT.
4. **Check each candidate**: `Space` animates the images, `Enter` zooms on
   the object, then `y` (real) or `n` (not real, and why). Next to it are its
   astrometry (position and magnitude in every image, rate, direction, how
   straight its path is) and the known objects within 10′ of it, compared
   with it — *★ new?* if none moves like it, *≈ known?* if one moves exactly
   like it a little ahead or behind on its path (an object whose orbit is
   imprecise), *= known* if one is predicted right there. `k` marks it as
   that known object: kept as a real asteroid for the training data, left
   out of the report. For images up to 30 days old it is also compared with
   the MPC's NEO Confirmation Page: *◆ NEOCP* if an object awaiting
   confirmation moves like it — someone else's recent find.
5. **Report** (`r`): the observations of the candidates you confirmed, or the
   blank report if there are none — copy it (`c`) or save it (`s`). Each
   confirmed candidate gets your next designation, and keeps it whenever the
   report is written again; one on the NEOCP is reported under the NEOCP's
   designation.

`v` on the start screen checks any MPC report file before it is submitted:
missing `COD` line, wrong observatory code, lines left over from another image
set, known asteroids measured by mistake, single observations — and gives a
corrected version.

Every key is listed at the bottom of each screen and under `?`.

## Telescope profiles

FITS headers name the telescope (`TELESCOP`, `INSTRUME`, `ORIGIN`...) and
often give its site, but not its MPC code, its observers or its full
description — the report's COD, OBS and TEL lines. A *profile* holds those
for one telescope:

- the first time a telescope's images are opened, the app suggests the code —
  the MPC observatories within a few km of the site in the headers, ranked by
  how well their names match the telescope's (Pan-STARRS 1 → F51 rather than
  F52, ~10 m away; ZTF → I41) — and asks you to confirm it and fill OBS and
  TEL;
- from then on, every image set from that telescope gets its profile without
  asking;
- Pan-STARRS 1 (IASC's images) comes built in: F51, with the observers and
  telescope IASC's own settings use.

Profiles can be edited or removed in Settings. They are kept in
`~/.config/asteroid-ui/profiles.json`.

## Training data

Every candidate you confirm or reject is saved, with its pixels, to
`~/.local/share/asteroid-ui/dataset/`:

- `labels.csv` — one row per candidate: `label` (`real` / `bogus`), the
  `reason` it was rejected, `source` (`human`, or `catalog` for known
  asteroids, added automatically as confirmed real objects), speed,
  magnitude, observatory, epoch;
- `stamps/<set>/tracklet_<id>.fits` — a 63 × 63 × frames cube of the
  calibrated frames' pixels centred on the object in every frame.

That is the material for a classifier that can, one day, make the human
check unnecessary. Analyzing an image set again replaces its labels (its
candidates are numbered anew).

## Accessibility

- Everything works from the keyboard; the focused element has a visible
  border; every action is also a button.
- States are told by symbol and word, not by colour alone (✓ ✗ ⚠, □ ○ ◌).
- Text contrast is about 11:1 (Gruvbox dark).
- The animation is off until asked for, and never faster than 3 images per
  second.

## Files

| Path | What |
|:--|:--|
| `~/.config/asteroid-ui/config.json` | settings |
| `~/.config/asteroid-ui/profiles.json` | telescope profiles |
| `~/.config/asteroid-ui/designations.json` | the designations given to candidates |
| `~/.cache/asteroid-ui/sets/<set>/` | each analysis: `result.json`, display images, cut-outs |
| `~/.cache/asteroid-ui/engine.log` | the analysis engine's log |
| `~/.local/share/asteroid-ui/engine/` | the engine's Julia environment |
| `~/.cache/AsteroidPipeline/MPCORB.DAT` | the MPC's orbit file, refreshed daily |
| `~/.local/share/asteroid-ui/dataset/` | training data (above) |
