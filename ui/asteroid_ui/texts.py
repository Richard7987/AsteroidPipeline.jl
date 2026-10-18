"""Every piece of text the user reads, in one place.

The engine (Julia) answers with codes (`refine`, `missing_cod`...); they are
worded here. Every state carries a symbol as well as a colour, so it reads the
same to someone who doesn't tell the colours apart.
"""

# Analysis steps, in order, keyed by the engine's stage codes.
STAGES: list[tuple[str, str]] = [
    ("extract", "Unpack the images"),
    ("clean", "Clean invalid pixels"),
    ("refine", "Calibrate positions against Gaia stars"),
    ("search", "Search for moving objects"),
    ("identify", "Identify known asteroids"),
    ("nearby", "Check the known objects near each candidate"),
    ("field", "Predict the known asteroids in the field"),
    ("report", "Write the MPC report"),
    ("images", "Prepare the images"),
]

STATUS_MESSAGES = {
    "mpcorb_download": "Updating the MPC's orbit file (MPCORB.DAT)…",
    "mpcorb_load": "Loading the MPC's orbits…",
}

ICON_DONE, ICON_RUNNING, ICON_PENDING = "✓", "▶", "·"
ICON_ERROR, ICON_WARNING, ICON_OK = "✗", "⚠", "✓"
ICON_CANDIDATE, ICON_KNOWN, ICON_KNOWN_MISSED = "□", "○", "◌"

# The nearby-object check of a candidate -> (tag, style, sentence)
CHECKS = {
    "new": ("★ new?", "bold magenta",
            "No known object within 10′ moves like it: it may be new."),
    "possibly_known": ("≈ known?", "bold cyan",
                       "A known object nearby moves exactly like it — probably that object, with an imprecise "
                       "orbit. If so, press k."),
    "known": ("= known", "bold green", "A known object is predicted right there: it is that object. Press k."),
}

# Verdict of an analysis -> (icon, title, style, explanation)
VERDICTS = {
    "blank": (ICON_OK, "No new objects", "success",
              "Nothing new was found. The report is the blank one: copy it and submit it."),
    "candidates": (ICON_WARNING, "Candidates to check", "warning",
                   "Look at each candidate in the animation (Space) and decide: y = real, n = not real. "
                   "The report includes the real ones only."),
    "checked": (ICON_OK, "Candidates checked", "success",
                "Every candidate has been decided. The report includes the ones marked real."),
    "unreliable": (ICON_ERROR, "Damaged image set", "error",
                   "The images have too many defects to search them. Look at them yourself, "
                   "and tell whoever coordinates your observations."),
}

REPORT_TITLES = {
    "blank": "Blank report: no new objects",
    "confirmed": "Report of the candidates you confirmed as real",
    "unreliable": "Blank report for a damaged image set",
}

REPORT_HELP = {
    "blank": "Submit this text for this image set. It starts with the COD line every MPC report needs.",
    "confirmed": "Submit this text for this image set: it has the observations of every candidate "
                 "you marked as real, with temporary designations (P…).",
    "unreliable": "Nothing could be searched in these images. If they can't be measured by eye either, "
                  "submit this blank report and say the set was damaged.",
}

# Problems the report check finds: code -> (is_error, text); {detail} is the
# specific value (a code, a designation...).
ISSUES = {
    "empty": (True, "The report is empty."),
    "missing_cod": (True, "The «COD» line is missing: every report must start with it. Use the corrected report."),
    "cod_not_first": (True, "The «COD» line isn't the first line; it has to be."),
    "wrong_code": (True, "The observatory code is «{detail}», not the one in Settings for these images."),
    "station_mismatch": (True, "An observation has the code «{detail}» instead of this observatory's."),
    "malformed_line": (True, "A line isn't in the MPC format: «{detail}»."),
    "other_set": (True, "The observation of «{detail}» belongs to another image set (another time or "
                        "another part of the sky): it was left over in the report file."),
    "known_object": (False, "«{detail}» is a known asteroid — probably measured by clicking on it. "
                            "Remove it unless known objects are meant to be reported."),
    "single_observation": (False, "«{detail}» has a single observation; the MPC needs at least two."),
    "missing_end": (False, "The final «----- end -----» line is missing."),
}

HELP_HOME = """\
[b]How it works[/b]

[b]1.[/b] Choose an image set: a [b].zip[/b] from the list, or press [b]o[/b] to open a .zip, a folder of FITS images, or the FITS files themselves. [b]t[/b] analyzes every new .zip at once.
[b]2.[/b] The images are searched for moving objects; known asteroids are identified.
[b]3.[/b] Check each new candidate in the animated images: [b]y[/b] real, [b]n[/b] not real.
[b]4.[/b] Copy the MPC report ([b]r[/b]) and submit it.

Every decision is saved, with the candidate's pixels, to train a classifier that can one day do step 3 on its own.
"""

HELP_KEYS = """\
[b]Keys in an image set[/b]

[b]← →[/b]      previous / next image
[b]Space[/b]    start or stop the animation
[b]+ −[/b]      animation faster / slower (up to 3 images/s)
[b]↑ ↓[/b]      choose an object in the list (Enter: zoom on it)
[b]y[/b]        the candidate is a real moving object
[b]n[/b]        the candidate is not real (then say why)
[b]k[/b]        the candidate is the known object found near it
[b]u[/b]        undo the decision on the candidate
[b]m[/b]        show or hide the marks
[b]z[/b]        zoom on the chosen object (again: whole image)
[b]r[/b]        MPC report
[b]Esc[/b]      back

[b]Marks[/b]
□ yellow square   new candidate (✓ real, ✗ not real once decided)
○ green circle    known asteroid, detected
◌ dashed circle   known asteroid, not detected (too faint, or hidden)
"""
