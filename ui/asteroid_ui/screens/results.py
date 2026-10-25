"""One analyzed image set: verdict, objects, the animated images, the human check."""

from __future__ import annotations

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import Screen
from textual.widgets import Button, Footer, Label, OptionList, Static
from textual.widgets.option_list import Option

from asteroid_ui.imaging import CANDIDATE, KNOWN, MISSED, SetImages
from asteroid_ui.labels import BOGUS, REAL, REASONS, labels_for, set_label
from asteroid_ui.profiles import profile_of
from asteroid_ui.sets import SetEntry
from asteroid_ui.texts import CHECKS, ICON_CANDIDATE, ICON_KNOWN, ICON_KNOWN_MISSED, VERDICTS
from asteroid_ui.widgets import SkyViewer, StatusBar

_KIND_ICON = {CANDIDATE: (ICON_CANDIDATE, "bold yellow"), KNOWN: (ICON_KNOWN, "bold green"),
              MISSED: (ICON_KNOWN_MISSED, "cyan")}
_DECISION_TEXT = {REAL: ("✓ real", "bold green"), BOGUS: ("✗ not real", "bold red"), None: ("? to check", "yellow")}
_MAX_FPS = 3.0  # no faster: flashing above 3 per second is an accessibility hazard

LEGEND = Text.assemble(("□ new candidate", "bold yellow"), "    ", ("○ known, detected", "bold green"), "    ",
                       ("◌ known, not detected", "cyan"))


def _hms(ra: float) -> str:
    h = ra / 15
    return f"{int(h):02d} {int(h * 60 % 60):02d} {h * 3600 % 60:06.3f}"


def _dms(dec: float) -> str:
    d = abs(dec)
    return f"{'-' if dec < 0 else '+'}{int(d):02d} {int(d * 60 % 60):02d} {d * 3600 % 60:05.2f}"


def _distance(arcsec: float) -> str:
    return f"{arcsec:.0f}\"" if arcsec < 120 else f"{arcsec / 60:.1f}′"


def verdict_kind(result: dict, decisions: dict[int, str]) -> str:
    if result.get("report_kind") == "unreliable":
        return "unreliable"
    candidates = [t["id"] for t in result["tracklets"] if not t.get("known")]
    if not candidates:
        return "blank"
    return "checked" if all(c in decisions for c in candidates) else "candidates"


class ResultsScreen(Screen):
    BINDINGS = [
        Binding("left", "prev", "Previous", priority=True),
        Binding("right", "next", "Next", priority=True),
        Binding("space", "toggle_animation", "Animate", priority=True),
        Binding("y", "real", "Real"),
        Binding("n", "bogus", "Not real"),
        Binding("k", "known", "Known object"),
        Binding("u", "undo", "Undo", show=False),
        Binding("m", "toggle_marks", "Marks"),
        Binding("z", "toggle_zoom", "Zoom"),
        Binding("plus,equals_sign", "faster", "Faster", show=False),
        Binding("minus", "slower", "Slower", show=False),
        Binding("r", "report", "MPC report"),
        Binding("question_mark", "help", "Help"),
        Binding("escape", "back", "Back"),
    ]

    def __init__(self, entry: SetEntry) -> None:
        super().__init__()
        self.entry = entry
        self.result = entry.result or {}
        self.images = SetImages(entry.cache_dir, self.result)
        self.decisions: dict[int, str] = {}
        self.reasons: dict[int, str] = {}
        self.fps = 1.5  # the saved speed is read once mounted
        self.animating = False

    def compose(self) -> ComposeResult:
        yield Static("", id="verdict")
        with Horizontal(id="results-body"):
            with Vertical(id="sidebar"):
                yield Label("Objects", classes="section-title")
                yield OptionList(id="objects")
                yield Static("", id="details")
                yield Label("Images", classes="section-title")
                yield Static("", id="frames-info")
                yield Static("", id="warnings")
            with Vertical(id="viewer-pane"):
                yield SkyViewer(self.images, id="viewer")
                yield Static("", id="caption")
                yield Static(LEGEND, id="legend")
                with Horizontal(id="controls", classes="buttons"):
                    yield Button("◀ Prev", id="prev")
                    yield Button("▶ Animate", id="play", variant="primary")
                    yield Button("Next ▶", id="next")
                    yield Button("Marks: on", id="marks")
                    yield Button("Zoom", id="zoom")
                    yield Button("✓ Real (y)", id="real", variant="success")
                    yield Button("✗ Not real (n)", id="bogus", variant="error")
                    yield Button("Report (r)", id="report")
        yield StatusBar(id="status")
        yield Footer()

    def on_mount(self) -> None:
        self.fps = min(_MAX_FPS, max(0.5, float(self.app.config.fps)))
        self.timer = self.set_interval(1 / self.fps, self._advance, pause=True)
        human = {tid: row for tid, row in labels_for(self.entry.name).items() if row["source"] == "human"}
        self.decisions = {tid: row["label"] for tid, row in human.items()}
        self.reasons = {tid: row.get("reason", "") for tid, row in human.items()}
        self.images.set_decisions(self.decisions)
        self._fill_objects()
        self._fill_frames()
        self._update_verdict()
        self._update_caption()
        self.query_one("#objects", OptionList).focus()

    # Side panel
    def _candidate_ids(self) -> list[int]:
        return [t["id"] for t in self.result.get("tracklets", []) if not t.get("known")]

    def _update_verdict(self) -> None:
        kind = verdict_kind(self.result, self.decisions)
        icon, title, style, explanation = VERDICTS[kind]
        ids = self._candidate_ids()
        counts = ""
        if ids:
            real = sum(1 for i in ids if self.decisions.get(i) == REAL)
            bogus = sum(1 for i in ids if self.decisions.get(i) == BOGUS)
            counts = f"   [{real} real · {bogus} not real · {len(ids) - real - bogus} to check]"
        verdict = self.query_one("#verdict", Static)
        verdict.update(Text.assemble((f" {icon} {self.entry.name} — {title}{counts} ", "bold"), "\n ", explanation))
        verdict.set_classes(f"verdict-{style}")

    def _fill_objects(self) -> None:
        options = self.query_one("#objects", OptionList)
        keep = options.highlighted
        options.clear_options()
        groups = [(CANDIDATE, "New candidates"), (KNOWN, "Known asteroids, detected"),
                  (MISSED, "Known asteroids, not detected")]
        first = None
        for kind, title in groups:
            objs = [o for o in self.images.objects if o.kind == kind]
            options.add_option(Option(Text(f"{title} ({len(objs)})", style="bold underline"), disabled=True))
            for o in objs:
                icon, style = _KIND_ICON[kind]
                parts = [(f" {icon} ", style), (o.data["label"] if kind == CANDIDATE else o.label, "bold"),
                         f"  {o.detail}"]
                if kind == CANDIDATE:
                    tid = o.data["id"]
                    if self.reasons.get(tid, "").startswith("known:"):
                        text, dstyle = "= known", "bold green"
                    elif tid in self.decisions:
                        text, dstyle = _DECISION_TEXT[self.decisions[tid]]
                    else:
                        check = o.data.get("check") or {}
                        neocp = (o.data.get("neocp") or {}).get("verdict") == "on_neocp"
                        verdict = "neocp" if neocp else check.get("verdict")
                        text, dstyle = CHECKS[verdict][:2] if verdict in CHECKS else _DECISION_TEXT[None]
                    parts += ["  ", (text, dstyle)]
                options.add_option(Option(Text.assemble(*parts), id=o.key))
                first = first or o.key
            if not objs:
                options.add_option(Option(Text("   none", style="dim italic"), disabled=True))
        if keep is not None:
            options.highlighted = keep
        elif first is not None:
            options.highlighted = self._index_of(first)

    def _index_of(self, key: str) -> int | None:
        options = self.query_one("#objects", OptionList)
        return next((i for i in range(options.option_count) if options.get_option_at_index(i).id == key), None)

    def _fill_frames(self) -> None:
        lines = Text()
        for f in self.result.get("frames", []):
            ok = f["used"]
            lines.append(f" {f['index']}  {f['time_utc'][11:]} UTC  ")
            lines.append("✓ good\n" if ok else "✗ defects — left out\n", style="green" if ok else "bold red")
        self.query_one("#frames-info", Static).update(lines)
        notes = Text()
        profile = profile_of(self.result)
        notes.append(f"Observatory: {profile.describe() if profile else self.result.get('observatory', '?')}\n",
                     style="dim")
        orbit = self.result.get("mpcorb") or {}
        if orbit.get("path"):
            notes.append(f"Orbits: MPC file of {orbit.get('updated', '')[:10]}\n", style="dim")
        for w in self.result.get("warnings") or []:
            notes.append(f"⚠ {w}\n", style="yellow")
        self.query_one("#warnings", Static).update(notes)

    def on_option_list_option_highlighted(self, event: OptionList.OptionHighlighted) -> None:
        if event.option_list.id != "objects" or event.option.id is None:
            return
        viewer = self.query_one(SkyViewer)
        viewer.selected = event.option.id
        if viewer.zoom_on is not None:
            viewer.zoom_on = event.option.id
        self._show_details(event.option.id)
        self._update_caption()

    def _show_details(self, key: str | None) -> None:
        """The chosen object's astrometry and, for a candidate, the known objects near it."""
        obj = self.images.find(key)
        details = self.query_one("#details", Static)
        if obj is None:
            details.update("")
            return
        text = Text()
        if obj.kind == CANDIDATE:
            t = obj.data
            text.append(f"{t['label']}", style="bold")
            text.append(f"  {t['speed_arcsec_h']:.1f}\"/h toward {t.get('pa_deg', 0):.0f}° (PA)"
                        f"  · straight line ±{t.get('rms_arcsec', 0):.2f}\"\n")
            frames = {f["index"]: f for f in self.result.get("frames", [])}
            for p in t["points"]:
                mag = f"{p['mag']:.1f}" if p.get("mag") is not None else " ? "
                text.append(f" {p['frame']}  {frames.get(p['frame'], {}).get('time_utc', '')[11:]}  "
                            f"{_hms(p['ra'])} {_dms(p['dec'])}  {mag}\n", style="dim")
            check = t.get("check")
            if check:
                tag, style, sentence = CHECKS[check["verdict"]]
                text.append(f"{tag} ", style=style)
                text.append(sentence + "\n")
                near = [o for o in check["objects"] if o["agreement"] != "none"] or check["objects"][:2]
                for o in near[:3]:
                    how = {"position": "right there", "motion": "same motion"}.get(o["agreement"], "different motion")
                    text.append(f"  {o['name']}  {_distance(o['separation_arcsec'])} away, V {o['mag']:.1f}, "
                                f"{o['rate_arcsec_h']:.1f}\"/h to {o['pa_deg']:.0f}° — {how}\n",
                                style="" if o["agreement"] != "none" else "dim")
            neocp = t.get("neocp")
            if neocp and neocp.get("objects"):
                if neocp["verdict"] == "on_neocp":
                    tag, style, sentence = CHECKS["neocp"]
                    text.append(f"{tag} ", style=style)
                    text.append(sentence + "\n")
                for o in neocp["objects"][:2]:
                    how = "same motion" if o["agreement"] == "motion" else "different motion"
                    flag = " (uncertain ephemeris)" if o["uncertain"] else ""
                    text.append(f"  NEOCP {o['designation']}  {_distance(o['separation_arcsec'])} away, "
                                f"{o['rate_arcsec_h']:.0f}\"/h to {o['pa_deg']:.0f}° — {how}{flag}\n",
                                style="" if o["agreement"] == "motion" else "dim")
        else:
            k = obj.data
            text.append(f"{k['name']}", style="bold")
            text.append(f"  V {k['mag']:.1f} predicted · {'detected' if k['detected'] else 'not detected'}"
                        f" · position from the {k.get('source', 'catalog')}\n")
        details.update(text)

    def on_option_list_option_selected(self, event: OptionList.OptionSelected) -> None:
        if event.option_list.id == "objects" and event.option.id is not None:
            self.action_toggle_zoom()

    # The human check
    def _selected_candidate(self) -> int | None:
        key = self.query_one(SkyViewer).selected
        obj = self.images.find(key)
        if obj is None or obj.kind != CANDIDATE:
            self.app.notify("Choose a new candidate in the list first (↑ ↓).")
            return None
        return obj.data["id"]

    def _decide(self, tracklet_id: int, label: str | None, reason: str = "") -> None:
        set_label(self.entry.name, self.result, self.entry.cache_dir, tracklet_id, label, reason)
        if label is None:
            self.decisions.pop(tracklet_id, None)
            self.reasons.pop(tracklet_id, None)
        else:
            self.decisions[tracklet_id] = label
            self.reasons[tracklet_id] = reason
        self.images.set_decisions(self.decisions)
        self._fill_objects()
        self._update_verdict()
        self.query_one(SkyViewer).redraw()
        if label is not None:
            undecided = [i for i in self._candidate_ids() if i not in self.decisions]
            if undecided:  # straight on to the next one
                index = self._index_of(f"t:{undecided[0]}")
                if index is not None:
                    self.query_one("#objects", OptionList).highlighted = index
            else:
                self.app.notify("Every candidate is decided. Press r for the MPC report.")

    def action_real(self) -> None:
        tid = self._selected_candidate()
        if tid is not None:
            self._decide(tid, REAL)

    def action_bogus(self) -> None:
        tid = self._selected_candidate()
        if tid is None:
            return
        from asteroid_ui.screens.reason import ReasonScreen

        label = next(t["label"] for t in self.result["tracklets"] if t["id"] == tid)

        def chosen(reason: str | None) -> None:
            if reason is not None:
                self._decide(tid, BOGUS, reason)
                self.app.notify(f"{label}: not real — {REASONS[reason].split(' (')[0].lower()}.")

        self.app.push_screen(ReasonScreen(label), chosen)

    def action_known(self) -> None:
        tid = self._selected_candidate()
        if tid is None:
            return
        t = next(t for t in self.result["tracklets"] if t["id"] == tid)
        near = [o for o in (t.get("check") or {}).get("objects", []) if o["agreement"] in ("position", "motion")]
        if not near:
            self.app.notify(f"No known object near {t['label']} moves like it.", severity="warning")
            return
        name = near[0]["name"]
        self._decide(tid, REAL, f"known:{name}")
        self.app.notify(f"{t['label']} is {name}: kept as a real asteroid, left out of the report.")

    def action_undo(self) -> None:
        tid = self._selected_candidate()
        if tid is not None and tid in self.decisions:
            self._decide(tid, None)

    # Viewer
    def _update_caption(self) -> None:
        viewer = self.query_one(SkyViewer)
        frames = self.result.get("frames", [])
        f = frames[viewer.frame - 1] if frames else {}
        parts = [f"Image {viewer.frame} of {self.images.n_frames}", f"{f.get('time_utc', '')} UTC"]
        if f and not f.get("used", True):
            parts.append("left out (defects)")
        parts.append(f"▶ animating, {self.fps:g}/s" if self.animating else "animation off")
        parts.append("marks on" if viewer.marks else "marks off")
        target = self.images.find(viewer.zoom_on)
        parts.append(f"zoom on {target.label}" if target else "whole image")
        self.query_one("#caption", Static).update("  ·  ".join(parts))
        self.query_one("#play", Button).label = "■ Stop" if self.animating else "▶ Animate"
        self.query_one("#marks", Button).label = "Marks: on" if viewer.marks else "Marks: off"
        self.query_one("#zoom", Button).label = "Whole image" if target else "Zoom"

    def _step(self, delta: int) -> None:
        viewer = self.query_one(SkyViewer)
        viewer.frame = (viewer.frame - 1 + delta) % self.images.n_frames + 1
        self._update_caption()

    def _advance(self) -> None:
        self._step(1)

    def action_prev(self) -> None:
        self._step(-1)

    def action_next(self) -> None:
        self._step(1)

    def action_toggle_animation(self) -> None:
        self.animating = not self.animating
        if self.animating:
            self.timer.resume()
        else:
            self.timer.pause()
        self._update_caption()

    def _set_speed(self, fps: float) -> None:
        self.fps = min(_MAX_FPS, max(0.5, fps))
        self.timer.stop()
        self.timer = self.set_interval(1 / self.fps, self._advance, pause=not self.animating)
        self.app.config.fps = self.fps
        self.app.config.save()
        self._update_caption()

    def action_faster(self) -> None:
        self._set_speed(self.fps + 0.5)

    def action_slower(self) -> None:
        self._set_speed(self.fps - 0.5)

    def action_toggle_marks(self) -> None:
        viewer = self.query_one(SkyViewer)
        viewer.marks = not viewer.marks
        self._update_caption()

    def action_toggle_zoom(self) -> None:
        viewer = self.query_one(SkyViewer)
        if viewer.zoom_on is not None:
            viewer.zoom_on = None
        elif viewer.selected is not None:
            viewer.zoom_on = viewer.selected
        else:
            self.app.notify("Choose an object in the list first.")
        self._update_caption()

    # Other screens
    def action_report(self) -> None:
        from asteroid_ui.screens.report import ReportScreen

        self.app.push_screen(ReportScreen(self.entry, self.decisions, self.reasons))

    def action_help(self) -> None:
        from asteroid_ui.screens.help import HelpScreen

        self.app.push_screen(HelpScreen())

    def action_back(self) -> None:
        self.timer.stop()
        self.app.pop_screen()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        actions = {"prev": self.action_prev, "next": self.action_next, "play": self.action_toggle_animation,
                   "marks": self.action_toggle_marks, "zoom": self.action_toggle_zoom, "real": self.action_real,
                   "bogus": self.action_bogus, "report": self.action_report}
        if event.button.id in actions:
            actions[event.button.id]()
