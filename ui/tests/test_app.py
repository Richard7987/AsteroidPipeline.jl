import asyncio
import json
import os
import shutil
import zipfile
from pathlib import Path

import pytest
from PIL import Image

from asteroid_ui import clipboard, config, designations, engine, labels, profiles, sets
from asteroid_ui.app import AsteroidApp
from asteroid_ui.backend import Backend, BackendError
from asteroid_ui.screens.analysis import AnalysisScreen
from asteroid_ui.screens.check import CheckScreen
from asteroid_ui.screens.confirm import ConfirmScreen
from asteroid_ui.screens.home import HomeScreen
from asteroid_ui.screens.report import ReportScreen, free_path
from asteroid_ui.screens.results import ResultsScreen, _dms, _hms
from asteroid_ui.screens.settings import SettingsScreen


@pytest.fixture(autouse=True)
def empty_folders():
    for folder in (config.CONFIG_FILE.parent, config.CACHE_DIR, config.DATA_DIR):
        shutil.rmtree(folder, ignore_errors=True)


def track(ra, epoch=2460000.5):
    return [{"frame": k, "ra": ra + k * 1e-3, "dec": 1.0, "epoch": epoch + k / 96, "mag": 19.0, "x": 40.0 + k, "y": 40.0}
            for k in (1, 2, 3)]


def tracklet(tid, ra=10.0, **kw):
    return {"id": tid, "label": f"P{tid}", "known": None, "speed_arcsec_h": 20.0, "pa_deg": 90.0, "rms_arcsec": 0.1,
            "mag": 19.0, "stamp": "", "check": None, "neocp": None, "points": track(ra), **kw}


def saved_result(name="a", tracklets=None, known=None, **kw):
    result = {"name": name, "observatory": "F51", "source": f"/x/{name}.zip", "image_size": [100, 100], "reliable": True,
              "report_kind": "candidates", "report": "", "warnings": [], "mpcorb": {},
              "frames": [{"index": k, "file": f"f{k}.fits", "epoch": 2460000.5 + k / 96, "time_utc": "2023-02-25 00:15:00",
                          "ratio": 1.0, "used": True, "image": f"frame_{k}.pgm"} for k in (1, 2, 3)],
              "tracklets": [tracklet(3)] if tracklets is None else tracklets, "known": known or [], "wcs_headers": [], **kw}
    folder = config.SETS_CACHE / name
    folder.mkdir(parents=True, exist_ok=True)
    for k in (1, 2, 3):
        Image.new("L", (100, 100), 30).save(folder / f"frame_{k}.pgm")
    (folder / "result.json").write_text(json.dumps(result))
    return result


def entry_of(result, source=Path("/x/a.zip")):
    return sets.SetEntry(result["name"], source, "2023-02-25", result)


class FakeBackend(Backend):
    """Answers as the engine would, without Julia."""

    def __init__(self, fail=None):
        super().__init__()
        self.ready.set()
        self.fail = fail
        self.calls = []

    async def start(self):
        pass

    async def stop(self):
        pass

    async def request(self, cmd, on_event=None, **args):
        self.calls.append((cmd, args))
        if cmd == self.fail:
            raise BackendError("the engine failed")
        if cmd == "inspect":
            return {"event": "inspection", "telescope": "PS1", "instrument": "gpc1", "candidates": []}
        if cmd == "build_report":
            return {"event": "report", "text": "COD F51\n"}
        if cmd == "check_report":
            return {"event": "report_check", "ok": True, "cleaned": args["text"], "issues": []}
        return {"event": "result"}


def run(scenario, backend=None, **settings):
    """Runs `scenario(app, pilot)` in the app; anything the app raises fails the test."""

    async def main():
        app = AsteroidApp()
        app.backend = backend or FakeBackend()
        app.config.measurers = "A. Name"
        for key, value in settings.items():
            setattr(app.config, key, value)
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await scenario(app, pilot)
            await pilot.pause()
        return app

    return asyncio.run(main())


def notifications(app):
    return [n.message for n in app._notifications]


class Settings:
    designation_prefix, next_designation = "NHU", 2

    def save(self):
        pass


def test_a_designation_follows_the_object_not_the_tracklet_number():
    c = Settings()
    assert designations.assign(c, "a", [tracklet(3, ra=10.0)]) == {"3": "NHU0002"}
    assert designations.assign(c, "a", [tracklet(3, ra=10.0)]) == {"3": "NHU0002"}
    # analyzed again: another object has tracklet 3, and the first one tracklet 5
    assert designations.assign(c, "a", [tracklet(3, ra=200.0)]) == {"3": "NHU0003"}
    assert designations.assign(c, "a", [tracklet(5, ra=10.0)]) == {"5": "NHU0002"}
    # another set's object in the same place isn't this one
    assert designations.assign(c, "b", [tracklet(5, ra=10.0)]) == {"5": "NHU0004"}


def test_designation_numbers_are_never_reused():
    c = Settings()
    designations.assign(c, "a", [tracklet(1, ra=10.0)])
    c.next_designation = 2
    assert designations.assign(c, "b", [tracklet(1, ra=50.0)]) == {"1": "NHU0003"}
    assert c.next_designation == 4


def test_designations_given_before_positions_were_kept_still_hold():
    designations._save({"assigned": {"a/3": "NHU0007"}})
    assert designations.assign(Settings(), "a", [tracklet(3)]) == {"3": "NHU0007"}


def test_too_long_a_number_uses_up_nothing():
    c = Settings()
    c.designation_prefix, c.next_designation = "ABCDE", 99
    assert designations.assign(c, "a", [tracklet(1)]) == {"1": "ABCDE99"}
    with pytest.raises(ValueError):
        designations.assign(c, "a", [tracklet(2, ra=50.0), tracklet(3, ra=60.0)])
    assert c.next_designation == 100


def test_decisions_survive_a_new_analysis_by_position():
    old = saved_result()
    labels.set_label("a", old, config.SETS_CACHE / "a", 3, labels.BOGUS, "noise")
    again = saved_result(tracklets=[tracklet(8, ra=10.0), tracklet(3, ra=200.0)])
    assert labels.add_catalog_labels("a", again, config.SETS_CACHE / "a") == 0
    found = labels.labels_for("a")
    assert list(found) == [8] and found[8]["label"] == labels.BOGUS and found[8]["reason"] == "noise"


def test_decisions_on_candidates_not_found_again_are_counted():
    old = saved_result()
    labels.set_label("a", old, config.SETS_CACHE / "a", 3, labels.REAL)
    assert labels.add_catalog_labels("a", saved_result(tracklets=[]), config.SETS_CACHE / "a") == 1


def test_analyzing_again_asks_when_there_are_decisions():
    result = saved_result()
    labels.set_label("a", result, config.SETS_CACHE / "a", 3, labels.REAL)

    async def scenario(app, pilot):
        home = app.screen
        home.entries = [entry_of(result)]
        from textual.widgets.option_list import Option

        home.query_one("#sets").clear_options()
        home.query_one("#sets").add_option(Option("a", id="0"))
        home.query_one("#sets").highlighted = 0
        await pilot.press("a")
        await pilot.pause()
        assert isinstance(app.screen, ConfirmScreen)
        await pilot.press("n")
        await pilot.pause()
        assert isinstance(app.screen, HomeScreen)

    run(scenario)


def test_images_without_a_telescope_name_ask_every_time():
    site = {"telescope": None, "instrument": None, "candidates": []}
    profile = profiles.Profile(profiles.profile_key(site), "679")
    profiles.save_profile(profile)
    assert profiles.find_profile(site) is None
    assert profiles.all_profiles() == [p for p in profiles.all_profiles() if p.key != profiles.UNNAMED]


def test_a_named_telescope_is_remembered():
    site = {"telescope": "T1", "candidates": []}
    profiles.save_profile(profiles.Profile("T1", "679"))
    assert profiles.find_profile(site).observatory == "679"


def test_the_report_has_the_profile_as_it_is_now():
    result = saved_result(profile={"key": "T1", "observatory": "F51", "observers": "", "telescope": "old", "builtin": False})
    profiles.save_profile(profiles.Profile("T1", "F51", "New Observer", "new TEL"))

    async def scenario(app, pilot):
        app.push_screen(ReportScreen(entry_of(result, None), {3: "real"}, {}))
        await pilot.pause(0.2)
        header = app.backend.calls[-1][1]["header"]
        assert header["observers"] == "New Observer" and header["telescope"] == "new TEL"

    run(scenario)


def test_a_failed_analysis_again_is_said_not_hidden():
    result = saved_result()

    async def scenario(app, pilot):
        app.push_screen(AnalysisScreen([entry_of(result)]))
        await pilot.pause(0.3)
        assert isinstance(app.screen, AnalysisScreen)
        assert "the engine failed" in str(app.screen.query_one("#analysis-error").render())

    run(scenario, FakeBackend(fail="analyze"))


def test_a_damaged_index_closes_nothing():
    config.DATASET_DIR.mkdir(parents=True)
    (config.DATASET_DIR / "labels.json").write_text("{broken")
    result = saved_result()

    async def scenario(app, pilot):
        await pilot.press("s")
        await pilot.pause(0.1)
        assert isinstance(app.screen, SettingsScreen)
        assert "damaged" in str(app.screen.query_one("#settings-note").render())
        await pilot.press("escape")
        await pilot.pause()
        app.push_screen(AnalysisScreen([sets.SetEntry("a", Path("/x/a.zip"), "2023", None)]))
        await pilot.pause(0.3)
        assert isinstance(app.screen, ResultsScreen)
        assert any("Labels not saved" in m for m in notifications(app))
        assert (config.DATASET_DIR / "labels.json").read_text() == "{broken"

    run(scenario)
    assert result


@pytest.mark.parametrize("encoding", ["utf-8", "latin-1"])
def test_the_check_reads_accents_in_either_encoding(tmp_path, encoding):
    path = tmp_path / "report.txt"
    path.write_text("COD F51\nOBS L. Olguín\n", encoding=encoding)

    async def scenario(app, pilot):
        app.push_screen(CheckScreen(None))
        await pilot.pause(0.3)
        assert app.screen.cleaned == "COD F51\nOBS L. Olguín\n"

    run(scenario, report_to_check=str(path))


def test_saving_a_report_never_replaces_a_different_one(tmp_path):
    path = tmp_path / "a_MPCReport.txt"
    assert free_path(path, "one") == path
    path.write_text("one")
    assert free_path(path, "one") == path
    assert free_path(path, "two") == tmp_path / "a_MPCReport-2.txt"


def test_missing_images_are_said(tmp_path):
    result = saved_result()
    (config.SETS_CACHE / "a" / "frame_2.pgm").unlink()

    async def scenario(app, pilot):
        await pilot.press("enter")
        await pilot.pause(0.2)
        assert isinstance(app.screen, HomeScreen)
        assert any("Analyze the set again" in m for m in notifications(app))

    run(scenario, downloads=str(tmp_path))
    assert result


def test_an_object_without_a_magnitude_is_shown():
    points = [{"frame": 1, "epoch": 2460000.5, "ra": 10, "dec": 1, "x": 5.0, "y": 5.0},
              {"frame": 2, "epoch": 2460000.51, "ra": 10, "dec": 1, "x": None, "y": None}]
    known = [{"name": "2020 AB", "packed": "", "mag": None, "source": "MPCORB", "detected": False, "tracklet": None, "points": points},
             {"name": "(1) Ceres", "packed": "", "mag": 9.0, "source": "MPCORB", "detected": False, "tracklet": None, "points": points}]
    result = saved_result(known=known, tracklets=[tracklet(3, mag=None, rms_arcsec=None)])

    async def scenario(app, pilot):
        app.push_screen(ResultsScreen(entry_of(result)))
        await pilot.pause(0.2)
        screen = app.screen
        assert [o.label for o in screen.images.objects][-2:] == ["(1) Ceres", "2020 AB"]
        screen._show_details("k:2020 AB")
        screen._show_details("t:3")

    run(scenario)


def test_faint_candidates_are_marked_and_listed_last_but_stay():
    result = saved_result(tracklets=[tracklet(1, ra=10.0, below_limit=True), tracklet(2, ra=50.0)], limiting_mag=21.5)

    async def scenario(app, pilot):
        app.push_screen(ResultsScreen(entry_of(result)))
        await pilot.pause(0.2)
        screen = app.screen
        assert [o.key for o in screen.images.objects] == ["t:2", "t:1"]
        assert screen._candidate_ids() == [1, 2]
        listed = [str(screen.query_one("#objects").get_option_at_index(i).prompt)
                  for i in range(screen.query_one("#objects").option_count)]
        assert any("P1" in line and "faint" in line for line in listed)
        assert not any("P2" in line and "faint" in line for line in listed)
        screen._show_details("t:1")
        assert "limit V 21.5" in str(screen.query_one("#details").render())

    run(scenario)


def test_results_saved_before_the_limit_still_open():
    result = saved_result()

    async def scenario(app, pilot):
        app.push_screen(ResultsScreen(entry_of(result)))
        await pilot.pause(0.2)
        app.screen._show_details("t:3")
        assert "faint" not in str(app.screen.query_one("#details").render())

    run(scenario)


@pytest.mark.parametrize("text", ["From: /home/me/[/x/a.zip", "[/red] done", "Installing: ⌅ [a93c6f00] DataFrames"])
def test_outside_text_is_not_read_as_markup(text):
    async def scenario(app, pilot):
        app.notify(text)
        app.engine_status = text
        await pilot.pause(0.1)
        assert text in str(app.screen.query_one("#status").render())

    run(scenario)


def test_a_folder_with_brackets_in_its_path_opens(tmp_path):
    folder = tmp_path / "[" / "x]"
    folder.mkdir(parents=True)
    weird = tmp_path / "[" / "[" / "a.zip"

    async def scenario(app, pilot):
        await pilot.pause()
        app.push_screen(AnalysisScreen([sets.SetEntry("a", weird, "2023", None)]))
        await pilot.pause(0.3)

    run(scenario, FakeBackend(fail="inspect"), downloads=str(folder))


def test_a_number_that_only_looks_like_one_is_refused(tmp_path):
    async def scenario(app, pilot):
        await pilot.press("s")
        await pilot.pause(0.1)
        app.screen.query_one("#set-next_designation").value = "²"
        await pilot.press("ctrl+s")
        await pilot.pause(0.1)
        assert isinstance(app.screen, SettingsScreen)

    run(scenario, downloads=str(tmp_path))


def test_a_config_edited_by_hand_starts_the_app():
    config.CONFIG_FILE.parent.mkdir(parents=True)
    config.CONFIG_FILE.write_text(json.dumps({"fps": "fast", "measurers": "A. Name", "next_designation": "two"}))
    loaded = config.load_config()
    assert loaded.fps == 1.5 and loaded.measurers == "A. Name" and loaded.next_designation == 1
    config.CONFIG_FILE.write_text("[1, 2]")
    assert config.load_config().measurers == ""


def test_seconds_never_round_to_sixty():
    assert _hms(15 * (1 + 59 / 60 + 59.9999 / 3600)) == "02 00 00.000"
    assert _dms(-(1 + 59 / 60 + 59.999 / 3600)) == "-02 00 00.00"
    assert _hms(15 * (23 + 59 / 60 + 59.9999 / 3600)) == "00 00 00.000"
    assert _hms(36.2344) == "02 24 56.256" and _dms(2.06997) == "+02 04 11.89"


def _zip(path, names):
    with zipfile.ZipFile(path, "w") as z:
        for name in names:
            z.writestr(name, "x")


def test_discovery_finds_any_case_and_ignores_macos_files(tmp_path):
    _zip(tmp_path / "SET.ZIP", ["f1.fits", "f2.fits", "f3.fits"])
    _zip(tmp_path / "too_few.zip", ["f1.fits", "__MACOSX/._f1.fits", "__MACOSX/._f2.fits", "__MACOSX/._f3.fits"])
    deep = tmp_path / "1" / "2" / "3" / "4"
    deep.mkdir(parents=True)
    _zip(deep / "too_deep.zip", ["f1.fits", "f2.fits", "f3.fits"])
    assert [e.name for e in sets.discover(str(tmp_path))] == ["SET"]


def test_same_track_across_ra_zero():
    a = [{"epoch": 2460000.5, "ra": 359.9999, "dec": 0.0}]
    b = [{"epoch": 2460000.5, "ra": 0.0001, "dec": 0.0}]
    assert sets.same_track(a, b)
    assert not sets.same_track(a, [{"epoch": 2460001.5, "ra": 359.9999, "dec": 0.0}])


def test_the_clipboard_tries_the_next_tool(tmp_path, monkeypatch):
    (tmp_path / "wl-copy").write_text("#!/bin/sh\nexit 1\n")
    (tmp_path / "xclip").write_text(f"#!/bin/sh\ncat > {tmp_path}/got\n")
    for tool in ("wl-copy", "xclip"):
        (tmp_path / tool).chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}:{os.environ['PATH']}")

    class App:
        def copy_to_clipboard(self, text):
            pass

    clipboard.copy_text(App(), "hello")
    assert (tmp_path / "got").read_text() == "hello"


def test_a_julia_that_does_not_run_is_said(monkeypatch):
    monkeypatch.setenv("ASTEROIDPIPELINE_JULIA", "/nonexistent/julia")

    async def main():
        backend = Backend()
        await backend.start()
        await asyncio.wait_for(backend.ready.wait(), 5)
        return backend.failed

    assert "/nonexistent/julia" in asyncio.run(main())
    with pytest.raises(engine.EngineError):
        engine.julia_version("/nonexistent/julia")


def test_setup_script_quotes_any_path():
    assert engine._julia_string("/home/o'brien/$x\\y") == '"/home/o\'brien/\\$x\\\\y"'
