"""The set's images, aligned, with the objects marked on them.

The engine leaves every image already aligned to the first one (the telescope's
small shifts between exposures) and in 8 bits (frame_k.pgm), and every object's
positions in FITS pixels of the first image. Here they are only cropped, scaled
to the panel and marked, so going from image to image (the animation) is fast.

Marks differ in shape, not only colour: square = new candidate, circle = known
asteroid detected, dashed circle = known asteroid not detected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont

CANDIDATE, KNOWN, MISSED = "candidate", "known", "missed"

COLORS = {
    CANDIDATE: (250, 189, 47),    # yellow
    KNOWN: (184, 187, 38),        # green
    MISSED: (131, 165, 152),      # grey-blue
    "selected": (251, 241, 199),  # near white
}
BACKGROUND = (29, 32, 33)
ZOOM_HALF_SIZE = 150  # image pixels on each side of the object when zoomed


@dataclass
class SkyObject:
    key: str
    kind: str
    label: str
    detail: str
    points: dict[int, tuple[float, float]]  # image (1..n) -> FITS (x, y)
    data: dict[str, Any] = field(default_factory=dict)

    def center(self) -> tuple[float, float]:
        xs, ys = zip(*self.points.values())
        return sum(xs) / len(xs), sum(ys) / len(ys)


DECISION_MARK = {"real": " ✓", "bogus": " ✗"}


def objects_from_result(result: dict[str, Any], decisions: dict[int, str] | None = None) -> list[SkyObject]:
    """New candidates first, then known asteroids detected, then those not detected.

    `decisions` (tracklet id -> "real" / "bogus") adds ✓ / ✗ to candidates' labels.
    A detected known asteroid is drawn where it was measured, not where its
    orbit predicts it (a two-body orbit can be off by an arcsecond or more).
    """
    decisions = decisions or {}
    tracklets = {t["id"]: t for t in result["tracklets"]}
    objects: list[SkyObject] = []
    for t in result["tracklets"]:
        if t.get("known"):
            continue
        mag = f"V {t['mag']:.1f}" if t.get("mag") else "V ?"
        objects.append(SkyObject(
            f"t:{t['id']}", CANDIDATE, t["label"] + DECISION_MARK.get(decisions.get(t["id"], ""), ""),
            f"{mag} · {t['speed_arcsec_h']:.0f}\"/h", {p["frame"]: (p["x"], p["y"]) for p in t["points"]}, t))
    known = sorted(result["known"], key=lambda k: (not k["detected"], k["mag"]))
    for k in known:
        kind = KNOWN if k["detected"] else MISSED
        points = {p["frame"]: (p["x"], p["y"]) for p in k["points"]}
        measured = tracklets.get(k.get("tracklet"))
        if measured is not None:
            points.update({p["frame"]: (p["x"], p["y"]) for p in measured["points"]})
        objects.append(SkyObject(f"k:{k['name']}", kind, k["name"], f"V {k['mag']:.1f}", points, k))
    return objects


@lru_cache(maxsize=8)
def _font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    try:
        return ImageFont.load_default(size=size)
    except TypeError:  # Pillow < 10.1
        return ImageFont.load_default()


class SetImages:
    def __init__(self, result_dir: Path, result: dict[str, Any]) -> None:
        self.nx, self.ny = result["image_size"]
        self.frames = [Image.open(result_dir / f["image"]).convert("L") for f in result["frames"]]
        self.result = result
        self.objects = objects_from_result(result)
        self._cache: dict[tuple, Image.Image] = {}

    def set_decisions(self, decisions: dict[int, str]) -> None:
        self.objects = objects_from_result(self.result, decisions)
        self._cache.clear()

    @property
    def n_frames(self) -> int:
        return len(self.frames)

    def find(self, key: str | None) -> SkyObject | None:
        return next((o for o in self.objects if o.key == key), None)

    def _to_image(self, x: float, y: float) -> tuple[float, float]:
        """FITS pixel (x, y; y up, from 1) -> pixel of the stored image."""
        return x - 1, self.ny - y

    def render(self, frame: int, size: tuple[int, int], *, marks: bool = True,
               zoom_on: str | None = None, selected: str | None = None) -> Image.Image:
        """Image `frame` (from 1) on a canvas of `size` pixels, undistorted."""
        key = (frame, size, marks, zoom_on, selected)
        if key in self._cache:
            return self._cache[key]
        if len(self._cache) > 64:
            self._cache.clear()
        width, height = max(size[0], 16), max(size[1], 16)

        box = (0, 0, self.nx, self.ny)
        target = self.find(zoom_on)
        if target is not None:
            cx, cy = self._to_image(*target.center())
            half = ZOOM_HALF_SIZE * max(1.0, _span(target) / (1.4 * ZOOM_HALF_SIZE))
            # the panel's proportions, so no empty bands are left
            hx, hy = (half * width / height, half) if width >= height else (half, half * height / width)
            box = (int(cx - hx), int(cy - hy), int(cx + hx), int(cy + hy))

        crop = self.frames[frame - 1].crop(box)
        scale = min(width / crop.width, height / crop.height)
        shown = crop.resize((max(1, round(crop.width * scale)), max(1, round(crop.height * scale))),
                            Image.Resampling.LANCZOS if scale < 1 else Image.Resampling.NEAREST)
        canvas = Image.new("RGB", (width, height), BACKGROUND)
        ox, oy = (width - shown.width) // 2, (height - shown.height) // 2
        canvas.paste(shown.convert("RGB"), (ox, oy))

        if marks:
            def to_canvas(x: float, y: float) -> tuple[float, float]:
                ix, iy = self._to_image(x, y)
                return ox + (ix - box[0] + 0.5) * scale, oy + (iy - box[1] + 0.5) * scale

            self._draw_marks(canvas, frame, to_canvas, selected, zoomed=target is not None)
        self._cache[key] = canvas
        return canvas

    def _draw_marks(self, canvas: Image.Image, frame: int, to_canvas, selected: str | None, zoomed: bool) -> None:
        draw = ImageDraw.Draw(canvas, "RGBA")
        short = min(canvas.size)
        radius = max(7, round(short / (14 if zoomed else 45)))
        width = max(2, radius // 5)
        font = _font(max(12, round(short / 45)))
        for obj in sorted(self.objects, key=lambda o: o.key == selected):
            color = COLORS["selected"] if obj.key == selected else COLORS[obj.kind]
            path = [to_canvas(*obj.points[k]) for k in sorted(obj.points)]
            if len(path) > 1:  # the whole track, faint
                draw.line(path, fill=color + (110,), width=max(1, width // 2))
            if frame not in obj.points:
                continue
            u, v = to_canvas(*obj.points[frame])
            if not (-radius < u < canvas.width + radius and -radius < v < canvas.height + radius):
                continue
            bbox = (u - radius, v - radius, u + radius, v + radius)
            if obj.kind == CANDIDATE:
                draw.rectangle(bbox, outline=color, width=width)
            elif obj.kind == KNOWN:
                draw.ellipse(bbox, outline=color, width=width)
            else:
                for start in range(0, 360, 30):
                    draw.arc(bbox, start, start + 15, fill=color, width=width)
            text = obj.label
            tx, ty = u + radius + 4, v - radius
            tw = draw.textlength(text, font=font)
            if tx + tw > canvas.width - 2:
                tx = u - radius - 4 - tw
            draw.rectangle((tx - 2, ty - 1, tx + tw + 2, ty + font.size + 3), fill=(29, 32, 33, 170))
            draw.text((tx, ty), text, fill=color, font=font)


def _span(obj: SkyObject) -> float:
    pts = list(obj.points.values())
    return max((math.dist(a, b) for a in pts for b in pts), default=0.0)
