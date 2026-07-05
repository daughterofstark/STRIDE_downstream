r"""Deterministic, dependency-free SVG primitives and a minimal chart scaffold.

Pure string builders — no plotting library, no randomness, no timestamps — so the
same prepared data always renders byte-for-byte identical SVG. All coordinates are
formatted to two decimals; text is XML-escaped; colours come from the frozen
palette in :mod:`stride_s7.models.schema`. The :class:`Chart` helper lays out a
plot area inside the fixed margins and maps data coordinates to pixels, and draws
axes, ticks, a title, and axis labels.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..models.schema import (
    SVG_AXIS_COLOR,
    SVG_BACKGROUND,
    SVG_FONT_FAMILY,
    SVG_FONT_SIZE,
    SVG_GRID_COLOR,
    SVG_HEIGHT,
    SVG_MARGIN_BOTTOM,
    SVG_MARGIN_LEFT,
    SVG_MARGIN_RIGHT,
    SVG_MARGIN_TOP,
    SVG_TEXT_COLOR,
    SVG_TITLE_FONT_SIZE,
    SVG_WIDTH,
)

_ESCAPE = {
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&apos;",
}


def escape(text: object) -> str:
    """XML-escape a value's string form."""
    out = str(text)
    for ch, rep in _ESCAPE.items():
        out = out.replace(ch, rep)
    return out


def fmt(value: float) -> str:
    """Format a coordinate to two decimals (``-0`` normalised to ``0``)."""
    v = round(float(value), 2)
    if v == 0.0:
        v = 0.0
    return f"{v:.2f}"


def fmt_num(value: float, decimals: int = 2) -> str:
    """Format a data value for a tick/label deterministically."""
    if value != value:  # NaN
        return "nan"
    v = round(float(value), decimals)
    if v == 0.0:
        v = 0.0
    return f"{v:.{decimals}f}"


def rect(
    x: float,
    y: float,
    w: float,
    h: float,
    fill: str,
    stroke: str = "none",
    stroke_width: float = 0.0,
) -> str:
    return (
        f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{fmt(stroke_width)}"/>'
    )


def line(
    x1: float, y1: float, x2: float, y2: float, stroke: str, stroke_width: float = 1.0
) -> str:
    return (
        f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" '
        f'stroke="{stroke}" stroke-width="{fmt(stroke_width)}"/>'
    )


def circle(cx: float, cy: float, r: float, fill: str, stroke: str = "none") -> str:
    return (
        f'<circle cx="{fmt(cx)}" cy="{fmt(cy)}" r="{fmt(r)}" '
        f'fill="{fill}" stroke="{stroke}"/>'
    )


def text(
    x: float,
    y: float,
    content: object,
    *,
    size: int = SVG_FONT_SIZE,
    color: str = SVG_TEXT_COLOR,
    anchor: str = "start",
    rotate: float | None = None,
) -> str:
    transform = (
        f' transform="rotate({fmt(rotate)} {fmt(x)} {fmt(y)})"'
        if rotate is not None
        else ""
    )
    return (
        f'<text x="{fmt(x)}" y="{fmt(y)}" font-family="{SVG_FONT_FAMILY}" '
        f'font-size="{size}" fill="{color}" text-anchor="{anchor}"{transform}>'
        f"{escape(content)}</text>"
    )


def polyline(points: list[tuple[float, float]], stroke: str, stroke_width: float = 1.5) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return (
        f'<polyline points="{pts}" fill="none" stroke="{stroke}" '
        f'stroke-width="{fmt(stroke_width)}"/>'
    )


@dataclass
class LinearScale:
    """Map a data interval ``[d0, d1]`` onto a pixel interval ``[p0, p1]``."""

    d0: float
    d1: float
    p0: float
    p1: float

    def to_px(self, value: float) -> float:
        if self.d1 == self.d0:
            return (self.p0 + self.p1) / 2.0
        frac = (value - self.d0) / (self.d1 - self.d0)
        return self.p0 + frac * (self.p1 - self.p0)


@dataclass
class Chart:
    """A fixed-size SVG chart with a bordered plot area and axis helpers."""

    title: str
    width: int = SVG_WIDTH
    height: int = SVG_HEIGHT
    elements: list[str] = field(default_factory=list)

    @property
    def plot_left(self) -> float:
        return float(SVG_MARGIN_LEFT)

    @property
    def plot_right(self) -> float:
        return float(self.width - SVG_MARGIN_RIGHT)

    @property
    def plot_top(self) -> float:
        return float(SVG_MARGIN_TOP)

    @property
    def plot_bottom(self) -> float:
        return float(self.height - SVG_MARGIN_BOTTOM)

    def add(self, element: str) -> None:
        self.elements.append(element)

    def draw_frame(self) -> None:
        """Background, title, and the plot-area border."""
        self.add(rect(0, 0, self.width, self.height, SVG_BACKGROUND))
        self.add(
            text(
                self.width / 2.0,
                SVG_MARGIN_TOP - 20,
                self.title,
                size=SVG_TITLE_FONT_SIZE,
                anchor="middle",
            )
        )
        self.add(
            rect(
                self.plot_left,
                self.plot_top,
                self.plot_right - self.plot_left,
                self.plot_bottom - self.plot_top,
                "none",
                stroke=SVG_AXIS_COLOR,
                stroke_width=1.0,
            )
        )

    def x_axis_label(self, label: str) -> None:
        self.add(
            text(
                (self.plot_left + self.plot_right) / 2.0,
                self.height - 16,
                label,
                anchor="middle",
            )
        )

    def y_axis_label(self, label: str) -> None:
        self.add(
            text(
                18,
                (self.plot_top + self.plot_bottom) / 2.0,
                label,
                anchor="middle",
                rotate=-90,
            )
        )

    def y_ticks(self, scale: LinearScale, values: list[float], decimals: int = 2) -> None:
        for v in values:
            py = scale.to_px(v)
            self.add(line(self.plot_left - 4, py, self.plot_left, py, SVG_AXIS_COLOR))
            self.add(
                line(self.plot_left, py, self.plot_right, py, SVG_GRID_COLOR, 0.5)
            )
            self.add(
                text(
                    self.plot_left - 8,
                    py + 4,
                    fmt_num(v, decimals),
                    anchor="end",
                )
            )

    def x_category_ticks(self, labels: list[str], centers: list[float]) -> None:
        for lab, cx in zip(labels, centers, strict=True):
            self.add(
                line(cx, self.plot_bottom, cx, self.plot_bottom + 4, SVG_AXIS_COLOR)
            )
            self.add(
                text(
                    cx,
                    self.plot_bottom + 16,
                    lab,
                    anchor="end",
                    rotate=-45,
                )
            )

    def note(self, message: str) -> None:
        """A centred placeholder note (used for empty inputs)."""
        self.add(
            text(
                (self.plot_left + self.plot_right) / 2.0,
                (self.plot_top + self.plot_bottom) / 2.0,
                message,
                anchor="middle",
                color="#888888",
            )
        )

    def render(self) -> str:
        body = "\n".join(self.elements)
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n{body}\n</svg>\n'
        )
