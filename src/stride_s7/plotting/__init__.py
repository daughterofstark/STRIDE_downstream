"""S7 plotting subpackage: deterministic, dependency-free SVG rendering.

The renderers consume the plotting-ready frames produced by :mod:`stride_s7.build`
and emit SVG strings. No plotting library is used; output is byte-reproducible.
"""
from __future__ import annotations

from .figures import FIGURE_RENDERERS, render_figure
from .svg import Chart, LinearScale

__all__ = [
    "FIGURE_RENDERERS",
    "render_figure",
    "Chart",
    "LinearScale",
]
