"""Tests for the SVG plotting layer (render_f1..render_f8)."""
from __future__ import annotations

import pandas as pd

from stride_s7.build import build_all_figures
from stride_s7.models.schema import FIGURE_IDS, FIGURE_TITLES
from stride_s7.plotting import render_figure
from stride_s7.plotting.svg import LinearScale, escape, fmt, fmt_num
from tests.s7.fixtures import make_empty_inputs


def test_every_renderer_produces_svg(inputs: dict[str, pd.DataFrame]) -> None:
    figs = build_all_figures(inputs)
    for fid in FIGURE_IDS:
        svg = render_figure(fid, figs[fid], FIGURE_TITLES[fid])
        assert svg.startswith("<svg")
        assert svg.rstrip().endswith("</svg>")
        assert "viewBox" in svg


def test_rendering_is_deterministic(inputs: dict[str, pd.DataFrame]) -> None:
    figs = build_all_figures(inputs)
    for fid in FIGURE_IDS:
        first = render_figure(fid, figs[fid], FIGURE_TITLES[fid])
        second = render_figure(fid, figs[fid], FIGURE_TITLES[fid])
        assert first == second


def test_renderers_are_empty_safe() -> None:
    figs = build_all_figures(make_empty_inputs())
    for fid in FIGURE_IDS:
        svg = render_figure(fid, figs[fid], FIGURE_TITLES[fid])
        assert svg.startswith("<svg")
        assert "</svg>" in svg


def test_title_is_escaped_into_svg(inputs: dict[str, pd.DataFrame]) -> None:
    figs = build_all_figures(inputs)
    svg = render_figure("F3", figs["F3"], "A & B <heatmap>")
    assert "&amp;" in svg
    assert "&lt;heatmap&gt;" in svg
    assert "<heatmap>" not in svg


def test_fmt_normalises_negative_zero() -> None:
    assert fmt(-0.0) == "0.00"
    assert fmt(1.005) in {"1.00", "1.01"}  # round-half handling, still deterministic
    assert fmt_num(float("nan")) == "nan"


def test_escape_covers_xml_specials() -> None:
    assert escape("<a>&'\"") == "&lt;a&gt;&amp;&apos;&quot;"


def test_linear_scale_degenerate_domain() -> None:
    scale = LinearScale(1.0, 1.0, 0.0, 100.0)
    assert scale.to_px(1.0) == 50.0
