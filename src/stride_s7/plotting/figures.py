r"""Render the prepared figure data (F1–F8) to deterministic SVG strings.

Each ``render_fN`` consumes the plotting-ready frame produced by the matching
``prepare_fN`` builder and returns an SVG string. Renderers never read a source
table and compute nothing — they lay out marks from already-prepared values. Empty
inputs render a valid, titled placeholder rather than failing.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    CONSERVATION_CLASS_ORDER,
    SVG_CATEGORICAL_PALETTE,
    SVG_MARK_COLOR,
    SVG_MARK_COLOR_ALT,
    SVG_SEQUENTIAL_RAMP,
)
from .svg import Chart, LinearScale, circle, line, polyline, rect, text


def _palette_for(categories: list[str]) -> dict[str, str]:
    """Deterministic category → colour by sorted index."""
    out: dict[str, str] = {}
    for i, cat in enumerate(sorted(dict.fromkeys(categories))):
        out[cat] = SVG_CATEGORICAL_PALETTE[i % len(SVG_CATEGORICAL_PALETTE)]
    return out


def _legend(chart: Chart, mapping: dict[str, str], x: float, y0: float) -> None:
    """A compact swatch legend at ``(x, y0)`` descending."""
    for i, (label, color) in enumerate(mapping.items()):
        yy = y0 + i * 16
        chart.add(rect(x, yy - 9, 10, 10, color))
        chart.add(text(x + 14, yy, label, size=11))


def _ramp_color(value: float, lo: float, hi: float) -> str:
    if hi <= lo or value != value:
        return SVG_SEQUENTIAL_RAMP[0]
    frac = (value - lo) / (hi - lo)
    idx = int(round(frac * (len(SVG_SEQUENTIAL_RAMP) - 1)))
    idx = max(0, min(len(SVG_SEQUENTIAL_RAMP) - 1, idx))
    return SVG_SEQUENTIAL_RAMP[idx]


def render_f1(df: pd.DataFrame, title: str) -> str:
    """Reproducibility landscape: ρ vs position, one polyline per serotype."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("canonical position (ordered)")
    chart.y_axis_label("rho (residue)")
    if df.empty:
        chart.note("no residue-landscape rows")
        return chart.render()
    positions = list(dict.fromkeys(df["canon_label"].tolist()))
    xs = {p: i for i, p in enumerate(positions)}
    xscale = LinearScale(0, max(1, len(positions) - 1), chart.plot_left + 10, chart.plot_right - 10)
    yscale = LinearScale(0.0, 1.0, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, 0.25, 0.5, 0.75, 1.0])
    serotypes = sorted(df["serotype"].unique())
    palette = _palette_for(serotypes)
    for s in serotypes:
        sub = df[df["serotype"] == s]
        pts = [
            (xscale.to_px(xs[p]), yscale.to_px(float(r)))
            for p, r in zip(sub["canon_label"], sub["rho_residue"], strict=True)
        ]
        chart.add(polyline(pts, palette[s], 1.2))
        for px, py in pts:
            chart.add(circle(px, py, 2.2, palette[s]))
    _legend(chart, palette, chart.plot_right - 90, chart.plot_top + 12)
    return chart.render()


def render_f2(df: pd.DataFrame, title: str) -> str:
    """Achieved-resolution census: stacked bar of loci by gated scale per serotype."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("serotype")
    chart.y_axis_label("n loci (gated)")
    if df.empty:
        chart.note("no resolution-census rows")
        return chart.render()
    serotypes = sorted(df["serotype"].unique())
    levels = [
        lv
        for _, lv in sorted(
            {
                (int(i), s)
                for i, s in zip(
                    df["gated_scale_index"], df["gated_scale_level"], strict=True
                )
            }
        )
    ]
    palette = _palette_for(levels)
    totals = df.groupby("serotype")["n_loci"].sum()
    ymax = float(totals.max()) if not totals.empty else 1.0
    yscale = LinearScale(0.0, ymax, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, ymax / 2.0, ymax], decimals=0)
    n = len(serotypes)
    slot = (chart.plot_right - chart.plot_left) / max(1, n)
    bar_w = slot * 0.5
    centers = []
    for i, s in enumerate(serotypes):
        cx = chart.plot_left + slot * (i + 0.5)
        centers.append(cx)
        acc = 0.0
        sub = df[df["serotype"] == s]
        stack = {lv: 0.0 for lv in levels}
        for lv, cnt in zip(sub["gated_scale_level"], sub["n_loci"], strict=True):
            stack[lv] = stack.get(lv, 0.0) + float(cnt)
        for lv in levels:
            cnt = stack[lv]
            if cnt <= 0:
                continue
            y_top = yscale.to_px(acc + cnt)
            y_bot = yscale.to_px(acc)
            chart.add(rect(cx - bar_w / 2, y_top, bar_w, y_bot - y_top, palette[lv]))
            acc += cnt
    chart.x_category_ticks(serotypes, centers)
    _legend(chart, palette, chart.plot_right - 150, chart.plot_top + 12)
    return chart.render()


def render_f3(df: pd.DataFrame, title: str) -> str:
    """Domain × serotype ρ heatmap."""
    chart = Chart(title)
    chart.draw_frame()
    if df.empty:
        chart.note("no domain-serotype matrix rows")
        return chart.render()
    domains = sorted(df["domain"].unique())
    serotypes = sorted(df["serotype"].unique())
    lo, hi = 0.0, 1.0
    cell_w = (chart.plot_right - chart.plot_left) / max(1, len(serotypes))
    cell_h = (chart.plot_bottom - chart.plot_top) / max(1, len(domains))
    lookup = {
        (str(s), str(d)): float(r)
        for s, d, r in zip(df["serotype"], df["domain"], df["rho_domain"], strict=True)
    }
    for di, d in enumerate(domains):
        cy = chart.plot_top + di * cell_h
        chart.add(text(chart.plot_left - 6, cy + cell_h / 2 + 4, d, anchor="end", size=10))
        for si, s in enumerate(serotypes):
            cx = chart.plot_left + si * cell_w
            val = lookup.get((s, d))
            fill = _ramp_color(val, lo, hi) if val is not None else "#eeeeee"
            chart.add(rect(cx, cy, cell_w, cell_h, fill, stroke="#ffffff", stroke_width=1.0))
            if val is not None:
                tcol = "#ffffff" if val >= 0.6 else "#111111"
                chart.add(
                    text(cx + cell_w / 2, cy + cell_h / 2 + 4, f"{val:.2f}", anchor="middle", size=10, color=tcol)
                )
    for si, s in enumerate(serotypes):
        cx = chart.plot_left + si * cell_w + cell_w / 2
        chart.add(text(cx, chart.plot_bottom + 16, s, anchor="middle", size=10))
    return chart.render()


def render_f4(df: pd.DataFrame, title: str) -> str:
    """Signed-effect forest: β_signed ± CI whiskers for coherent mechanisms."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("beta_signed (+/- CI)")
    if df.empty:
        chart.note("no signed (coherent) mechanisms")
        return chart.render()
    rows = [
        (f"{s}:{c}", float(b), float(lo), float(hi), bool(sig))
        for s, c, b, lo, hi, sig in zip(
            df["serotype"], df["canon_label"], df["beta_signed"],
            df["beta_ci_lower"], df["beta_ci_upper"], df["significant_fdr"], strict=True,
        )
    ]
    lo_all = min([r[2] for r in rows] + [0.0])
    hi_all = max([r[3] for r in rows] + [0.0])
    span = (hi_all - lo_all) or 1.0
    xscale = LinearScale(lo_all - 0.05 * span, hi_all + 0.05 * span, chart.plot_left + 10, chart.plot_right - 10)
    n = len(rows)
    step = (chart.plot_bottom - chart.plot_top) / max(1, n)
    zero_px = xscale.to_px(0.0)
    chart.add(line(zero_px, chart.plot_top, zero_px, chart.plot_bottom, "#999999", 0.8))
    for tick in (xscale.d0, 0.0, xscale.d1):
        px = xscale.to_px(tick)
        chart.add(text(px, chart.plot_bottom + 16, f"{tick:.2f}", anchor="middle", size=10))
    for i, (label, beta, clo, chi, sig) in enumerate(rows):
        cy = chart.plot_top + step * (i + 0.5)
        color = SVG_MARK_COLOR if sig else SVG_MARK_COLOR_ALT
        chart.add(line(xscale.to_px(clo), cy, xscale.to_px(chi), cy, color, 1.2))
        chart.add(circle(xscale.to_px(beta), cy, 2.6, color))
        chart.add(text(chart.plot_left - 6, cy + 3, label, anchor="end", size=9))
    return chart.render()


def render_f5(df: pd.DataFrame, title: str) -> str:
    """Cross-serotype conservation: per-position bar coloured by conservation class."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("shared position (ordered)")
    chart.y_axis_label("frac reproducible")
    if df.empty:
        chart.note("no shared positions")
        return chart.render()
    classes = [c for c in CONSERVATION_CLASS_ORDER if c in set(df["conservation_class"])]
    extra = sorted(set(df["conservation_class"]) - set(CONSERVATION_CLASS_ORDER))
    ordered_classes = classes + extra
    palette = _palette_for(ordered_classes)
    yscale = LinearScale(0.0, 1.0, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, 0.5, 1.0])
    n = len(df)
    slot = (chart.plot_right - chart.plot_left) / max(1, n)
    bar_w = slot * 0.7
    for i, row in enumerate(df.to_dict("records")):
        cx = chart.plot_left + slot * (i + 0.5)
        frac = float(row["frac_reproducible"])
        y_top = yscale.to_px(frac)
        color = palette.get(str(row["conservation_class"]), SVG_MARK_COLOR)
        chart.add(rect(cx - bar_w / 2, y_top, bar_w, chart.plot_bottom - y_top, color))
        if bool(row["is_serotype_divergent"]):
            chart.add(text(cx, y_top - 3, "*", anchor="middle", size=11))
    _legend(chart, palette, chart.plot_right - 175, chart.plot_top + 12)
    return chart.render()


def render_f6(df: pd.DataFrame, title: str) -> str:
    """Variance composition: stacked τ²/σ̄² fractions per (serotype, domain)."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("serotype / domain")
    chart.y_axis_label("variance fraction")
    if df.empty:
        chart.note("no variance-budget rows")
        return chart.render()
    yscale = LinearScale(0.0, 1.0, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, 0.5, 1.0])
    n = len(df)
    slot = (chart.plot_right - chart.plot_left) / max(1, n)
    bar_w = slot * 0.6
    labels, centers = [], []
    for i, row in enumerate(df.to_dict("records")):
        cx = chart.plot_left + slot * (i + 0.5)
        centers.append(cx)
        labels.append(f"{row['serotype']}/{row['domain']}")
        ftau = float(row["frac_tau2"])
        fsig = float(row["frac_sigma2"])
        y_tau_top = yscale.to_px(ftau)
        chart.add(rect(cx - bar_w / 2, y_tau_top, bar_w, chart.plot_bottom - y_tau_top, SVG_MARK_COLOR))
        y_sig_top = yscale.to_px(ftau + fsig)
        chart.add(rect(cx - bar_w / 2, y_sig_top, bar_w, y_tau_top - y_sig_top, SVG_MARK_COLOR_ALT))
    chart.x_category_ticks(labels, centers)
    _legend(chart, {"tau^2": SVG_MARK_COLOR, "sigma-bar^2": SVG_MARK_COLOR_ALT}, chart.plot_right - 130, chart.plot_top + 12)
    return chart.render()


def render_f7(df: pd.DataFrame, title: str) -> str:
    """ρ-vs-scale trajectories for the catalytic regions."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("scale index (0=residue .. 6=complex)")
    chart.y_axis_label("rho")
    if df.empty:
        chart.note("no catalytic-region scale curves")
        return chart.render()
    xscale = LinearScale(0.0, 6.0, chart.plot_left + 10, chart.plot_right - 10)
    yscale = LinearScale(0.0, 1.0, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, 0.5, 1.0])
    for xi in range(7):
        px = xscale.to_px(xi)
        chart.add(text(px, chart.plot_bottom + 16, str(xi), anchor="middle", size=10))
    serotypes = sorted(df["serotype"].unique())
    palette = _palette_for(serotypes)
    for (s, _canon), sub in df.groupby(["serotype", "canon_label"], sort=True):
        pts = [
            (xscale.to_px(int(i)), yscale.to_px(float(r)))
            for i, r in zip(sub["scale_index"], sub["rho"], strict=True)
        ]
        chart.add(polyline(pts, palette[str(s)], 1.2))
    _legend(chart, palette, chart.plot_right - 90, chart.plot_top + 12)
    return chart.render()


def render_f8(df: pd.DataFrame, title: str) -> str:
    """Coherence vs ρ scatter at the domain scale."""
    chart = Chart(title)
    chart.draw_frame()
    chart.x_axis_label("rho (domain)")
    chart.y_axis_label("coherence (domain)")
    if df.empty:
        chart.note("no domain-reproducibility rows")
        return chart.render()
    xscale = LinearScale(0.0, 1.0, chart.plot_left + 10, chart.plot_right - 10)
    yscale = LinearScale(0.0, 1.0, chart.plot_bottom, chart.plot_top)
    chart.y_ticks(yscale, [0.0, 0.5, 1.0])
    for t in (0.0, 0.5, 1.0):
        chart.add(text(xscale.to_px(t), chart.plot_bottom + 16, f"{t:.1f}", anchor="middle", size=10))
    for row in df.to_dict("records"):
        cx = xscale.to_px(float(row["rho_domain"]))
        cy = yscale.to_px(float(row["coherence_domain"]))
        color = SVG_MARK_COLOR if bool(row["is_coherent"]) else SVG_MARK_COLOR_ALT
        chart.add(circle(cx, cy, 3.0, color))
    _legend(chart, {"coherent": SVG_MARK_COLOR, "mixed": SVG_MARK_COLOR_ALT}, chart.plot_right - 110, chart.plot_top + 12)
    return chart.render()


#: id → renderer
FIGURE_RENDERERS = {
    "F1": render_f1,
    "F2": render_f2,
    "F3": render_f3,
    "F4": render_f4,
    "F5": render_f5,
    "F6": render_f6,
    "F7": render_f7,
    "F8": render_f8,
}


def render_figure(figure_id: str, df: pd.DataFrame, title: str) -> str:
    """Render one figure by id from its prepared data."""
    return FIGURE_RENDERERS[figure_id](df, title)
