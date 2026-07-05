"""Writers for the S7 output artifacts.

Figures are written as prepared-data ``.csv`` + ``.parquet`` (the plotting-ready
table) alongside the rendered ``.svg``. Manuscript tables are written as ``.csv`` +
``.parquet`` + a GitHub-flavoured Markdown ``.md``. A machine-readable
``report_summary.json`` and an ``artifact_manifest.parquet`` record provenance and
every generated file. All writers are deterministic (fixed float formatting, no
timestamps).
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S7Report
from ..models.schema import (
    CSV_SUFFIX,
    MD_SUFFIX,
    OUT_MANIFEST,
    OUT_SUMMARY,
    PARQUET_SUFFIX,
    ROUND_DECIMALS,
    SVG_SUFFIX,
)


def _format_cell(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if value != value:  # NaN
            return ""
        return f"{round(value, ROUND_DECIMALS):g}"
    return str(value)


def write_dataframe(df: pd.DataFrame, output_dir: str | Path, slug: str) -> list[str]:
    """Write ``slug.csv`` and ``slug.parquet`` deterministically; return filenames."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_name = slug + CSV_SUFFIX
    pq_name = slug + PARQUET_SUFFIX
    # deterministic CSV: fixed quoting, explicit line terminator, formatted cells
    with (out / csv_name).open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(list(df.columns))
        for row in df.itertuples(index=False, name=None):
            writer.writerow([_format_cell(v) for v in row])
    df.to_parquet(out / pq_name, index=False)
    return [csv_name, pq_name]


def write_svg(svg: str, output_dir: str | Path, slug: str) -> str:
    """Write ``slug.svg``; return the filename."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    name = slug + SVG_SUFFIX
    (out / name).write_text(svg)
    return name


def write_markdown_table(
    df: pd.DataFrame, output_dir: str | Path, slug: str, title: str
) -> str:
    """Write a deterministic GitHub-flavoured Markdown table; return the filename."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    name = slug + MD_SUFFIX
    lines = [f"# {title}", ""]
    cols = list(df.columns)
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join("---" for _ in cols) + " |")
    for row in df.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(_format_cell(v) for v in row) + " |")
    (out / name).write_text("\n".join(lines) + "\n")
    return name


def write_manifest(report: S7Report, output_dir: str | Path) -> str:
    """Write ``artifact_manifest.parquet`` listing every generated artifact."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for rec in [*report.figures, *report.tables]:
        rows.append(
            {
                "artifact_id": rec.artifact_id,
                "kind": rec.kind,
                "title": rec.title,
                "slug": rec.slug,
                "sources": ";".join(rec.sources),
                "files": ";".join(rec.files),
                "n_rows": rec.n_rows,
                "tier": rec.tier,
            }
        )
    manifest = pd.DataFrame(
        rows,
        columns=[
            "artifact_id",
            "kind",
            "title",
            "slug",
            "sources",
            "files",
            "n_rows",
            "tier",
        ],
    )
    manifest.to_parquet(out / OUT_MANIFEST, index=False)
    return OUT_MANIFEST


def write_summary(report: S7Report, output_dir: str | Path) -> Path:
    """Write the machine-readable ``report_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S7",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "provenance": report.provenance,
        "facts": report.facts,
        "figures": [asdict(r) for r in report.figures],
        "tables": [asdict(r) for r in report.tables],
        "limitations": report.limitations,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
