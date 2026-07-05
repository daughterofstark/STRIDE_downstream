"""Structural validation of the S7 reporting layer.

Structural checks only — no scientific conclusions are validated (S7 asserts
nothing about the numbers, only that the right artifacts, columns, filenames, and
provenance are present and on disk). Each check appends a
:class:`~stride_s7.models.ValidationCheck` and raises
:class:`~stride_s7.models.errors.ConsistencyError` on failure.

Checks:

- **completeness** — every design figure (F1–F8) and table (T1–T5) has a record;
- **columns** — each prepared figure-data frame and each manuscript table carries
  exactly its declared schema;
- **filenames** — every recorded output filename is the deterministic
  ``<slug><suffix>`` and every figure has ``.svg`` + ``.csv`` + ``.parquet`` while
  every table has ``.csv`` + ``.parquet`` + ``.md``;
- **on-disk** — every recorded file exists in the output directory;
- **provenance** — the header carries ``calibrated`` (False), the provisional ρ\\*,
  the K note, and a per-input SHA-256 digest.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..models import S7Report
from ..models.errors import ConsistencyError
from ..models.schema import (
    CSV_SUFFIX,
    FIGURE_DATA_COLUMNS,
    FIGURE_IDS,
    FIGURE_SLUGS,
    MD_SUFFIX,
    PARQUET_SUFFIX,
    SVG_SUFFIX,
    TABLE_COLUMNS,
    TABLE_IDS,
    TABLE_SLUGS,
)


def validate_completeness(report: S7Report) -> None:
    """Every design figure and table is present exactly once."""
    fig_ids = [r.artifact_id for r in report.figures]
    tab_ids = [r.artifact_id for r in report.tables]
    missing_f = [f for f in FIGURE_IDS if f not in fig_ids]
    missing_t = [t for t in TABLE_IDS if t not in tab_ids]
    if missing_f:
        raise ConsistencyError(f"missing figure record(s): {missing_f}")
    if missing_t:
        raise ConsistencyError(f"missing table record(s): {missing_t}")
    if len(fig_ids) != len(set(fig_ids)):
        raise ConsistencyError(f"duplicate figure records: {fig_ids}")
    if len(tab_ids) != len(set(tab_ids)):
        raise ConsistencyError(f"duplicate table records: {tab_ids}")
    report.add(
        "completeness",
        "artifacts",
        True,
        f"{len(fig_ids)} figures + {len(tab_ids)} tables recorded",
    )


def validate_columns(
    figure_data: dict[str, pd.DataFrame],
    tables: dict[str, pd.DataFrame],
    report: S7Report,
) -> None:
    """Prepared figure frames and manuscript tables carry their declared schema."""
    for fid in FIGURE_IDS:
        expected = list(FIGURE_DATA_COLUMNS[fid])
        got = list(figure_data[fid].columns)
        if got != expected:
            raise ConsistencyError(
                f"figure {fid} data columns {got} != expected {expected}"
            )
    for tid in TABLE_IDS:
        expected = list(TABLE_COLUMNS[tid])
        got = list(tables[tid].columns)
        if got != expected:
            raise ConsistencyError(
                f"table {tid} columns {got} != expected {expected}"
            )
    report.add(
        "columns",
        "artifacts",
        True,
        "every figure-data frame and table has its declared schema",
    )


def validate_filenames(report: S7Report) -> None:
    """Recorded filenames are the deterministic slug names with the right suffixes."""
    for rec in report.figures:
        slug = FIGURE_SLUGS[rec.artifact_id]
        expected = {slug + SVG_SUFFIX, slug + CSV_SUFFIX, slug + PARQUET_SUFFIX}
        if set(rec.files) != expected or rec.slug != slug:
            raise ConsistencyError(
                f"figure {rec.artifact_id} files {rec.files} != {sorted(expected)}"
            )
    for rec in report.tables:
        slug = TABLE_SLUGS[rec.artifact_id]
        expected = {slug + CSV_SUFFIX, slug + PARQUET_SUFFIX, slug + MD_SUFFIX}
        if set(rec.files) != expected or rec.slug != slug:
            raise ConsistencyError(
                f"table {rec.artifact_id} files {rec.files} != {sorted(expected)}"
            )
    report.add(
        "filenames",
        "artifacts",
        True,
        "all artifact filenames are deterministic slug names",
    )


def validate_on_disk(report: S7Report, output_dir: str | Path) -> None:
    """Every recorded artifact file exists in ``output_dir``."""
    out = Path(output_dir)
    missing = []
    for rec in [*report.figures, *report.tables]:
        for name in rec.files:
            if not (out / name).is_file():
                missing.append(name)
    if missing:
        raise ConsistencyError(f"recorded artifact(s) not on disk: {missing}")
    report.add(
        "on_disk",
        "artifacts",
        True,
        f"{len(report.figures) + len(report.tables)} artifacts present on disk",
    )


def validate_provenance(report: S7Report) -> None:
    """The provenance header is complete."""
    prov = report.provenance
    for key in ("calibrated", "provisional_rho_star", "n_replicates_note", "inputs"):
        if key not in prov:
            raise ConsistencyError(f"provenance missing '{key}'")
    if prov.get("calibrated") is not False:
        raise ConsistencyError("provenance 'calibrated' must be False")
    inputs = prov.get("inputs", {})
    if not isinstance(inputs, dict) or not inputs:
        raise ConsistencyError("provenance 'inputs' must be a non-empty mapping")
    for name, meta in inputs.items():
        if "sha256" not in meta:
            raise ConsistencyError(f"provenance input '{name}' missing sha256")
    report.add(
        "provenance",
        "summary",
        True,
        f"provenance complete over {len(inputs)} inputs, calibrated=false",
    )


def validate_all(
    figure_data: dict[str, pd.DataFrame],
    tables: dict[str, pd.DataFrame],
    report: S7Report,
    output_dir: str | Path,
) -> None:
    """Run every structural check, populating ``report``."""
    validate_completeness(report)
    validate_columns(figure_data, tables, report)
    validate_filenames(report)
    validate_on_disk(report, output_dir)
    validate_provenance(report)
