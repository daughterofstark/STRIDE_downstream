r"""Stage-S7 orchestration (the reporting layer).

Thin composition of the reusable subpackages: resolve the S2–S6 output directories
→ load the reduction tables (never raw STRIDE) → prepare plotting-ready figure data
(F1–F8) and assemble manuscript tables (T1–T5) → render the figures to
deterministic SVG → write artifacts → validate structurally → write the summary.
It computes **no new statistics**: every value is carried through from a prior
stage. Figures are emitted as dependency-free SVG (plus prepared-data CSV/Parquet);
manuscript tables as CSV/Parquet/Markdown. The replicate layer (S6) feeds no design
figure/table, so its blocked-analysis ledger is surfaced under ``limitations`` and
its inputs are digested for provenance — documented, not dropped.

Public entry points:

- :func:`build_s7` — resolve + load + build + in-memory validate; returns the
  prepared artifacts + report, no file writes.
- :func:`run_s7`   — additionally renders/writes artifacts and runs the on-disk
  validation, then writes ``report_summary.json``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .build import build_all_figures, build_all_tables
from .io import (
    file_digest,
    load_inputs,
    resolve_input_paths,
    write_dataframe,
    write_manifest,
    write_markdown_table,
    write_summary,
    write_svg,
)
from .io.loaders import INPUT_REQUIRED_COLUMNS
from .models import ArtifactRecord, S7Report
from .models.schema import (
    CSV_SUFFIX,
    FIGURE_IDS,
    FIGURE_SLUGS,
    FIGURE_SOURCES,
    FIGURE_TITLES,
    IN_REPLICATE_BLOCKED_ANALYSES,
    IN_REPLICATE_REGIME,
    IN_SEROTYPE_SUMMARY,
    MD_SUFFIX,
    N_REPLICATES,
    PARQUET_SUFFIX,
    PROVISIONAL_RHO_STAR,
    SVG_SUFFIX,
    TABLE_IDS,
    TABLE_SLUGS,
    TABLE_SOURCES,
    TABLE_TITLES,
)
from .plotting import render_figure
from .validation import validate_on_disk
from .validation.checks import (
    validate_columns,
    validate_completeness,
    validate_filenames,
    validate_provenance,
)

#: default per-stage input directories (relative), mirroring prior stages' CLIs
DEFAULT_STAGE_DIRS: dict[str, str] = {
    "s2": "outputs_s2",
    "s3": "outputs_s3",
    "s4": "outputs_s4",
    "s5": "outputs_s5",
    "s6": "outputs_s6",
}


@dataclass
class S7Artifacts:
    """The prepared, plotting-ready figure data, rendered SVG, and tables."""

    figure_data: dict[str, pd.DataFrame] = field(default_factory=dict)
    figure_svg: dict[str, str] = field(default_factory=dict)
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)


def _dominant_tier(df: pd.DataFrame) -> str:
    if "tier" not in df.columns or df.empty:
        return ""
    counts = df["tier"].astype(str).value_counts()
    return str(counts.index[0]) if not counts.empty else ""


def build_s7(stage_dirs: dict[str, Path]) -> tuple[S7Artifacts, S7Report]:
    """Load the S2–S6 inputs, build the artifacts, and validate in memory.

    Parameters
    ----------
    stage_dirs
        Mapping of stage key (``"s2"``..``"s6"``) → that stage's output directory.

    Returns
    -------
    tuple[S7Artifacts, S7Report]
        The prepared figure data / rendered SVG / manuscript tables, and the run
        report (provenance, artifact records, checks). No files are written.

    Raises
    ------
    stride_s7.models.errors.InputError
        If a required prior-stage table is missing or malformed.
    stride_s7.models.errors.ConsistencyError
        If a structural invariant fails.
    """
    paths = resolve_input_paths(stage_dirs)
    frames = load_inputs(paths)

    figure_data = build_all_figures(frames)
    tables = build_all_tables(frames)
    figure_svg = {
        fid: render_figure(fid, figure_data[fid], FIGURE_TITLES[fid])
        for fid in FIGURE_IDS
    }
    artifacts = S7Artifacts(figure_data=figure_data, figure_svg=figure_svg, tables=tables)

    report = S7Report()
    report.serotypes = sorted(
        frames[IN_SEROTYPE_SUMMARY]["serotype"].astype(str).unique().tolist()
    )

    # ---- artifact records (files are the deterministic slug names) -----------
    for fid in FIGURE_IDS:
        slug = FIGURE_SLUGS[fid]
        report.figures.append(
            ArtifactRecord(
                artifact_id=fid,
                kind="figure",
                title=FIGURE_TITLES[fid],
                slug=slug,
                sources=list(FIGURE_SOURCES[fid]),
                files=[slug + SVG_SUFFIX, slug + CSV_SUFFIX, slug + PARQUET_SUFFIX],
                n_rows=int(len(figure_data[fid])),
                tier=_dominant_tier(figure_data[fid]),
            )
        )
    for tid in TABLE_IDS:
        slug = TABLE_SLUGS[tid]
        report.tables.append(
            ArtifactRecord(
                artifact_id=tid,
                kind="table",
                title=TABLE_TITLES[tid],
                slug=slug,
                sources=list(TABLE_SOURCES[tid]),
                files=[slug + CSV_SUFFIX, slug + PARQUET_SUFFIX, slug + MD_SUFFIX],
                n_rows=int(len(tables[tid])),
                tier=_dominant_tier(tables[tid]),
            )
        )

    # ---- provenance header (design §5.4) ------------------------------------
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": PROVISIONAL_RHO_STAR,
        "n_replicates_note": (
            f"Figures and tables inherit the uncalibrated gate (calibrated=false, "
            f"provisional rho*={PROVISIONAL_RHO_STAR}) and the K={N_REPLICATES} "
            "replicate regime from their source stages; S7 makes no calibrated "
            "pass/fail claim and recomputes nothing."
        ),
        "figures_format": "svg",
        "figures_format_note": (
            "Figures are emitted as deterministic, dependency-free SVG (vector) "
            "plus their prepared-data CSV/Parquet. PNG/PDF rasterisation is "
            "intentionally not produced: it would require a heavyweight, "
            "version-sensitive rendering dependency whose bytes are not "
            "reproducible, breaking the framework's determinism guarantee. This is "
            "a format choice, not a data-availability limitation."
        ),
        "determinism_note": (
            "Outputs are byte-reproducible: rows are sorted, floats are rounded, "
            "and no wall-clock timestamp is embedded (the design's provenance "
            "'date' field is intentionally omitted to preserve reproducibility)."
        ),
        "inputs": {
            name: {
                "path": str(paths[name]),
                "sha256": file_digest(paths[name]),
            }
            for name in INPUT_REQUIRED_COLUMNS
        },
    }

    # ---- limitations: the replicate layer feeds no design figure/table -------
    report.limitations = _replicate_limitations(
        frames[IN_REPLICATE_REGIME], frames[IN_REPLICATE_BLOCKED_ANALYSES]
    )

    report.facts = _facts(figure_data, tables, report)

    # ---- in-memory structural validation ------------------------------------
    validate_completeness(report)
    validate_columns(figure_data, tables, report)
    validate_filenames(report)
    validate_provenance(report)

    return artifacts, report


def run_s7(
    stage_dirs: dict[str, Path],
    output_dir: str | Path,
) -> tuple[S7Artifacts, S7Report]:
    """Build, then render/write every artifact and run the on-disk validation.

    Writes, per figure: ``<slug>.svg`` + ``<slug>.csv`` + ``<slug>.parquet``; per
    table: ``<slug>.csv`` + ``<slug>.parquet`` + ``<slug>.md``; plus
    ``artifact_manifest.parquet`` and ``report_summary.json``.
    """
    artifacts, report = build_s7(stage_dirs)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for fid in FIGURE_IDS:
        slug = FIGURE_SLUGS[fid]
        write_dataframe(artifacts.figure_data[fid], out, slug)
        write_svg(artifacts.figure_svg[fid], out, slug)
    for tid in TABLE_IDS:
        slug = TABLE_SLUGS[tid]
        write_dataframe(artifacts.tables[tid], out, slug)
        write_markdown_table(artifacts.tables[tid], out, slug, TABLE_TITLES[tid])

    write_manifest(report, out)
    validate_on_disk(report, out)
    write_summary(report, out)
    return artifacts, report


def _replicate_limitations(
    regime: pd.DataFrame, blocked: pd.DataFrame
) -> list[dict[str, object]]:
    """Assemble the replicate-layer status/limitations from S6 (pure passthrough)."""
    limitations: list[dict[str, object]] = []
    for row in blocked.to_dict("records"):
        limitations.append(
            {
                "source": "S6",
                "analysis_id": str(row.get("analysis_id", "")),
                "status": str(row.get("status", "")),
                "available": bool(row.get("available", False)),
                "reason": str(row.get("reason", "")),
            }
        )
    if not regime.empty and "per_replicate_effects_available" in regime.columns:
        any_effects = bool(regime["per_replicate_effects_available"].astype(bool).any())
        limitations.append(
            {
                "source": "S6",
                "analysis_id": "replicate_regime",
                "status": "reported",
                "available": any_effects,
                "reason": (
                    "No design figure/table draws on the replicate (per-run) axis; "
                    "the replicate regime is recorded here for provenance. "
                    f"per_replicate_effects_available={any_effects}."
                ),
            }
        )
    return limitations


def _facts(
    figure_data: dict[str, pd.DataFrame],
    tables: dict[str, pd.DataFrame],
    report: S7Report,
) -> dict[str, object]:
    facts: dict[str, object] = {
        "n_figures": len(figure_data),
        "n_tables": len(tables),
        "n_serotypes": len(report.serotypes),
        "figure_rows": {fid: int(len(df)) for fid, df in figure_data.items()},
        "table_rows": {tid: int(len(df)) for tid, df in tables.items()},
    }
    return facts
