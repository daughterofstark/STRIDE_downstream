"""Stage-S3 orchestration.

Thin composition of the reusable subpackages: load the S0 STRIDE table → build
the four hierarchy-reduction tables → validate the reduction layer →
(optionally) write artifacts. No figures, no cross-serotype tests (that is S5),
no calibrated pass/fail claims (the gate is uncalibrated, §0.1).

Public entry points:

- :func:`build_s3` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s3`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_chain_contrast,
    build_monotonicity_audit,
    build_resolution_gap,
    build_scale_curve,
)
from .io import (
    file_digest,
    load_stride_table,
    write_hierarchy_summary,
    write_tables,
)
from .models import S3Report
from .models.schema import PROVISIONAL_RHO_STAR
from .validation import (
    validate_chain_contrast_totals,
    validate_gap_consistency,
    validate_monotonicity_audit_consistency,
    validate_scale_curve_completeness,
    validate_unique_keys,
)


@dataclass
class S3Tables:
    """The four S3 hierarchy-reduction tables, returned together."""

    scale_curve: pd.DataFrame
    resolution_gap: pd.DataFrame
    monotonicity_audit: pd.DataFrame
    chain_contrast: pd.DataFrame


def build_s3(
    stride_table_path: str | Path,
) -> tuple[S3Tables, S3Report]:
    """Load the STRIDE table, build the four reduction tables, and validate them.

    Parameters
    ----------
    stride_table_path
        Path to the S0 ``stride_table.parquet`` (the profile).

    Returns
    -------
    (tables, report)
        The four reduction tables and the run report (with provenance header).

    Raises
    ------
    S3Error
        Any subclass on the first input or consistency problem.
    """
    report = S3Report()

    stride_table = load_stride_table(stride_table_path)
    report.serotypes = sorted(
        stride_table["serotype"].astype(str).unique().tolist()
    )
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": PROVISIONAL_RHO_STAR,
        "n_replicates_note": (
            "K per serotype is a Level-1 fact; at K=3 only domain-scale and "
            "coarser claims are licensed, residue-scale is exploratory"
        ),
        "rho_star_band_note": (
            "S3 hierarchy products (scale curves, Δρ, monotonicity, chain "
            "contrast) are ρ*-independent descriptions of the profile; the "
            "provisional ρ* is used only to label the profile's gated scale and "
            "the distributed-effect flag"
        ),
        "inputs": {
            "stride_table": {
                "path": str(stride_table_path),
                "sha256": file_digest(stride_table_path),
            },
        },
    }
    report.add(
        "inputs loaded",
        "global",
        True,
        f"stride_table={len(stride_table)} rows",
    )

    # -- build --------------------------------------------------------------
    scale_curve = build_scale_curve(stride_table)
    resolution_gap = build_resolution_gap(stride_table)
    monotonicity_audit = build_monotonicity_audit(stride_table)
    chain_contrast = build_chain_contrast(stride_table)

    # -- validate (structural / arithmetic only) ----------------------------
    validate_unique_keys(
        scale_curve,
        resolution_gap,
        monotonicity_audit,
        chain_contrast,
        report,
    )
    validate_scale_curve_completeness(scale_curve, resolution_gap, report)
    validate_gap_consistency(resolution_gap, report)
    validate_monotonicity_audit_consistency(monotonicity_audit, report)
    validate_chain_contrast_totals(chain_contrast, report)

    # -- report facts -------------------------------------------------------
    report.n_scale_curve = int(len(scale_curve))
    report.n_resolution_gap = int(len(resolution_gap))
    report.n_monotonicity_audit = int(len(monotonicity_audit))
    report.n_chain_contrast = int(len(chain_contrast))
    report.facts = _facts(resolution_gap, monotonicity_audit)
    report.add(
        "s3 hierarchy-reduction tables built",
        "global",
        True,
        f"{report.n_scale_curve} scale-curve, "
        f"{report.n_resolution_gap} gap, "
        f"{report.n_monotonicity_audit} monotonicity, "
        f"{report.n_chain_contrast} chain-contrast rows",
    )

    tables = S3Tables(
        scale_curve=scale_curve,
        resolution_gap=resolution_gap,
        monotonicity_audit=monotonicity_audit,
        chain_contrast=chain_contrast,
    )
    return tables, report


def run_s3(
    stride_table_path: str | Path,
    output_dir: str | Path,
) -> tuple[S3Tables, S3Report]:
    """Full S3: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``scale_curve.parquet``, ``resolution_gap.parquet``,
    ``monotonicity_audit.parquet``, ``chain_contrast.parquet`` and
    ``hierarchy_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s3(stride_table_path)
    write_tables(
        tables.scale_curve,
        tables.resolution_gap,
        tables.monotonicity_audit,
        tables.chain_contrast,
        output_dir,
    )
    write_hierarchy_summary(report, output_dir)
    return tables, report


def _facts(
    resolution_gap: pd.DataFrame, monotonicity_audit: pd.DataFrame
) -> dict[str, object]:
    """Compact hierarchy facts for the summary JSON."""
    facts: dict[str, object] = {}
    if not resolution_gap.empty:
        facts["n_distributed_effects"] = int(
            resolution_gap["is_distributed"].sum()
        )
    if not monotonicity_audit.empty:
        facts["n_monotone_loci"] = int(monotonicity_audit["is_monotone"].sum())
        facts["n_non_monotone_loci"] = int(
            (~monotonicity_audit["is_monotone"]).sum()
        )
    return facts
