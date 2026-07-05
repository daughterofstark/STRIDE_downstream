"""Stage-S5 orchestration.

Thin composition of the reusable subpackages: load the S0 STRIDE table and the
S1A conservation table → build the four cross-serotype tables → validate the
layer → (optionally) write artifacts. No figures, no calibrated pass/fail claims
(the gate is uncalibrated, §0.1), and no test that treats residues as biological
replicates (serotype is the unit of replication at n = 4, §5.2).

Public entry points:

- :func:`build_s5` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s5`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_cross_serotype_scorecard,
    build_direction_concordance,
    build_domain_serotype_matrix,
    build_position_conservation,
)
from .io import (
    file_digest,
    load_conservation_table,
    load_stride_table,
    write_conservation_summary,
    write_tables,
)
from .models import S5Report
from .models.schema import PROVISIONAL_RHO_STAR
from .validation import (
    validate_cross_serotype_scorecard,
    validate_direction_concordance,
    validate_domain_serotype_matrix,
    validate_position_conservation,
    validate_unique_keys,
)


@dataclass
class S5Tables:
    """The four S5 cross-serotype tables, returned together."""

    position_conservation: pd.DataFrame
    direction_concordance: pd.DataFrame
    domain_serotype_matrix: pd.DataFrame
    cross_serotype_scorecard: pd.DataFrame


def build_s5(
    stride_table_path: str | Path,
    conservation_table_path: str | Path,
    rho_star: float = PROVISIONAL_RHO_STAR,
) -> tuple[S5Tables, S5Report]:
    """Load the inputs, build the four cross-serotype tables, and validate.

    Parameters
    ----------
    stride_table_path
        Path to the S0 ``stride_table.parquet`` (the tidy profile).
    conservation_table_path
        Path to the S1A ``conservation_table.parquet`` (the shared-position index).
    rho_star
        The provisional gate threshold used for the descriptive reproducibility
        determination (default :data:`PROVISIONAL_RHO_STAR`).

    Returns
    -------
    (tables, report)
        The four cross-serotype tables and the run report (with provenance
        header).

    Raises
    ------
    S5Error
        Any subclass on an input or consistency problem.
    """
    report = S5Report()

    stride_table = load_stride_table(stride_table_path)
    conservation_table = load_conservation_table(conservation_table_path)
    report.serotypes = sorted(
        stride_table["serotype"].astype(str).unique().tolist()
    )
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": rho_star,
        "n_replicates_note": (
            "Serotype is the unit of biological replication (n=4). Per-serotype "
            "values are aggregated first, then compared across serotypes; "
            "residues are never treated as independent samples, and results are "
            "descriptive (counts / effect sizes), not p-values across residues."
        ),
        "rho_star_note": (
            "The gate is uncalibrated; 'reproducible' means the residue-scale ρ "
            "clears the provisional ρ* (equivalently, the locus is gated at the "
            "residue scale). Residue-scale products are exploratory; the "
            "domain × serotype matrix is the licensed claim level at K=3."
        ),
        "inputs": {
            "stride_table": {
                "path": str(stride_table_path),
                "sha256": file_digest(stride_table_path),
            },
            "conservation_table": {
                "path": str(conservation_table_path),
                "sha256": file_digest(conservation_table_path),
            },
        },
    }
    report.add(
        "inputs loaded",
        "global",
        True,
        f"stride_table={len(stride_table)} rows, "
        f"conservation_table={len(conservation_table)} rows",
    )

    # -- build --------------------------------------------------------------
    position_conservation = build_position_conservation(
        stride_table, conservation_table, rho_star
    )
    direction_concordance = build_direction_concordance(stride_table, rho_star)
    domain_serotype_matrix = build_domain_serotype_matrix(stride_table)
    cross_serotype_scorecard = build_cross_serotype_scorecard(
        stride_table, conservation_table, rho_star
    )

    # -- validate (structural / arithmetic only) ----------------------------
    validate_unique_keys(
        position_conservation,
        direction_concordance,
        domain_serotype_matrix,
        cross_serotype_scorecard,
        report,
    )
    validate_position_conservation(position_conservation, report)
    validate_direction_concordance(direction_concordance, report)
    validate_domain_serotype_matrix(domain_serotype_matrix, report)
    validate_cross_serotype_scorecard(cross_serotype_scorecard, report)

    # -- report facts -------------------------------------------------------
    report.n_position_conservation = int(len(position_conservation))
    report.n_direction_concordance = int(len(direction_concordance))
    report.n_domain_serotype_matrix = int(len(domain_serotype_matrix))
    report.n_cross_serotype_scorecard = int(len(cross_serotype_scorecard))
    report.facts = _facts(
        position_conservation,
        direction_concordance,
        domain_serotype_matrix,
    )
    report.add(
        "s5 cross-serotype tables built",
        "global",
        True,
        f"{report.n_position_conservation} position-conservation, "
        f"{report.n_direction_concordance} direction-concordance, "
        f"{report.n_domain_serotype_matrix} domain-serotype-matrix, "
        f"{report.n_cross_serotype_scorecard} scorecard rows",
    )

    tables = S5Tables(
        position_conservation=position_conservation,
        direction_concordance=direction_concordance,
        domain_serotype_matrix=domain_serotype_matrix,
        cross_serotype_scorecard=cross_serotype_scorecard,
    )
    return tables, report


def run_s5(
    stride_table_path: str | Path,
    conservation_table_path: str | Path,
    output_dir: str | Path,
    rho_star: float = PROVISIONAL_RHO_STAR,
) -> tuple[S5Tables, S5Report]:
    """Full S5: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``position_conservation.parquet``, ``direction_concordance.parquet``,
    ``domain_serotype_matrix.parquet``, ``cross_serotype_scorecard.parquet`` and
    ``conservation_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s5(
        stride_table_path, conservation_table_path, rho_star
    )
    write_tables(
        tables.position_conservation,
        tables.direction_concordance,
        tables.domain_serotype_matrix,
        tables.cross_serotype_scorecard,
        output_dir,
    )
    write_conservation_summary(report, output_dir)
    return tables, report


def _facts(
    position_conservation: pd.DataFrame,
    direction_concordance: pd.DataFrame,
    domain_serotype_matrix: pd.DataFrame,
) -> dict[str, object]:
    """Compact cross-serotype facts for the summary JSON."""
    facts: dict[str, object] = {}
    if not position_conservation.empty:
        class_counts = (
            position_conservation["conservation_class"].value_counts().to_dict()
        )
        facts["conservation_class_counts"] = {
            str(k): int(v) for k, v in class_counts.items()
        }
        facts["n_serotype_divergent"] = int(
            position_conservation["is_serotype_divergent"].astype(bool).sum()
        )
        facts["n_catalytic_triad_positions"] = int(
            position_conservation["is_catalytic_triad"].astype(bool).sum()
        )
    if not direction_concordance.empty:
        conc_counts = (
            direction_concordance["concordance_class"].value_counts().to_dict()
        )
        facts["concordance_class_counts"] = {
            str(k): int(v) for k, v in conc_counts.items()
        }
    if not domain_serotype_matrix.empty:
        facts["n_domain_serotype_cells"] = int(len(domain_serotype_matrix))
        facts["n_catalytic_domain_cells"] = int(
            domain_serotype_matrix["is_catalytic_domain"].astype(bool).sum()
        )
    return facts
