"""Stage-S4 orchestration.

Thin composition of the reusable subpackages: load the S0 STRIDE table → build
the four uncertainty-layer tables → validate the layer → (optionally) write
artifacts. No figures, no cross-serotype tests (that is S5), no calibrated
pass/fail claims (the gate is uncalibrated, §0.1).

Public entry points:

- :func:`build_s4` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s4`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_domain_effect_summary,
    build_residue_variance,
    build_significance_screen,
    build_variance_budget,
)
from .io import (
    file_digest,
    load_stride_table,
    write_tables,
    write_uncertainty_summary,
)
from .models import S4Report
from .models.schema import GATE_ALPHA, PROVISIONAL_RHO_STAR
from .validation import (
    validate_domain_effect_totals,
    validate_significance_screen,
    validate_unique_keys,
    validate_variance_fractions,
)


@dataclass
class S4Tables:
    """The four S4 uncertainty-layer tables, returned together."""

    variance_budget: pd.DataFrame
    residue_variance: pd.DataFrame
    significance_screen: pd.DataFrame
    domain_effect_summary: pd.DataFrame


def build_s4(stride_table_path: str | Path) -> tuple[S4Tables, S4Report]:
    """Load the STRIDE table, build the four uncertainty tables, and validate.

    Parameters
    ----------
    stride_table_path
        Path to the S0 ``stride_table.parquet`` (the tidy profile).

    Returns
    -------
    (tables, report)
        The four uncertainty tables and the run report (with provenance header).

    Raises
    ------
    S4Error
        Any subclass on the first input or consistency problem.
    """
    report = S4Report()

    stride_table = load_stride_table(stride_table_path)
    report.serotypes = sorted(
        stride_table["serotype"].astype(str).unique().tolist()
    )
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": PROVISIONAL_RHO_STAR,
        "fdr_alpha": GATE_ALPHA,
        "n_replicates_note": (
            "K per serotype is a Level-1 fact; τ² and σ̄² summarise the K=3 "
            "replicates but do not expose per-run values. Serotype is the unit "
            "of replication and the FDR family."
        ),
        "rho_star_band_note": (
            "S4 uncertainty products (variance budgets, τ² ranking, CI screen, "
            "β_se-weighted summaries) are ρ*-independent descriptions of the "
            "profile and the emitted mechanisms; the provisional ρ* only labels "
            "the profile's gated scale."
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
    variance_budget = build_variance_budget(stride_table)
    residue_variance = build_residue_variance(stride_table)
    significance_screen = build_significance_screen(stride_table)
    domain_effect_summary = build_domain_effect_summary(significance_screen)

    # -- validate (structural / arithmetic only) ----------------------------
    validate_unique_keys(
        variance_budget,
        residue_variance,
        significance_screen,
        domain_effect_summary,
        report,
    )
    validate_variance_fractions(variance_budget, residue_variance, report)
    validate_significance_screen(significance_screen, report)
    validate_domain_effect_totals(domain_effect_summary, report)

    # -- report facts -------------------------------------------------------
    report.n_variance_budget = int(len(variance_budget))
    report.n_residue_variance = int(len(residue_variance))
    report.n_significance_screen = int(len(significance_screen))
    report.n_domain_effect_summary = int(len(domain_effect_summary))
    report.facts = _facts(
        variance_budget, significance_screen, domain_effect_summary
    )
    report.add(
        "s4 uncertainty-layer tables built",
        "global",
        True,
        f"{report.n_variance_budget} variance-budget, "
        f"{report.n_residue_variance} residue-variance, "
        f"{report.n_significance_screen} significance-screen, "
        f"{report.n_domain_effect_summary} domain-effect rows",
    )

    tables = S4Tables(
        variance_budget=variance_budget,
        residue_variance=residue_variance,
        significance_screen=significance_screen,
        domain_effect_summary=domain_effect_summary,
    )
    return tables, report


def run_s4(
    stride_table_path: str | Path,
    output_dir: str | Path,
) -> tuple[S4Tables, S4Report]:
    """Full S4: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``variance_budget.parquet``, ``residue_variance.parquet``,
    ``significance_screen.parquet``, ``domain_effect_summary.parquet`` and
    ``uncertainty_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s4(stride_table_path)
    write_tables(
        tables.variance_budget,
        tables.residue_variance,
        tables.significance_screen,
        tables.domain_effect_summary,
        output_dir,
    )
    write_uncertainty_summary(report, output_dir)
    return tables, report


def _facts(
    variance_budget: pd.DataFrame,
    significance_screen: pd.DataFrame,
    domain_effect_summary: pd.DataFrame,
) -> dict[str, object]:
    """Compact uncertainty facts for the summary JSON."""
    facts: dict[str, object] = {}
    if not variance_budget.empty:
        regime_counts = (
            variance_budget["variance_regime"].value_counts().to_dict()
        )
        facts["domain_variance_regime_counts"] = {
            str(k): int(v) for k, v in regime_counts.items()
        }
    if not significance_screen.empty:
        facts["n_signed_mechanisms"] = int(
            significance_screen["is_signed"].astype(bool).sum()
        )
        facts["n_ci_excludes_zero"] = int(
            significance_screen["ci_excludes_zero"].astype(bool).sum()
        )
        facts["n_significant_fdr"] = int(
            significance_screen["significant_fdr"].astype(bool).sum()
        )
    if not domain_effect_summary.empty:
        facts["n_domains_with_mechanisms"] = int(len(domain_effect_summary))
    return facts
