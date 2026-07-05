"""Stage-S6 orchestration.

Thin composition of the reusable subpackages: load the S1A replicate inventory
(required) and the optional S0 replicate table → tidy the per-run effects → build
the four replicate-layer tables → validate the layer → (optionally) write
artifacts. No figures, no calibrated pass/fail claims (the gate is uncalibrated,
§0.1), no τ²-based products (those are S4's aggregate axis), and no reconstruction
of per-run θ when the replicate table is absent (the design's blocked state, §4.1)
— the per-run analyses are simply recorded as blocked.

Public entry points:

- :func:`build_s6` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s6`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_replicate_blocked_analyses,
    build_replicate_concordance,
    build_replicate_effect_spread,
    build_replicate_regime,
    per_replicate_effects_available,
    per_run_effects,
)
from .io import (
    file_digest,
    load_replicate_inventory,
    load_replicate_table,
    write_replicate_summary,
    write_tables,
)
from .models import S6Report
from .validation import validate_all


@dataclass
class S6Tables:
    """The four S6 replicate-layer tables, returned together."""

    replicate_regime: pd.DataFrame
    replicate_effect_spread: pd.DataFrame
    replicate_concordance: pd.DataFrame
    replicate_blocked_analyses: pd.DataFrame


def build_s6(
    replicate_inventory_path: str | Path,
    replicate_table_path: str | Path | None = None,
) -> tuple[S6Tables, S6Report]:
    """Load the inputs, build the four replicate-layer tables, and validate.

    Parameters
    ----------
    replicate_inventory_path
        Path to the S1A ``replicate_inventory.parquet`` (the replicate-structure
        index). Required.
    replicate_table_path
        Path to the S0 ``replicate_table.parquet`` (the Level-1 per-run
        observations). Optional: when ``None`` or absent, the per-run analyses are
        recorded as blocked rather than approximated (design §4.1).

    Returns
    -------
    (tables, report)
        The four replicate-layer tables and the run report (with provenance
        header and the blocked-analysis ledger echoed into the summary).

    Raises
    ------
    S6Error
        Any subclass on an input or consistency problem.
    """
    report = S6Report()

    replicate_inventory = load_replicate_inventory(replicate_inventory_path)
    replicate_table = (
        load_replicate_table(replicate_table_path)
        if replicate_table_path is not None
        else None
    )
    effects = per_run_effects(replicate_table)
    available = per_replicate_effects_available(effects)

    report.serotypes = sorted(
        replicate_inventory["serotype"].astype(str).unique().tolist()
    )
    report.per_replicate_effects_available = bool(available)
    report.provenance = {
        "calibrated": False,
        "n_replicates_note": (
            "Replicates are the sampling unit for the per-run observations "
            "(K=3 in the design dataset). STRIDE licenses residue-scale claims "
            "only at K>=5, so the per-run spread and concordance are exploratory; "
            "per-run θ is read, never recomputed, and is not reconstructible from "
            "the K-aggregate τ²/σ̄²."
        ),
        "per_replicate_effects_note": (
            "The per-run rank concordance and effect spread require the Level-1 "
            "per-run correlation CSVs (ingested by S0 into replicate_table). When "
            "those are absent (the real dengue upload's state, §4.1) the "
            "replicate table is not written and both analyses are recorded as "
            "blocked. Leave-one-replicate-out stability is always blocked (it "
            "needs a STRIDE re-run, §4.2). The τ²-based replicate-disagreement "
            "mapping and τ²/σ̄² regime diagnostic (§3.1) read the aggregate and "
            "are produced by S4, not S6."
        ),
        "per_replicate_effects_available": bool(available),
        "inputs": {
            "replicate_inventory": {
                "path": str(replicate_inventory_path),
                "sha256": file_digest(replicate_inventory_path),
            },
            "replicate_table": {
                "path": (
                    str(replicate_table_path)
                    if replicate_table_path is not None
                    else ""
                ),
                "present": replicate_table is not None,
                "sha256": (
                    file_digest(replicate_table_path)
                    if replicate_table_path is not None
                    else ""
                ),
            },
        },
    }
    report.add(
        "inputs loaded",
        "global",
        True,
        f"replicate_inventory={len(replicate_inventory)} rows, "
        f"replicate_table={'absent' if replicate_table is None else len(replicate_table)}, "
        f"per_replicate_effects_available={available}",
    )

    # -- build --------------------------------------------------------------
    replicate_regime = build_replicate_regime(replicate_inventory, effects)
    replicate_effect_spread = build_replicate_effect_spread(
        effects, replicate_inventory, replicate_table
    )
    replicate_concordance = build_replicate_concordance(effects)
    replicate_blocked_analyses = build_replicate_blocked_analyses(available)

    # -- validate (structural / arithmetic only) ----------------------------
    validate_all(
        replicate_regime,
        replicate_effect_spread,
        replicate_concordance,
        replicate_blocked_analyses,
        report,
    )

    # -- report facts -------------------------------------------------------
    report.n_replicate_regime = int(len(replicate_regime))
    report.n_replicate_effect_spread = int(len(replicate_effect_spread))
    report.n_replicate_concordance = int(len(replicate_concordance))
    report.n_replicate_blocked_analyses = int(len(replicate_blocked_analyses))
    report.blocked_analyses = replicate_blocked_analyses.to_dict("records")
    report.facts = _facts(
        replicate_regime, replicate_concordance, replicate_blocked_analyses
    )
    report.add(
        "s6 replicate tables built",
        "global",
        True,
        f"{report.n_replicate_regime} regime, "
        f"{report.n_replicate_effect_spread} effect-spread, "
        f"{report.n_replicate_concordance} concordance, "
        f"{report.n_replicate_blocked_analyses} blocked-analysis rows",
    )

    tables = S6Tables(
        replicate_regime=replicate_regime,
        replicate_effect_spread=replicate_effect_spread,
        replicate_concordance=replicate_concordance,
        replicate_blocked_analyses=replicate_blocked_analyses,
    )
    return tables, report


def run_s6(
    replicate_inventory_path: str | Path,
    output_dir: str | Path,
    replicate_table_path: str | Path | None = None,
) -> tuple[S6Tables, S6Report]:
    """Full S6: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``replicate_regime.parquet``, ``replicate_effect_spread.parquet``,
    ``replicate_concordance.parquet``, ``replicate_blocked_analyses.parquet`` and
    ``replicate_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s6(replicate_inventory_path, replicate_table_path)
    write_tables(
        tables.replicate_regime,
        tables.replicate_effect_spread,
        tables.replicate_concordance,
        tables.replicate_blocked_analyses,
        output_dir,
    )
    write_replicate_summary(report, output_dir)
    return tables, report


def _facts(
    replicate_regime: pd.DataFrame,
    replicate_concordance: pd.DataFrame,
    replicate_blocked_analyses: pd.DataFrame,
) -> dict[str, object]:
    """Compact replicate-layer facts for the summary JSON."""
    facts: dict[str, object] = {}
    if not replicate_regime.empty:
        facts["n_replicates_by_serotype"] = {
            str(row["serotype"]): int(row["n_replicates"])
            for row in replicate_regime.to_dict("records")
        }
        facts["n_serotypes_residue_licensed"] = int(
            replicate_regime["residue_claims_licensed"].astype(bool).sum()
        )
        facts["n_serotypes_with_per_run_effects"] = int(
            replicate_regime["per_replicate_effects_available"].astype(bool).sum()
        )
    if not replicate_concordance.empty:
        class_counts = (
            replicate_concordance["concordance_class"].value_counts().to_dict()
        )
        facts["concordance_class_counts"] = {
            str(k): int(v) for k, v in class_counts.items()
        }
    if not replicate_blocked_analyses.empty:
        facts["n_blocked_analyses"] = int(
            (~replicate_blocked_analyses["available"].astype(bool)).sum()
        )
        facts["n_available_analyses"] = int(
            replicate_blocked_analyses["available"].astype(bool).sum()
        )
    return facts
