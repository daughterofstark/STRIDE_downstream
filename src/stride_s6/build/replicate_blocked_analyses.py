"""Build the blocked-analysis ledger.

A first-class, machine-readable record of the replicate-level analyses whose
availability depends on the per-run inputs, so the design's blocked subset (§3.1
☞, §4.1, §4.2) is *documented* rather than approximated. One row per analysis:

- **leave-one-replicate-out stability** — always blocked: recomputing ρ and the
  gated scale on replicate subsets requires re-running STRIDE, which this
  read-only pipeline never does (only the K-aggregate is available, §4.2);
- **per-run rank concordance** and **per-run effect spread** — available iff the
  S0 replicate table exposes per-run θ (i.e. the Level-1 correlation CSVs were
  supplied); blocked otherwise, with the required (absent) input named.

Deterministic (fixed row order); depends only on the availability flag.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    ANALYSIS_EFFECT_SPREAD,
    ANALYSIS_LORO_STABILITY,
    ANALYSIS_RANK_CONCORDANCE,
    REPLICATE_BLOCKED_ANALYSES_COLUMNS,
    REQUIRED_INPUT_PER_RUN_CSVS,
    REQUIRED_INPUT_STRIDE_RERUN,
    STATUS_AVAILABLE,
    STATUS_BLOCKED,
)

_PER_RUN_DESIGN_REF = "§3.1 (☞), §4.1"
_LORO_DESIGN_REF = "§3.1 (☞), §4.2"


def build_replicate_blocked_analyses(
    per_replicate_effects_available: bool,
) -> pd.DataFrame:
    """Return the blocked-analysis ledger.

    Parameters
    ----------
    per_replicate_effects_available
        Whether the S0 replicate table exposed per-run θ for at least the minimum
        number of runs. Governs the status of the two per-run analyses; the
        leave-one-replicate-out row is always blocked.
    """
    per_run_status = (
        STATUS_AVAILABLE
        if per_replicate_effects_available
        else STATUS_BLOCKED
    )
    if per_replicate_effects_available:
        concordance_reason = (
            "per-run θ present in the replicate table; computed over the "
            "complete-case positions of each serotype"
        )
        spread_reason = (
            "per-run θ present in the replicate table; summarised descriptively "
            "per observed position"
        )
    else:
        concordance_reason = (
            "no per-run θ in the pipeline (replicate table absent or empty); "
            "concordance is not reconstructible from the K-aggregate τ²/σ̄²"
        )
        spread_reason = (
            "no per-run θ in the pipeline (replicate table absent or empty); "
            "per-run effects cannot be recovered from the K-aggregate"
        )

    rows = [
        {
            "analysis_id": ANALYSIS_RANK_CONCORDANCE,
            "description": (
                "per-run θ rank concordance across replicates "
                "(Kendall's W / mean pairwise Spearman)"
            ),
            "status": per_run_status,
            "available": bool(per_replicate_effects_available),
            "reason": concordance_reason,
            "required_input": REQUIRED_INPUT_PER_RUN_CSVS,
            "design_ref": _PER_RUN_DESIGN_REF,
        },
        {
            "analysis_id": ANALYSIS_EFFECT_SPREAD,
            "description": (
                "descriptive across-run spread of each position's per-run θ"
            ),
            "status": per_run_status,
            "available": bool(per_replicate_effects_available),
            "reason": spread_reason,
            "required_input": REQUIRED_INPUT_PER_RUN_CSVS,
            "design_ref": _PER_RUN_DESIGN_REF,
        },
        {
            "analysis_id": ANALYSIS_LORO_STABILITY,
            "description": (
                "leave-one-replicate-out stability of ρ and the gated scale"
            ),
            "status": STATUS_BLOCKED,
            "available": False,
            "reason": (
                "recomputing ρ and re-gating on replicate subsets requires "
                "re-running STRIDE; only the K-aggregate is in the pipeline"
            ),
            "required_input": REQUIRED_INPUT_STRIDE_RERUN,
            "design_ref": _LORO_DESIGN_REF,
        },
    ]

    return pd.DataFrame(
        rows, columns=list(REPLICATE_BLOCKED_ANALYSES_COLUMNS)
    ).reset_index(drop=True)
