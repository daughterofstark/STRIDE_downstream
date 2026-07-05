r"""Frozen schema constants for Stage S6 outputs.

S6 is the **replicate layer** — the per-replicate-observation stage. It is the
only stage that consumes the replicate-specific inputs: the S1A
``replicate_inventory`` (the replicate-structure index — how many runs contain
each ``canon_label`` per serotype, and whether the replicate set is complete) and,
when it exists, the S0 ``replicate_table`` (the Level-1 per-run observations — one
row per ``(serotype, replicate, canon_label)`` carrying that run's effect ``r``,
i.e. the per-run θ). From these it builds the replicate-axis products the design
calls for in §3.1:

- **replicate-regime accounting** — per serotype, the replicate count *K*, the
  completeness of the replicate set, whether residue-scale claims are licensed
  (K ≥ 5), and whether genuine per-run effects are available at all;
- **per-run effect spread** — a *descriptive* across-run summary of each shared
  position's per-run θ (mean / sd / range), available only when the per-run
  observations exist;
- **per-replicate rank concordance** — Kendall's *W* and the mean pairwise
  Spearman correlation of the per-run effect rankings across the runs of a
  serotype (design §3.1: "per-run θ correlation / Kendall's W across the runs"),
  available only when the per-run observations exist; and
- a **blocked-analysis ledger** — a first-class, machine-readable record of the
  replicate-level analyses the design flags as unavailable, with the reason and
  the required (absent) input for each, so the blocked subset is documented rather
  than approximated.

What S6 deliberately does **not** do (see :mod:`stride_s6` and ``docs/s6.md``):

- It does **not** re-implement the τ²-based replicate-disagreement mapping or the
  τ²/σ̄² ratio ("replicate-dominated vs sampling-dominated") diagnostic of design
  §3.1 — those read the *aggregate* variance components and are already produced by
  **S4** (``residue_variance`` / ``variance_budget``). S6 never reads
  ``stride_table``; it owns the *per-replicate* axis, not the aggregate one.
- It does **not** attempt **leave-one-replicate-out** stability of ρ or the gated
  scale — that requires re-running STRIDE on replicate subsets to recompute ρ and
  re-gate (only the K-aggregate is in the pipeline), so it is permanently blocked
  and is recorded as such (design §3.1 ☞, §4.2).
- It does **not** reconstruct per-run θ when the per-run correlation CSVs are
  absent (the real dengue upload's state, design §4.1): with no
  ``replicate_table``, the per-run spread and concordance are empty and are
  recorded as blocked, never approximated from the K-aggregate.

This module is the single source of truth for the input columns S6 depends on,
the availability / licensing thresholds, the concordance vocabulary, the blocked
analysis identifiers, and the output-table schemas / keys / filenames.

Design invariants honoured (authoritative design document):

- **Replicates are the sampling unit** for the per-run observations; S6 aggregates
  each position to per-run values and never pools frames across runs as if
  independent (§5.2).
- **Residue-scale outputs are exploratory** at K = 3 (outside the licensed
  operating range); the per-run spread and concordance are residue-scale and
  labelled exploratory (§0.1, §5.3).
- **The gate is uncalibrated**; S6 makes no calibrated claim and stamps
  ``calibrated = false`` in provenance (§0.1, §5.4).
- Per-run θ is **read, never recomputed**, and is **not reconstructible** from the
  K-aggregate (§4.2).
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S1A replicate-inventory columns S6 depends on (the replicate-structure index)
# ===========================================================================
#: Columns S6 reads from the S1A ``replicate_inventory`` — the authoritative
#: per-``(serotype, canon_label)`` replicate-structure map. This input is always
#: present (S1A always emits it, with ``n_replicates == 0`` everywhere when no
#: per-run correlation CSVs were supplied to S0).
REPLICATE_INVENTORY_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "n_replicates",
    "in_all_replicates",
)

# ===========================================================================
# S0 replicate-table columns S6 depends on (the Level-1 per-run observations)
# ===========================================================================
#: Columns S6 reads from the S0 ``replicate_table`` when it exists. This input is
#: **optional**: S0 only writes it when the per-run correlation CSVs (Level-1
#: inputs) were supplied. ``replicate_index`` identifies the run; ``r`` is that
#: run's effect (the per-run θ).
REPLICATE_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "replicate_index",
    "canon_label",
    "r",
)

#: Optional column carried through to the spread table's context when present.
REPLICATE_TABLE_DOMAIN_COLUMN: Final[str] = "domain_label"

# ===========================================================================
# Thresholds and rounding
# ===========================================================================
#: Number of decimals every emitted floating-point statistic is rounded to, so
#: outputs are byte-for-byte deterministic across runs and platforms.
ROUND_DECIMALS: Final[int] = 6

#: STRIDE's operating-range guidance licenses residue-scale claims only at
#: K ≥ 5 replicates; at K = 3 residue-scale outputs are exploratory (§0.1).
MIN_K_FOR_RESIDUE_LICENSE: Final[int] = 5

#: A serotype has usable per-run effects (and so can support a concordance /
#: spread computation) only when at least this many runs carry a finite θ.
MIN_REPLICATES_FOR_CONCORDANCE: Final[int] = 2

#: Rank concordance (Kendall's W) is only meaningful over at least this many
#: complete-case positions; below it the coefficient is reported as undefined.
MIN_POSITIONS_FOR_CONCORDANCE: Final[int] = 3

# ===========================================================================
# Two-tier vocabulary (design §5.4)
# ===========================================================================
#: Residue-scale, single-serotype, caveated — the per-run spread and concordance.
TIER_EXPLORATORY: Final[str] = "exploratory"
#: Domain-scale / licensed — unused by S6's residue-scale claim tables, kept for
#: vocabulary parity with the rest of the pipeline.
TIER_LICENSED: Final[str] = "licensed"
TIERS: Final[tuple[str, ...]] = (TIER_LICENSED, TIER_EXPLORATORY)

# ===========================================================================
# Concordance-class vocabulary (per-run rank concordance, design §3.1)
# ===========================================================================
#: Kendall's W ≥ ``CONCORDANCE_STRONG_W`` → strong; ≥ ``CONCORDANCE_MODERATE_W``
#: → moderate; below → weak; too few complete positions / runs → insufficient.
CONCORDANCE_STRONG_W: Final[float] = 0.7
CONCORDANCE_MODERATE_W: Final[float] = 0.4

CONCORDANCE_STRONG: Final[str] = "strong"
CONCORDANCE_MODERATE: Final[str] = "moderate"
CONCORDANCE_WEAK: Final[str] = "weak"
CONCORDANCE_INSUFFICIENT: Final[str] = "insufficient"
CONCORDANCE_CLASSES: Final[tuple[str, ...]] = (
    CONCORDANCE_STRONG,
    CONCORDANCE_MODERATE,
    CONCORDANCE_WEAK,
    CONCORDANCE_INSUFFICIENT,
)

# ===========================================================================
# Blocked-analysis ledger (design §3.1 ☞, §4.1, §4.2)
# ===========================================================================
#: Analysis identifiers recorded in ``replicate_blocked_analyses``.
ANALYSIS_RANK_CONCORDANCE: Final[str] = "per_replicate_rank_concordance"
ANALYSIS_EFFECT_SPREAD: Final[str] = "per_replicate_effect_spread"
ANALYSIS_LORO_STABILITY: Final[str] = "leave_one_replicate_out_stability"

#: Availability status vocabulary for the ledger.
STATUS_AVAILABLE: Final[str] = "available"
STATUS_BLOCKED: Final[str] = "blocked"
LEDGER_STATUSES: Final[tuple[str, ...]] = (STATUS_AVAILABLE, STATUS_BLOCKED)

#: The per-run correlation CSVs (Level-1 inputs) that the per-run analyses need.
REQUIRED_INPUT_PER_RUN_CSVS: Final[str] = (
    "per-run correlation tables ({SEROTYPE}_correlations_v5.csv, one per replicate)"
)
#: What leave-one-replicate-out additionally requires.
REQUIRED_INPUT_STRIDE_RERUN: Final[str] = (
    "re-running STRIDE on replicate subsets to recompute ρ and re-gate"
)

# ===========================================================================
# Output-table schemas, keys, and filenames
# ===========================================================================
#: ``replicate_regime`` — one row per serotype: the replicate-structure ledger.
REPLICATE_REGIME_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "n_replicates",
    "n_positions",
    "n_positions_in_all_replicates",
    "frac_complete",
    "residue_claims_licensed",
    "per_replicate_effects_available",
    "n_replicates_with_effects",
)
REPLICATE_REGIME_KEY: Final[tuple[str, ...]] = ("serotype",)

#: ``replicate_effect_spread`` — one row per observed ``(serotype, canon_label)``:
#: the descriptive across-run spread of the per-run θ (Tier B, exploratory).
REPLICATE_EFFECT_SPREAD_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain_label",
    "n_obs",
    "theta_mean",
    "theta_sd",
    "theta_min",
    "theta_max",
    "theta_range",
    "abs_theta_mean",
    "max_pairwise_abs_diff",
    "in_all_replicates",
    "tier",
)
REPLICATE_EFFECT_SPREAD_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

#: ``replicate_concordance`` — one row per serotype with usable per-run effects:
#: the rank concordance of the per-run effect rankings (Tier B, exploratory).
REPLICATE_CONCORDANCE_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "n_replicates_with_effects",
    "n_positions_complete",
    "kendalls_w",
    "mean_pairwise_spearman",
    "concordance_class",
    "tier",
)
REPLICATE_CONCORDANCE_KEY: Final[tuple[str, ...]] = ("serotype",)

#: ``replicate_blocked_analyses`` — one row per replicate-level analysis whose
#: availability depends on the (possibly absent) per-run inputs.
REPLICATE_BLOCKED_ANALYSES_COLUMNS: Final[tuple[str, ...]] = (
    "analysis_id",
    "description",
    "status",
    "available",
    "reason",
    "required_input",
    "design_ref",
)
REPLICATE_BLOCKED_ANALYSES_KEY: Final[tuple[str, ...]] = ("analysis_id",)

# ---- input filenames (defaults) ------------------------------------------
IN_REPLICATE_TABLE: Final[str] = "replicate_table.parquet"
IN_REPLICATE_INVENTORY: Final[str] = "replicate_inventory.parquet"

# ---- output filenames -----------------------------------------------------
OUT_REPLICATE_REGIME: Final[str] = "replicate_regime.parquet"
OUT_REPLICATE_EFFECT_SPREAD: Final[str] = "replicate_effect_spread.parquet"
OUT_REPLICATE_CONCORDANCE: Final[str] = "replicate_concordance.parquet"
OUT_REPLICATE_BLOCKED_ANALYSES: Final[str] = "replicate_blocked_analyses.parquet"
OUT_REPLICATE_SUMMARY: Final[str] = "replicate_summary.json"
