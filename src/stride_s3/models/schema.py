"""Frozen schema constants for Stage S3 outputs.

S3 is the **hierarchy reduction** layer. It consumes the S0 STRIDE table (the
profile: every locus evaluated at every one of the seven hierarchy scales) and
reduces it *along the scale axis* into per-locus and per-chain hierarchy
products:

- the ρ-vs-scale curve for every locus (how ρ climbs as the region coarsens);
- the domain-vs-residue reproducibility gap Δρ = ρ(domain) − ρ(residue), which
  flags effects that are reproducible only once aggregated;
- the monotonicity / upward-closure (I2) audit — whether ρ is non-decreasing as
  the scale coarsens, per locus; and
- the chain-level contrast (e.g. NS2B cofactor vs NS3 protease) of aggregate
  reproducibility and signed direction.

This module is the single source of truth for:

- the S0 STRIDE-table columns S3 depends on,
- the scale grammar and the two-tier (licensed / exploratory) vocabulary,
- the provisional gate constant, and
- the four output-table schemas, their keys, and the artifact filenames.

Every design decision here follows the authoritative design document:

- **the gate is uncalibrated** — S3 reports ρ as a continuous score across the
  scale axis and never emits a calibrated pass/fail verdict (§0.1, §5.3);
- **domain scale is the licensed claim level at K=3**; residue-scale products are
  labelled exploratory / outside the operating range (§0.1, §5.3, §5.4);
- ρ is **never recomputed** from its variance components — the ``rho`` column is
  read directly (§1.2).

S3 never re-reads the raw STRIDE files, never sweeps a ρ* band (that is S2/S5 —
S3's hierarchy products are ρ*-independent descriptions of the profile), and
never modifies any earlier stage.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S0 STRIDE-table columns S3 depends on (a subset of the S0 frozen schema)
# ===========================================================================
#: Columns S3 reads from the S0 STRIDE table. The profile payload carries ρ at
#: every scale; the parsed hierarchy columns carry chain/domain identity; the
#: mechanism payload (gated rows only) carries the signed direction used by the
#: chain contrast.
STRIDE_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "scale_level",
    "scale_index",
    "region_id",
    "rho",
    "gated",
    "beta",
    "coherence",
    "h_chain",
    "h_domain",
    "is_gated_scale",
    "mech_direction",
)

# ===========================================================================
# Scale grammar (mirrors the S0 grammar; re-declared so S3 is self-contained)
# ===========================================================================
RESIDUE_SCALE_LEVEL: Final[str] = "residue"
DOMAIN_SCALE_LEVEL: Final[str] = "domain"

#: Scale ordering from finest (0 = residue) to coarsest (6 = complex). ρ is
#: expected to be non-decreasing along this axis (the I2 upward-closure property).
SCALE_LEVEL_TO_INDEX: Final[dict[str, int]] = {
    "residue": 0,
    "secondary_structure": 1,
    "motif": 2,
    "domain": 3,
    "chain": 4,
    "protein": 5,
    "complex": 6,
}
#: The inverse map, index → level.
SCALE_INDEX_TO_LEVEL: Final[dict[int, str]] = {
    idx: level for level, idx in SCALE_LEVEL_TO_INDEX.items()
}
#: The full ordered scale-level tuple, finest → coarsest.
SCALE_LEVELS_ORDERED: Final[tuple[str, ...]] = tuple(
    SCALE_INDEX_TO_LEVEL[i] for i in range(len(SCALE_INDEX_TO_LEVEL))
)
#: Number of hierarchy scales (grammar depth): every locus has exactly this many
#: profile rows, one per scale.
N_SCALES: Final[int] = len(SCALE_LEVEL_TO_INDEX)

# ===========================================================================
# The uncalibrated gate (provisional; never treated as a verdict)
# ===========================================================================
#: The provisional configured threshold shipped in every mechanism file. The
#: design is emphatic that this is *provisional* and uncalibrated (§0.1). S3 uses
#: it only to label the profile's own gated scale and to define the
#: "distributed effect" flag; it makes no calibrated claim.
PROVISIONAL_RHO_STAR: Final[float] = 0.5

#: Decimals ρ-derived quantities are rounded to (keeps emitted rows deterministic
#: and free of binary-float noise).
RHO_DECIMALS: Final[int] = 6

# ===========================================================================
# Two-tier output vocabulary (design §5.4)
# ===========================================================================
#: Tier A — licensed at K=3: domain-scale and coarser.
TIER_LICENSED: Final[str] = "licensed"
#: Tier B — exploratory: finer than domain (residue / SS / motif), uncalibrated.
TIER_EXPLORATORY: Final[str] = "exploratory"
TIERS: Final[tuple[str, ...]] = (TIER_LICENSED, TIER_EXPLORATORY)

# ===========================================================================
# Signed-direction vocabulary (mirrors the mechanism ``direction`` field)
# ===========================================================================
DIRECTION_INCREASE: Final[str] = "increase"
DIRECTION_DECREASE: Final[str] = "decrease"
DIRECTION_MIXED: Final[str] = "mixed"
DIRECTIONS: Final[tuple[str, ...]] = (
    DIRECTION_INCREASE,
    DIRECTION_DECREASE,
    DIRECTION_MIXED,
)

# ===========================================================================
# Output table 1 — scale_curve.parquet
# ===========================================================================
#: One row per (serotype, canon_label, scale_index): the ρ-vs-scale trajectory of
#: a locus, with per-step and cumulative gains. This is the mechanical basis for
#: "how spread out is the reproducible signal" (§3.4, figure F7).
#: Unique key: (serotype, canon_label, scale_index).
SCALE_CURVE_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "scale_index",
    "scale_level",
    "rho",
    "rho_prev",             # ρ at the next-finer scale (NaN at residue)
    "rho_step_gain",        # rho - rho_prev (NaN at residue)
    "rho_cumulative_gain",  # rho - rho(residue)
    "is_gated_scale",       # profile's gated flag: finest ρ≥provisional scale
    "tier",                 # licensed (domain+coarser) | exploratory (finer)
)
SCALE_CURVE_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "scale_index",
)

# ===========================================================================
# Output table 2 — resolution_gap.parquet
# ===========================================================================
#: One row per (serotype, canon_label): the domain-vs-residue reproducibility gap.
#: ``is_distributed`` flags a locus that is not reproducible at residue scale but
#: is at domain scale under the provisional gate — a distributed effect visible
#: only when aggregated (§3.4).
#: Unique key: (serotype, canon_label).
RESOLUTION_GAP_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "rho_residue",                 # ρ at residue scale
    "rho_domain",                  # ρ at domain scale
    "delta_rho_domain_residue",    # rho_domain - rho_residue
    "rho_min",                     # min ρ across all scales
    "rho_max",                     # max ρ across all scales
    "gated_scale_level",           # the profile's gated scale for this locus
    "gated_scale_index",
    "rho_at_gated",                # ρ at the gated scale
    "is_distributed",              # bool: residue < ρ* ≤ domain (aggregated-only)
    "domain_tier",                 # always licensed (the aggregated claim level)
)
RESOLUTION_GAP_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

# ===========================================================================
# Output table 3 — monotonicity_audit.parquet
# ===========================================================================
#: One row per (serotype, canon_label): the upward-closure (I2) audit. ρ should
#: be non-decreasing as the scale coarsens; this table records whether it is and,
#: if not, where and by how much it drops (§3.4).
#: Unique key: (serotype, canon_label).
MONOTONICITY_AUDIT_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "n_scales",                    # always N_SCALES (7)
    "is_monotone",                 # bool: ρ non-decreasing finest → coarsest
    "n_violations",                # count of adjacent steps where ρ decreases
    "max_decrease",                # largest single-step ρ drop (0.0 if monotone)
    "first_violation_scale_index", # finer index of the first drop, or -1
    "rho_residue",                 # ρ at residue scale (for context)
    "rho_complex",                 # ρ at complex scale (for context)
)
MONOTONICITY_AUDIT_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

# ===========================================================================
# Output table 4 — chain_contrast.parquet
# ===========================================================================
#: One row per (serotype, chain): the chain-level contrast (e.g. NS2B vs NS3) of
#: aggregate reproducibility and signed direction (§3.4). Aggregates are taken
#: over the chain's member loci; no region-constant assumption is made.
#: Unique key: (serotype, chain). Chain scale is coarser than domain → Tier A.
CHAIN_CONTRAST_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "n_loci",                  # distinct residue loci in the chain
    "rho_residue_mean",        # aggregate residue-scale ρ (Tier B, exploratory)
    "rho_residue_median",
    "rho_residue_min",
    "rho_residue_max",
    "rho_chain_mean",          # aggregate chain-scale ρ (Tier A, licensed)
    "rho_chain_median",
    "beta_residue_mean",       # aggregate reproducible magnitude at residue scale
    "n_mechanisms",            # gated mechanisms in the chain
    "n_increase",              # gated mechanisms with direction == increase
    "n_decrease",              # gated mechanisms with direction == decrease
    "n_mixed",                 # gated mechanisms with direction == mixed
    "n_signed",                # n_increase + n_decrease
    "tier",                    # always licensed (chain scale is coarser than domain)
)
CHAIN_CONTRAST_KEY: Final[tuple[str, ...]] = ("serotype", "chain")

# ===========================================================================
# input + output artifact filenames
# ===========================================================================
IN_STRIDE_TABLE: Final[str] = "stride_table.parquet"

OUT_SCALE_CURVE: Final[str] = "scale_curve.parquet"
OUT_RESOLUTION_GAP: Final[str] = "resolution_gap.parquet"
OUT_MONOTONICITY_AUDIT: Final[str] = "monotonicity_audit.parquet"
OUT_CHAIN_CONTRAST: Final[str] = "chain_contrast.parquet"
OUT_HIERARCHY_SUMMARY: Final[str] = "hierarchy_summary.json"
