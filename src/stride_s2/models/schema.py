"""Frozen schema constants for Stage S2 outputs.

S2 is the **per-serotype reduction** layer. It consumes the S0 STRIDE table
(the Level-2 profile+mechanism payload across all seven hierarchy scales) and
the S1B ``residue_annotation`` / ``domain_annotation`` structural-label tables,
and derives per-serotype reproducibility summaries: the reproducibility
landscape, the achieved-resolution census, domain-scale reproducibility, the
signed/significant screen, and a per-serotype scorecard.

This module is the single source of truth for:

- the S0/S1B columns S2 depends on,
- the provisional gate threshold and the ρ* band S2 sweeps,
- the two-tier (licensed / exploratory) output vocabulary, and
- the five output-table schemas, their keys, and the artifact filenames.

Every design decision here follows the authoritative design document:

- **the gate is uncalibrated** — every resolution statement is reported as a
  function of ρ* over a band, never as a single pass/fail verdict (§0.1, §5.3);
- **domain scale is the licensed claim level at K=3**; residue-scale outputs are
  labelled exploratory / outside the operating range (§0.1, §5.3, §5.4);
- **serotype is the unit of biological replication** — S2 stays strictly
  *within* a serotype and performs no cross-serotype test (that is S5) (§5.1).

S2 never re-reads the raw STRIDE files, never recomputes ρ from its variance
components (the design forbids it — §1.2), and never modifies any earlier stage.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S0 STRIDE-table columns S2 depends on (a subset of the S0 frozen schema)
# ===========================================================================
#: Columns S2 reads from the S0 STRIDE table. The profile payload carries the
#: reproducibility numbers at every scale; the mechanism payload (gated rows
#: only) carries the signed effect and its CI.
STRIDE_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "scale_level",
    "scale_index",
    "region_id",
    "rho",
    "gated",
    "beta",
    "beta_se",
    "tau2",
    "sigma2_bar",
    "coherence",
    "h_chain",
    "h_domain",
    "is_gated_scale",
    "mech_direction",
    "mech_beta_signed",
    "mech_beta_ci_lower",
    "mech_beta_ci_upper",
    "mech_reproducible_magnitude_energy",
)

# ===========================================================================
# S1B annotation columns S2 depends on
# ===========================================================================
RESIDUE_ANNOTATION_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "domain_status",
    "conservation_class",
)
DOMAIN_ANNOTATION_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "domain_status",
    "n_residues",
)

# ===========================================================================
# Scale levels (mirrors the S0 grammar; re-declared so S2 is self-contained)
# ===========================================================================
RESIDUE_SCALE_LEVEL: Final[str] = "residue"
DOMAIN_SCALE_LEVEL: Final[str] = "domain"

#: Scale ordering from finest (0 = residue) to coarsest (6 = complex). Re-gating
#: picks the finest scale (smallest index) at which ρ ≥ ρ*.
SCALE_LEVEL_TO_INDEX: Final[dict[str, int]] = {
    "residue": 0,
    "secondary_structure": 1,
    "motif": 2,
    "domain": 3,
    "chain": 4,
    "protein": 5,
    "complex": 6,
}

# ===========================================================================
# The uncalibrated gate and the ρ* band
# ===========================================================================
#: The provisional configured threshold shipped in every mechanism file. The
#: design is emphatic that this is *provisional* and uncalibrated (§0.1); S2
#: never treats a single ρ* as a verdict.
PROVISIONAL_RHO_STAR: Final[float] = 0.5

#: Directional-coherence gate threshold (mechanism gate parameter). A region is
#: directionally *coherent* at/above this; below it the signed claim is "mixed".
COHERENCE_THRESHOLD: Final[float] = 0.6

#: The default ρ* band S2 sweeps: 0.50 → 0.90 inclusive in steps of 0.05, per
#: the design's stated example band (§5.3, §3.2). The band is a first-class
#: input — callers may pass their own — but this is the frozen default.
DEFAULT_RHO_STAR_BAND: Final[tuple[float, ...]] = (
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
    0.85,
    0.90,
)

#: Number of decimals ρ* values are rounded to (keeps the swept band and the
#: emitted rows deterministic and free of binary-float noise).
RHO_STAR_DECIMALS: Final[int] = 4

# ===========================================================================
# Two-tier output vocabulary (design §5.4)
# ===========================================================================
#: Tier A — licensed at K=3: domain-scale, uncertainty-aware.
TIER_LICENSED: Final[str] = "licensed"
#: Tier B — exploratory: residue-scale, uncalibrated, outside the operating range.
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

#: Sentinel used in the census for loci that gate at no scale under a given ρ*
#: (i.e. ρ < ρ* at every scale, so no resolution is achieved). Not observed in
#: the current files at ρ* ≤ 0.9 but handled explicitly so the census is total.
SCALE_UNRESOLVED: Final[str] = "unresolved"

# ===========================================================================
# Output table 1 — resolution_census.parquet
# ===========================================================================
#: One row per (serotype, rho_star, gated_scale_level): how many loci re-gate at
#: that scale under that ρ*. The achieved-resolution census swept over the band.
#: Unique key: (serotype, rho_star, gated_scale_level).
RESOLUTION_CENSUS_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "rho_star",
    "is_provisional_rho_star",  # bool: this ρ* == the provisional 0.5
    "gated_scale_level",        # residue | ... | complex | unresolved
    "gated_scale_index",        # 0..6, or -1 for the unresolved sentinel
    "tier",                     # licensed (domain+coarser) | exploratory (finer)
    "n_loci",                   # loci whose finest ρ≥ρ* scale is this level
)
RESOLUTION_CENSUS_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "rho_star",
    "gated_scale_level",
)

# ===========================================================================
# Output table 2 — residue_landscape.parquet  (Tier B — exploratory)
# ===========================================================================
#: One row per (serotype, canon_label): the residue-scale reproducibility
#: landscape. ρ, unsigned magnitude, coherence and direction at residue scale,
#: plus the finest scale the residue re-gates at under the provisional ρ*.
#: Unique key: (serotype, canon_label). Every row is Tier B (exploratory).
RESIDUE_LANDSCAPE_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "domain_status",
    "conservation_class",
    "rho_residue",                    # ρ at residue scale
    "beta_residue",                   # unsigned reproducible magnitude, residue
    "coherence_residue",              # directional coherence at residue scale
    "tau2_residue",
    "sigma2_bar_residue",
    "gated_scale_level_provisional",  # finest ρ≥0.5 scale for this locus
    "gated_scale_index_provisional",
    "gates_at_residue_provisional",   # bool: finest scale is residue at ρ*=0.5
    "tier",                           # always exploratory
)
RESIDUE_LANDSCAPE_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

# ===========================================================================
# Output table 3 — domain_reproducibility.parquet  (Tier A — licensed)
# ===========================================================================
#: One row per (serotype, chain, domain): domain-scale reproducibility. ρ,
#: unsigned magnitude, coherence and the variance components read from the
#: domain-scale profile row of the domain's region, plus a directional-coherence
#: label. Unique key: (serotype, chain, domain). Every row is Tier A (licensed).
DOMAIN_REPRODUCIBILITY_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "domain_status",
    "region_id",
    "n_residues",             # member residue count (from S1B domain annotation)
    "rho_domain",             # ρ at domain scale for the domain region
    "beta_domain",            # unsigned reproducible magnitude at domain scale
    "coherence_domain",       # directional coherence of the domain region
    "tau2_domain",
    "sigma2_bar_domain",
    "is_coherent",            # bool: coherence ≥ COHERENCE_THRESHOLD
    "tier",                   # always licensed
)
DOMAIN_REPRODUCIBILITY_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
)

# ===========================================================================
# Output table 4 — signed_screen.parquet
# ===========================================================================
#: One row per (serotype, canon_label, rho_star) over the gated mechanisms: the
#: signed/significant screen. A mechanism passes when direction ≠ mixed AND its
#: β CI excludes 0 AND ρ ≥ ρ*. Swept across the band so the screen is reported
#: as a function of ρ*, never at a single threshold.
#: Unique key: (serotype, canon_label, rho_star).
SIGNED_SCREEN_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "rho_star",
    "is_provisional_rho_star",
    "gated_scale_level",       # the mechanism's own gated scale
    "tier",                    # licensed if gated scale is domain+coarser
    "rho",                     # ρ of the gated region
    "direction",               # increase | decrease | mixed
    "is_signed",               # bool: direction != mixed
    "beta_signed",             # signed effect (null when mixed)
    "beta_ci_lower",
    "beta_ci_upper",
    "ci_excludes_zero",        # bool: signed AND CI strictly excludes 0
    "meets_rho_star",          # bool: rho >= rho_star
    "passes_screen",           # bool: signed AND ci_excludes_zero AND meets_rho_star
)
SIGNED_SCREEN_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "rho_star",
)

# ===========================================================================
# Output table 5 — serotype_summary.parquet
# ===========================================================================
#: One row per (serotype, rho_star): the per-serotype scorecard swept over the
#: band. Loci counts, resolution census rollups, direction counts, and the ρ
#: distribution (median/IQR over residue-scale ρ). Mirrors design T1 / §3.6.
#: Unique key: (serotype, rho_star).
SEROTYPE_SUMMARY_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "rho_star",
    "is_provisional_rho_star",
    "n_loci",                       # distinct residue loci in the serotype
    "n_gated_residue",              # loci re-gating at residue scale under ρ*
    "n_gated_domain_or_coarser",   # loci re-gating at domain scale or coarser
    "n_unresolved",                 # loci gating at no scale under ρ*
    "frac_gated_residue",           # n_gated_residue / n_loci
    "n_mechanisms",                 # gated mechanisms (one per gated region)
    "n_signed",                     # mechanisms with direction != mixed
    "n_mixed",                      # mechanisms with direction == mixed
    "n_signed_significant",         # signed mechanisms passing the screen at ρ*
    "frac_mixed",                   # n_mixed / n_mechanisms
    "rho_residue_median",           # median residue-scale ρ
    "rho_residue_q1",               # 25th percentile residue-scale ρ
    "rho_residue_q3",               # 75th percentile residue-scale ρ
    "rho_residue_min",
    "rho_residue_max",
)
SEROTYPE_SUMMARY_KEY: Final[tuple[str, ...]] = ("serotype", "rho_star")

# ===========================================================================
# input + output artifact filenames
# ===========================================================================
IN_STRIDE_TABLE: Final[str] = "stride_table.parquet"
IN_RESIDUE_ANNOTATION: Final[str] = "residue_annotation.parquet"
IN_DOMAIN_ANNOTATION: Final[str] = "domain_annotation.parquet"

OUT_RESOLUTION_CENSUS: Final[str] = "resolution_census.parquet"
OUT_RESIDUE_LANDSCAPE: Final[str] = "residue_landscape.parquet"
OUT_DOMAIN_REPRODUCIBILITY: Final[str] = "domain_reproducibility.parquet"
OUT_SIGNED_SCREEN: Final[str] = "signed_screen.parquet"
OUT_SEROTYPE_SUMMARY: Final[str] = "serotype_summary.parquet"
OUT_REDUCTION_SUMMARY: Final[str] = "reduction_summary.json"
