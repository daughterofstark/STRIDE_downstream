r"""Frozen schema constants for Stage S4 outputs.

S4 is the **uncertainty layer**. It consumes the S0 STRIDE table (the tidy
master profile: every locus at every one of the seven hierarchy scales, carrying
the variance components τ² and σ̄² and the gated mechanism payload) and derives
per-domain and per-residue uncertainty products:

- the variance-component budget per domain — how each region's non-reproducibility
  splits between τ² (replicate disagreement) and σ̄² (sampling noise);
- the τ²/σ̄² ratio diagnostic (replicate-dominated vs sampling-dominated regimes);
- the per-residue replicate-disagreement map, ranking positions by τ²;
- the CI-based significance screen over the gated mechanisms, with
  Benjamini–Hochberg FDR control across the positions of a serotype; and
- the β_se-weighted effect summary per domain (weighting by 1/β_se²).

This module is the single source of truth for:

- the S0 STRIDE-table columns S4 depends on,
- the scale grammar and the two-tier (licensed / exploratory) vocabulary,
- the provisional gate constants (ρ\* and α) and the variance-regime thresholds,
- and the four output-table schemas, their keys, and the artifact filenames.

Every design decision here follows the authoritative design document:

- **the gate is uncalibrated** — S4 reports uncertainty descriptively and makes
  no calibrated pass/fail resolution claim (§0.1, §5.3);
- **domain scale is the licensed claim level at K=3**; residue-scale products are
  labelled exploratory / outside the operating range (§0.1, §5.3, §5.4);
- ρ and the variance components are **read**, never recomputed (§1.2);
- **serotype is the unit of biological replication** — S4 stays strictly within a
  serotype and performs no cross-serotype test (that is S5); the FDR family is
  the positions of a single serotype (§5.2).

S4 never re-reads the raw STRIDE files, never sweeps a ρ\* band (its uncertainty
products are ρ\*-independent descriptions of the profile and the emitted
mechanisms), and never modifies any earlier stage.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S0 STRIDE-table columns S4 depends on (a subset of the S0 frozen schema)
# ===========================================================================
#: Columns S4 reads from the S0 STRIDE table. The profile payload carries the
#: variance components at every scale; the mechanism payload (gated rows only)
#: carries the signed effect, its SE and its CI, used by the significance screen.
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
    "mech_beta_se",
)

# ===========================================================================
# Scale grammar (mirrors the S0 grammar; re-declared so S4 is self-contained)
# ===========================================================================
RESIDUE_SCALE_LEVEL: Final[str] = "residue"
DOMAIN_SCALE_LEVEL: Final[str] = "domain"

#: Scale ordering from finest (0 = residue) to coarsest (6 = complex).
SCALE_LEVEL_TO_INDEX: Final[dict[str, int]] = {
    "residue": 0,
    "secondary_structure": 1,
    "motif": 2,
    "domain": 3,
    "chain": 4,
    "protein": 5,
    "complex": 6,
}
#: The domain scale index — the boundary between exploratory (finer) and licensed
#: (this and coarser) tiers.
DOMAIN_SCALE_INDEX: Final[int] = SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]

# ===========================================================================
# The uncalibrated gate (provisional; never treated as a verdict)
# ===========================================================================
#: The provisional configured threshold shipped in every mechanism file. The
#: design is emphatic that this is *provisional* and uncalibrated (§0.1). S4 uses
#: it only to label the profile's gated scale; it makes no calibrated claim.
PROVISIONAL_RHO_STAR: Final[float] = 0.5

#: The gate's significance level (mechanism gate parameter ``alpha``). Used as the
#: Benjamini–Hochberg FDR level and as the two-sided level whose critical value
#: defines "CI excludes 0" for a normal (Wald) interval.
GATE_ALPHA: Final[float] = 0.05

#: Number of decimals ρ/variance-derived quantities are rounded to (keeps emitted
#: rows deterministic and free of binary-float noise).
RHO_DECIMALS: Final[int] = 6

#: Number of decimals p-values are rounded to.
P_DECIMALS: Final[int] = 8

# ===========================================================================
# Variance-regime vocabulary (τ²/σ̄² ratio diagnostic, design §3.1)
# ===========================================================================
#: A region is *replicate-dominated* when τ² carries at least this fraction of
#: the unreproduced variance (high run-to-run disagreement; candidate
#: metastability), *sampling-dominated* when it carries at most
#: ``SAMPLING_DOMINATED_FRACTION``, and *balanced* in between.
REPLICATE_DOMINATED_FRACTION: Final[float] = 0.6
SAMPLING_DOMINATED_FRACTION: Final[float] = 0.4

VARIANCE_REGIME_REPLICATE: Final[str] = "replicate_dominated"
VARIANCE_REGIME_SAMPLING: Final[str] = "sampling_dominated"
VARIANCE_REGIME_BALANCED: Final[str] = "balanced"
#: Used when the unreproduced variance (τ² + σ̄²) is zero, so no split is defined.
VARIANCE_REGIME_UNDEFINED: Final[str] = "undefined"
VARIANCE_REGIMES: Final[tuple[str, ...]] = (
    VARIANCE_REGIME_REPLICATE,
    VARIANCE_REGIME_SAMPLING,
    VARIANCE_REGIME_BALANCED,
    VARIANCE_REGIME_UNDEFINED,
)

# ===========================================================================
# Two-tier output vocabulary (design §5.4)
# ===========================================================================
#: Tier A — licensed at K=3: domain-scale and coarser.
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

# ===========================================================================
# Output table 1 — variance_budget.parquet  (Tier A — licensed, per domain)
# ===========================================================================
#: One row per (serotype, chain, domain): the domain-scale variance-component
#: budget. τ² and σ̄² are read from the domain region's domain-scale profile row
#: (region-level constants), and split into fractions with a regime label
#: (§3.1, §3.5, table T5). Unique key: (serotype, chain, domain).
VARIANCE_BUDGET_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "region_id",
    "rho_domain",
    "beta_domain",
    "beta_se_domain",
    "tau2",
    "sigma2_bar",
    "total_unreproduced",   # tau2 + sigma2_bar
    "frac_tau2",            # tau2 / total_unreproduced (NaN if total == 0)
    "frac_sigma2",          # sigma2_bar / total_unreproduced (NaN if total == 0)
    "tau2_sigma2_ratio",    # tau2 / sigma2_bar (NaN if sigma2_bar == 0)
    "variance_regime",      # replicate_dominated | sampling_dominated | balanced | undefined
    "tier",                 # always licensed
)
VARIANCE_BUDGET_KEY: Final[tuple[str, ...]] = ("serotype", "chain", "domain")

# ===========================================================================
# Output table 2 — residue_variance.parquet  (Tier B — exploratory, per residue)
# ===========================================================================
#: One row per (serotype, canon_label): the residue-scale variance decomposition
#: and the replicate-disagreement ranking. ``tau2_rank`` orders positions within
#: a serotype by τ² descending (rank 1 = the run-to-run disagreement is largest),
#: per the design's "rank positions by tau2" diagnostic (§3.1).
#: Unique key: (serotype, canon_label). Every row is Tier B (exploratory).
RESIDUE_VARIANCE_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "rho_residue",
    "beta_residue",
    "beta_se_residue",
    "tau2",
    "sigma2_bar",
    "total_unreproduced",
    "frac_tau2",
    "frac_sigma2",
    "tau2_sigma2_ratio",
    "variance_regime",
    "tau2_rank",            # 1-based rank within serotype by τ² descending
    "tier",                 # always exploratory
)
RESIDUE_VARIANCE_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

# ===========================================================================
# Output table 3 — significance_screen.parquet  (per gated mechanism, FDR)
# ===========================================================================
#: One row per gated mechanism (keyed by the residue locus that carries the gated
#: row): the CI-based significance screen with Benjamini–Hochberg FDR control.
#: ``significant_raw`` is the design's "CI excludes 0" test; ``significant_fdr``
#: is the BH-adjusted decision at level α across the *signed* mechanisms of the
#: serotype (the FDR family; §3.5, §5.2). Unique key: (serotype, canon_label).
SIGNIFICANCE_SCREEN_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "gated_scale_level",
    "tier",                 # licensed if gated scale is domain+coarser else exploratory
    "direction",
    "is_signed",            # direction != mixed
    "beta_signed",
    "beta_se",
    "beta_ci_lower",
    "beta_ci_upper",
    "ci_excludes_zero",     # signed AND CI strictly excludes 0
    "z_score",              # beta_signed / beta_se (NaN when unavailable)
    "p_value",              # two-sided Wald p-value (NaN when unavailable)
    "p_value_bh",           # Benjamini–Hochberg adjusted p (within serotype)
    "significant_raw",      # == ci_excludes_zero
    "significant_fdr",      # signed AND p_value_bh <= alpha
)
SIGNIFICANCE_SCREEN_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")

# ===========================================================================
# Output table 4 — domain_effect_summary.parquet  (Tier A — licensed, per domain)
# ===========================================================================
#: One row per (serotype, chain, domain) carrying at least one gated mechanism:
#: the β_se-weighted effect summary and the per-domain CI-exclusion fraction.
#: The weighted mean weights each signed effect by 1/β_se² (§3.5); the fraction
#: is the design's "fraction of signed mechanisms whose beta CI excludes 0, per
#: domain" (§3.5). Unique key: (serotype, chain, domain). Every row is Tier A.
DOMAIN_EFFECT_SUMMARY_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "n_mechanisms",         # gated mechanisms whose gated region sits in this domain
    "n_signed",             # of those, direction != mixed
    "n_mixed",
    "n_ci_excludes_zero",   # signed mechanisms whose CI excludes 0
    "frac_ci_excludes_zero",  # n_ci_excludes_zero / n_signed (NaN if n_signed == 0)
    "n_significant_fdr",    # signed mechanisms significant after BH-FDR
    "beta_weighted_mean",   # 1/β_se²-weighted mean of signed effects (NaN if none)
    "beta_weighted_se",     # sqrt(1 / Σ(1/β_se²)) (NaN if none)
    "beta_unweighted_mean",  # plain mean of signed effects (NaN if none)
    "tier",                 # always licensed
)
DOMAIN_EFFECT_SUMMARY_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
)

# ===========================================================================
# input + output artifact filenames
# ===========================================================================
IN_STRIDE_TABLE: Final[str] = "stride_table.parquet"

OUT_VARIANCE_BUDGET: Final[str] = "variance_budget.parquet"
OUT_RESIDUE_VARIANCE: Final[str] = "residue_variance.parquet"
OUT_SIGNIFICANCE_SCREEN: Final[str] = "significance_screen.parquet"
OUT_DOMAIN_EFFECT_SUMMARY: Final[str] = "domain_effect_summary.parquet"
OUT_UNCERTAINTY_SUMMARY: Final[str] = "uncertainty_summary.json"
