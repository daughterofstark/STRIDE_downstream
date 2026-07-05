r"""Frozen schema constants for Stage S7 outputs.

S7 is the **reporting layer** — the terminal stage. It consumes the already-reduced
outputs of S2–S6 and assembles the design's publication figures (F1–F8, design
table §3.7) and manuscript tables (T1–T5, design §3.8), deterministically and
provenance-stamped. It performs **no new statistics and no biological inference**:
every number it emits is read from a prior stage's table and merely selected,
ordered, and formatted. It never reads the raw STRIDE ``profile.csv`` /
``mechanism.json`` and never recomputes a quantity an earlier stage already
produced.

This module is the single source of truth for:

- the figure and table identifiers and their human-readable slugs,
- the S2–S6 input tables S7 reads and the columns it depends on from each,
- the output-artifact filenames (SVG figures; CSV + Parquet prepared data;
  Markdown + CSV + Parquet manuscript tables; the summary JSON),
- the deterministic ordering keys, rounding, and the minimal SVG style constants,
- the catalytic-machinery labels, and
- the mapping of each figure/table to its supporting source table(s).

Design decisions honoured (authoritative design document):

- **Reporting only.** S7 draws every value from S2–S6; it introduces no analysis to
  produce a figure (§5.1 stage table: "all above → F1–F8, T1–T5. Deterministic,
  provenance-stamped").
- **Deterministic + provenance-stamped.** Outputs are byte-reproducible (sorted
  rows, fixed rounding, hand-rendered SVG with no embedded timestamps); each run
  writes a summary carrying the input digests, the provisional ρ\*, ``calibrated =
  false`` and the K = 3 note (§5.4).
- **Two tiers / uncalibrated gate.** Figures and tables carry through the
  ``licensed`` / ``exploratory`` tier labels their source rows already have and make
  no calibrated pass/fail claim (§0.1, §5.3).
- **Vector, dependency-free figures.** Figures are emitted as deterministic SVG
  built from primitives — no plotting library is added, preserving the framework's
  determinism guarantee and no-unnecessary-dependencies convention. PNG/PDF
  rasterisation is intentionally not produced (it would require a heavyweight,
  version-sensitive renderer whose bytes are not reproducible); SVG is the
  publication-quality vector deliverable. See ``docs/s7.md``.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# Rounding / formatting
# ===========================================================================
#: Decimals every emitted float is rounded to, so outputs are byte-reproducible.
ROUND_DECIMALS: Final[int] = 6

#: Provisional gate threshold (uncalibrated); the value the prior stages swept
#: their bands around and the one S7's single-threshold views select.
PROVISIONAL_RHO_STAR: Final[float] = 0.5

#: Number of replicates in the design dataset (the K = 3 provenance note).
N_REPLICATES: Final[int] = 3

# ===========================================================================
# Catalytic machinery (design §3.3, F7 / T3)
# ===========================================================================
CATALYTIC_TRIAD_DOMAIN: Final[str] = "Catalytic Triad"
OXYANION_LOOP_DOMAIN: Final[str] = "Oxyanion Loop"
CATALYTIC_DOMAINS: Final[tuple[str, ...]] = (
    CATALYTIC_TRIAD_DOMAIN,
    OXYANION_LOOP_DOMAIN,
)

# ===========================================================================
# Two-tier vocabulary (carried through from source rows)
# ===========================================================================
TIER_LICENSED: Final[str] = "licensed"
TIER_EXPLORATORY: Final[str] = "exploratory"

# ===========================================================================
# Input tables (produced by S2–S6) and the columns S7 depends on
# ===========================================================================
#: filename → (producing stage, required columns). S7 asserts these columns are
#: present on load; it reads no other prior-stage tables.
IN_RESIDUE_LANDSCAPE: Final[str] = "residue_landscape.parquet"
IN_RESOLUTION_CENSUS: Final[str] = "resolution_census.parquet"
IN_SEROTYPE_SUMMARY: Final[str] = "serotype_summary.parquet"
IN_DOMAIN_REPRODUCIBILITY: Final[str] = "domain_reproducibility.parquet"
IN_SCALE_CURVE: Final[str] = "scale_curve.parquet"
IN_SIGNIFICANCE_SCREEN: Final[str] = "significance_screen.parquet"
IN_VARIANCE_BUDGET: Final[str] = "variance_budget.parquet"
IN_DOMAIN_EFFECT_SUMMARY: Final[str] = "domain_effect_summary.parquet"
IN_DOMAIN_SEROTYPE_MATRIX: Final[str] = "domain_serotype_matrix.parquet"
IN_POSITION_CONSERVATION: Final[str] = "position_conservation.parquet"
IN_REPLICATE_REGIME: Final[str] = "replicate_regime.parquet"
IN_REPLICATE_BLOCKED_ANALYSES: Final[str] = "replicate_blocked_analyses.parquet"

RESIDUE_LANDSCAPE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "rho_residue",
    "coherence_residue",
    "tier",
)
RESOLUTION_CENSUS_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "rho_star",
    "is_provisional_rho_star",
    "gated_scale_level",
    "gated_scale_index",
    "n_loci",
)
SEROTYPE_SUMMARY_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "rho_star",
    "is_provisional_rho_star",
    "n_loci",
    "n_gated_residue",
    "n_mechanisms",
    "n_signed",
    "n_mixed",
    "n_signed_significant",
    "frac_mixed",
    "rho_residue_median",
    "rho_residue_q1",
    "rho_residue_q3",
    "rho_residue_min",
    "rho_residue_max",
)
DOMAIN_REPRODUCIBILITY_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "rho_domain",
    "beta_domain",
    "coherence_domain",
    "is_coherent",
    "tier",
)
SCALE_CURVE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "scale_index",
    "scale_level",
    "rho",
    "tier",
)
SIGNIFICANCE_SCREEN_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "direction",
    "is_signed",
    "beta_signed",
    "beta_ci_lower",
    "beta_ci_upper",
    "significant_fdr",
    "tier",
)
VARIANCE_BUDGET_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "tau2",
    "sigma2_bar",
    "frac_tau2",
    "frac_sigma2",
    "tau2_sigma2_ratio",
    "variance_regime",
    "tier",
)
DOMAIN_EFFECT_SUMMARY_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "n_mechanisms",
    "n_signed",
    "n_ci_excludes_zero",
    "frac_ci_excludes_zero",
    "n_significant_fdr",
    "beta_weighted_mean",
    "beta_weighted_se",
    "tier",
)
DOMAIN_SEROTYPE_MATRIX_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "rho_domain",
    "beta_domain",
    "is_catalytic_domain",
    "tier",
)
POSITION_CONSERVATION_REQUIRED: Final[tuple[str, ...]] = (
    "canon_label",
    "chain",
    "domain",
    "n_serotypes_present",
    "n_serotypes_reproducible",
    "n_serotypes_signed_reproducible",
    "frac_reproducible",
    "conservation_class",
    "is_serotype_divergent",
    "is_catalytic_triad",
    "rho_residue_median",
    "tier",
)
REPLICATE_REGIME_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "n_replicates",
    "per_replicate_effects_available",
    "residue_claims_licensed",
)
REPLICATE_BLOCKED_ANALYSES_REQUIRED: Final[tuple[str, ...]] = (
    "analysis_id",
    "status",
    "available",
    "reason",
)

# ===========================================================================
# Figure identifiers, slugs, and source mapping (design §3.7 / table F1–F8)
# ===========================================================================
FIGURE_IDS: Final[tuple[str, ...]] = (
    "F1",
    "F2",
    "F3",
    "F4",
    "F5",
    "F6",
    "F7",
    "F8",
)
#: id → filename slug (used for ``<slug>.svg`` / ``.csv`` / ``.parquet``)
FIGURE_SLUGS: Final[dict[str, str]] = {
    "F1": "F1_reproducibility_landscape",
    "F2": "F2_resolution_census",
    "F3": "F3_domain_serotype_rho_heatmap",
    "F4": "F4_signed_effect_forest",
    "F5": "F5_cross_serotype_conservation",
    "F6": "F6_variance_composition",
    "F7": "F7_rho_vs_scale_catalytic",
    "F8": "F8_coherence_vs_rho",
}
#: id → human title (rendered into the SVG and recorded in the summary)
FIGURE_TITLES: Final[dict[str, str]] = {
    "F1": "Reproducibility landscape (rho vs residue, per serotype)",
    "F2": "Achieved-resolution census (gated scale per serotype)",
    "F3": "Domain x serotype rho heatmap",
    "F4": "Signed-effect forest (beta_signed +/- CI, coherent mechanisms)",
    "F5": "Cross-serotype conservation over shared positions",
    "F6": "Variance composition (tau^2 vs sigma-bar^2) by domain",
    "F7": "rho-vs-scale trajectories for catalytic regions",
    "F8": "Coherence vs rho (directional cleanliness of domains)",
}
#: id → the prior-stage source table(s) it is assembled from
FIGURE_SOURCES: Final[dict[str, tuple[str, ...]]] = {
    "F1": (IN_RESIDUE_LANDSCAPE,),
    "F2": (IN_RESOLUTION_CENSUS,),
    "F3": (IN_DOMAIN_SEROTYPE_MATRIX,),
    "F4": (IN_SIGNIFICANCE_SCREEN,),
    "F5": (IN_POSITION_CONSERVATION,),
    "F6": (IN_VARIANCE_BUDGET,),
    "F7": (IN_SCALE_CURVE,),
    "F8": (IN_DOMAIN_REPRODUCIBILITY,),
}

# ===========================================================================
# Table identifiers, slugs, and source mapping (design §3.8 / table T1–T5)
# ===========================================================================
TABLE_IDS: Final[tuple[str, ...]] = ("T1", "T2", "T3", "T4", "T5")
TABLE_SLUGS: Final[dict[str, str]] = {
    "T1": "T1_per_serotype_summary",
    "T2": "T2_domain_rho_signed_effect",
    "T3": "T3_catalytic_cross_serotype",
    "T4": "T4_top_shared_signed_positions",
    "T5": "T5_variance_component_budget",
}
TABLE_TITLES: Final[dict[str, str]] = {
    "T1": "Per-serotype summary",
    "T2": "Domain-level rho and signed effect, per serotype",
    "T3": "Catalytic Triad / Oxyanion Loop cross-serotype behaviour",
    "T4": "Top reproducible signed positions shared across serotypes",
    "T5": "Variance-component budget per domain",
}
TABLE_SOURCES: Final[dict[str, tuple[str, ...]]] = {
    "T1": (IN_SEROTYPE_SUMMARY,),
    "T2": (IN_DOMAIN_REPRODUCIBILITY, IN_DOMAIN_EFFECT_SUMMARY),
    "T3": (IN_DOMAIN_SEROTYPE_MATRIX,),
    "T4": (IN_POSITION_CONSERVATION,),
    "T5": (IN_VARIANCE_BUDGET,),
}

#: How many rows T4 (top shared signed positions) retains after ranking.
T4_TOP_N: Final[int] = 25

# ===========================================================================
# Output-table (manuscript-table) column schemas
# ===========================================================================
T1_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "n_loci",
    "n_gated_residue",
    "n_mechanisms",
    "n_signed",
    "n_mixed",
    "n_signed_significant",
    "frac_mixed",
    "rho_residue_median",
    "rho_residue_q1",
    "rho_residue_q3",
    "rho_residue_min",
    "rho_residue_max",
    "rho_star",
)
T2_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "rho_domain",
    "beta_domain",
    "coherence_domain",
    "is_coherent",
    "beta_weighted_mean",
    "beta_weighted_se",
    "n_signed",
    "n_significant_fdr",
    "tier",
)
T3_COLUMNS: Final[tuple[str, ...]] = (
    "domain",
    "serotype",
    "chain",
    "rho_domain",
    "beta_domain",
    "is_catalytic_domain",
    "tier",
)
T4_COLUMNS: Final[tuple[str, ...]] = (
    "canon_label",
    "chain",
    "domain",
    "n_serotypes_present",
    "n_serotypes_reproducible",
    "n_serotypes_signed_reproducible",
    "frac_reproducible",
    "conservation_class",
    "is_serotype_divergent",
    "is_catalytic_triad",
    "rho_residue_median",
    "tier",
)
T5_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "tau2",
    "sigma2_bar",
    "frac_tau2",
    "frac_sigma2",
    "tau2_sigma2_ratio",
    "variance_regime",
    "tier",
)
TABLE_COLUMNS: Final[dict[str, tuple[str, ...]]] = {
    "T1": T1_COLUMNS,
    "T2": T2_COLUMNS,
    "T3": T3_COLUMNS,
    "T4": T4_COLUMNS,
    "T5": T5_COLUMNS,
}

# ===========================================================================
# Prepared-figure-data column schemas (plotting-ready tables the builders emit)
# ===========================================================================
F1_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "rho_residue",
    "tier",
)
F2_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "gated_scale_level",
    "gated_scale_index",
    "n_loci",
)
F3_COLUMNS: Final[tuple[str, ...]] = (
    "domain",
    "serotype",
    "chain",
    "rho_domain",
    "is_catalytic_domain",
    "tier",
)
F4_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "domain",
    "beta_signed",
    "beta_ci_lower",
    "beta_ci_upper",
    "significant_fdr",
    "tier",
)
F5_COLUMNS: Final[tuple[str, ...]] = (
    "canon_label",
    "chain",
    "domain",
    "n_serotypes_reproducible",
    "frac_reproducible",
    "conservation_class",
    "is_serotype_divergent",
    "tier",
)
F6_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "frac_tau2",
    "frac_sigma2",
    "variance_regime",
    "tier",
)
F7_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "domain",
    "scale_index",
    "scale_level",
    "rho",
    "tier",
)
F8_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "rho_domain",
    "coherence_domain",
    "is_coherent",
    "tier",
)
FIGURE_DATA_COLUMNS: Final[dict[str, tuple[str, ...]]] = {
    "F1": F1_COLUMNS,
    "F2": F2_COLUMNS,
    "F3": F3_COLUMNS,
    "F4": F4_COLUMNS,
    "F5": F5_COLUMNS,
    "F6": F6_COLUMNS,
    "F7": F7_COLUMNS,
    "F8": F8_COLUMNS,
}

# ===========================================================================
# Conservation-class ordering / palette keys (stable categorical order)
# ===========================================================================
CONSERVATION_CLASS_ORDER: Final[tuple[str, ...]] = (
    "reproducible_all",
    "reproducible_majority",
    "reproducible_some",
    "reproducible_none",
)

# ===========================================================================
# Output filenames
# ===========================================================================
OUT_SUMMARY: Final[str] = "report_summary.json"
OUT_MANIFEST: Final[str] = "artifact_manifest.parquet"

# ---- figure/table artifact suffixes --------------------------------------
SVG_SUFFIX: Final[str] = ".svg"
CSV_SUFFIX: Final[str] = ".csv"
PARQUET_SUFFIX: Final[str] = ".parquet"
MD_SUFFIX: Final[str] = ".md"

# ===========================================================================
# SVG style constants (minimal, fixed — no unnecessary styling)
# ===========================================================================
SVG_WIDTH: Final[int] = 720
SVG_HEIGHT: Final[int] = 420
SVG_MARGIN_LEFT: Final[int] = 70
SVG_MARGIN_RIGHT: Final[int] = 24
SVG_MARGIN_TOP: Final[int] = 48
SVG_MARGIN_BOTTOM: Final[int] = 64
SVG_FONT_FAMILY: Final[str] = "monospace"
SVG_FONT_SIZE: Final[int] = 12
SVG_TITLE_FONT_SIZE: Final[int] = 15
SVG_AXIS_COLOR: Final[str] = "#333333"
SVG_GRID_COLOR: Final[str] = "#dddddd"
SVG_TEXT_COLOR: Final[str] = "#111111"
SVG_MARK_COLOR: Final[str] = "#1f5c99"
SVG_MARK_COLOR_ALT: Final[str] = "#b8531a"
SVG_BACKGROUND: Final[str] = "#ffffff"

#: A fixed, colour-blind-friendly categorical palette (stable order), used for
#: per-serotype series and categorical fills. Deterministic assignment by index.
SVG_CATEGORICAL_PALETTE: Final[tuple[str, ...]] = (
    "#1f5c99",  # blue
    "#b8531a",  # orange
    "#2e7d32",  # green
    "#6a1b9a",  # purple
    "#8d6e00",  # olive
    "#00838f",  # teal
)

#: A fixed 5-stop sequential ramp (low→high) for the ρ heatmap fills.
SVG_SEQUENTIAL_RAMP: Final[tuple[str, ...]] = (
    "#f7fbff",
    "#c6dbef",
    "#6baed6",
    "#2171b5",
    "#08306b",
)
