"""Frozen schema constants for STRIDE data levels.

This module is the single source of truth for the *shape* of everything the
framework ingests. Two data levels are described (Level 0 raw MD inputs are never
read by this framework):

- **Level 1 — replicate observations** (``*_correlations_v5.csv``): the primary
  per-replicate, per-residue measurements written by the STRIDE production
  engine.
- **Level 2 — STRIDE summaries** (``*_profile.csv`` / ``*_mechanism.json``): the
  triplicate summaries derived from Level 1.

Every constant here was cross-checked against the STRIDE production source (which
emits Level 1) and against real Level-2 outputs (DENV1–4). The framework
validates incoming files against these constants and refuses to guess on drift.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# Level 1 — replicate observations: *_correlations_v5.csv
# ===========================================================================
# The STRIDE engine writes a wide, additive schema. Downstream analysis only
# *needs* a stable core; later milestones (M1/M2/M3) append more columns. We
# therefore split the schema into a REQUIRED core (validated strictly) and a
# set of KNOWN-optional columns (validated when present, tolerated when absent).

#: Columns any replicate table MUST contain (the residue identity + the effect
#: field). Missing any of these is a hard error.
CORRELATIONS_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "file_resid",     # int  — residue id as it appears in the structure file
    "canon_resid",    # int  — canonical residue number (file_resid - offset)
    "name",           # str  — residue name (e.g. HIS)
    "label",          # str  — canonical label, e.g. "NS3:51" (join key to L2)
    "r",              # float — signed effect θ (the per-replicate effect field)
    "abs_r",          # float — |θ|
)

#: Known-optional columns the STRIDE engine may append (M0 stats, M1 n_eff,
#: M2 bootstrap, M3 structural metadata). Validated for type when present.
CORRELATIONS_KNOWN_OPTIONAL_COLUMNS: Final[tuple[str, ...]] = (
    # M0 base statistics
    "domain_label", "p_raw", "p_bonf", "sig_bonf", "abs_r_com",
    "block_std", "block_sem", "block_ci_hw", "converged", "rmsf",
    "p_fdr", "sig_fdr",
    # M1 effective sample size / uncertainty
    "tau_int", "n_eff", "neff_status", "theta_se",
    # M2 block-bootstrap CIs
    "theta_bootstrap_se", "theta_bootstrap_ci_lower", "theta_bootstrap_ci_upper",
    # M3 structural metadata
    "chain", "domain", "motif", "secondary_structure", "region_id",
)

#: Replicate columns that, when present, must be numeric.
CORRELATIONS_FLOAT_COLUMNS: Final[frozenset[str]] = frozenset({
    "r", "abs_r", "p_raw", "p_bonf", "abs_r_com", "block_std", "block_sem",
    "block_ci_hw", "rmsf", "p_fdr", "tau_int", "n_eff", "theta_se",
    "theta_bootstrap_se", "theta_bootstrap_ci_lower", "theta_bootstrap_ci_upper",
})

#: Replicate columns that, when present, must be integer-valued.
CORRELATIONS_INT_COLUMNS: Final[frozenset[str]] = frozenset({
    "file_resid", "canon_resid",
})

#: Replicate string columns that must be non-empty when present.
CORRELATIONS_STR_COLUMNS: Final[frozenset[str]] = frozenset({
    "name", "label", "domain_label", "neff_status", "chain", "domain",
    "motif", "secondary_structure", "region_id",
})

#: Suffix used to discover a replicate table inside a serotype's analysis dir.
CORRELATIONS_SUFFIX: Final[str] = "_correlations_v5.csv"

#: Directory (inside each ``<run>/<serotype>/``) that holds the replicate table.
ANALYSIS_OUTPUT_DIRNAME: Final[str] = "analysis_output"

#: Default run-directory names (order defines replicate index 1..K). Discovery
#: does not hardcode these — it discovers whatever run dirs exist — but they are
#: the conventional names and are used by the synthetic example.
DEFAULT_RUN_DIR_NAMES: Final[tuple[str, ...]] = ("1st_run", "2nd_run", "3rd_run")


# ===========================================================================
# Level 2 — STRIDE summaries: *_profile.csv
# ===========================================================================
PROFILE_COLUMNS: Final[tuple[str, ...]] = (
    "protein", "locus", "canon_label", "scale_index", "scale_level",
    "region_id", "rho", "gated", "beta", "beta_se", "tau2", "sigma2_bar",
    "a_signed", "coherence", "method", "status",
)
PROFILE_FLOAT_COLUMNS: Final[tuple[str, ...]] = (
    "rho", "beta", "beta_se", "tau2", "sigma2_bar", "a_signed", "coherence",
)
PROFILE_INT_COLUMNS: Final[tuple[str, ...]] = ("scale_index",)
PROFILE_BOOL_COLUMNS: Final[tuple[str, ...]] = ("gated",)
PROFILE_STR_COLUMNS: Final[tuple[str, ...]] = (
    "protein", "locus", "canon_label", "scale_level", "region_id",
    "method", "status",
)
PROFILE_SUFFIX: Final[str] = "_profile.csv"


# ===========================================================================
# Level 2 — STRIDE summaries: *_mechanism.json
# ===========================================================================
MECHANISM_TOP_KEYS: Final[tuple[str, ...]] = (
    "schema_version", "calibrated", "uncalibrated_note", "gate", "summary",
    "mechanisms", "unresolved_loci",
)
MECHANISM_GATE_KEYS: Final[tuple[str, ...]] = (
    "rho_star", "alpha", "coherence_threshold",
)
MECHANISM_SUMMARY_KEYS: Final[tuple[str, ...]] = (
    "n_loci", "n_mechanisms", "n_unresolved", "n_gate_uncertain",
)
MECHANISM_MECH_FIELDS: Final[tuple[str, ...]] = (
    "region_id", "label", "scale_level", "scale_index", "n_loci", "loci",
    "rho", "rho_star", "calibrated", "direction", "beta_signed",
    "beta_ci_lower", "beta_ci_upper", "beta_se", "coherence",
    "reproducible_magnitude_energy", "method", "gate_uncertain", "status",
)
MECHANISM_NULLABLE_WHEN_MIXED: Final[tuple[str, ...]] = (
    "beta_signed", "beta_ci_lower", "beta_ci_upper", "beta_se",
)
VALID_DIRECTIONS: Final[frozenset[str]] = frozenset(
    {"increase", "decrease", "mixed"}
)
MECHANISM_SUFFIX: Final[str] = "_mechanism.json"


# ===========================================================================
# hierarchy grammar (shared by Level 2)
# ===========================================================================
HIERARCHY_LEVELS: Final[tuple[str, ...]] = (
    "complex", "protein", "chain", "domain", "motif",
    "secondary_structure", "residue",
)
SCALE_INDEX_TO_LEVEL: Final[dict[int, str]] = {
    0: "residue", 1: "secondary_structure", 2: "motif", 3: "domain",
    4: "chain", 5: "protein", 6: "complex",
}
SCALE_LEVEL_TO_INDEX: Final[dict[str, int]] = {
    v: k for k, v in SCALE_INDEX_TO_LEVEL.items()
}
N_SCALES: Final[int] = len(HIERARCHY_LEVELS)
SCALE_LEVEL_PATH_DEPTH: Final[dict[str, int]] = {
    "complex": 1, "protein": 2, "chain": 3, "domain": 4, "motif": 5,
    "secondary_structure": 6, "residue": 7,
}
HIERARCHY_COLUMNS: Final[tuple[str, ...]] = tuple(
    f"h_{level}" for level in HIERARCHY_LEVELS
)


# ===========================================================================
# canonical tables — the two frozen output schemas
# ===========================================================================
#: Replicate canonical table. Unique key: (serotype, replicate, canon_label).
#: Preserves the original per-run quantities; never merged with Level 2.
REPLICATE_TABLE_IDENTITY: Final[tuple[str, ...]] = (
    "serotype", "replicate", "replicate_index", "canon_label",
    "file_resid", "canon_resid", "name",
)
REPLICATE_TABLE_PROVENANCE: Final[tuple[str, ...]] = (
    "source_path", "run_dir",
)

#: STRIDE canonical table. Unique key: (serotype, canon_label, scale_level).
MECHANISM_MERGE_COLUMNS: Final[tuple[str, ...]] = (
    "mech_label", "mech_direction", "mech_beta_signed", "mech_beta_ci_lower",
    "mech_beta_ci_upper", "mech_beta_se", "mech_coherence",
    "mech_reproducible_magnitude_energy", "mech_rho_star", "mech_calibrated",
    "mech_gate_uncertain", "mech_status", "mech_region_id", "mech_n_loci",
)
STRIDE_TABLE_COLUMNS: Final[tuple[str, ...]] = (
    # identity
    "serotype", "canon_label", "scale_level", "scale_index", "locus",
    "region_id",
    # profile payload
    "rho", "gated", "beta", "beta_se", "tau2", "sigma2_bar", "a_signed",
    "coherence", "method", "status",
    # explicit hierarchy
    *HIERARCHY_COLUMNS,
    # mechanism payload (gated rows only)
    "is_gated_scale", *MECHANISM_MERGE_COLUMNS,
    # provenance
    "profile_source", "mechanism_source", "gate_rho_star", "gate_alpha",
    "gate_coherence_threshold", "mechanism_calibrated",
    "mechanism_schema_version",
)

#: Numerical tolerance for float equality in consistency checks.
FLOAT_TOL: Final[float] = 1e-9


# ===========================================================================
# output artifact filenames
# ===========================================================================
OUT_STRIDE_PARQUET: Final[str] = "stride_table.parquet"
OUT_STRIDE_CSV: Final[str] = "stride_table.csv"
OUT_REPLICATE_PARQUET: Final[str] = "replicate_table.parquet"
OUT_REPLICATE_CSV: Final[str] = "replicate_table.csv"
OUT_SCHEMA_REPORT: Final[str] = "schema_report.json"
OUT_VALIDATION_REPORT: Final[str] = "validation_report.md"
