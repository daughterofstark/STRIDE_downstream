r"""Frozen schema constants for Stage S5 outputs.

S5 is the **cross-serotype layer** (serotype is the unit of biological
replication, n = 4). It consumes the S0 STRIDE table (the tidy master profile:
every locus at every one of the seven hierarchy scales, carrying ρ, the gated
mechanism direction, and the variance components) together with the S1A
``conservation_table`` (the shared-position index: which serotypes contain each
``canon_label``). It then compares the four serotypes at the level of shared
canonical positions and named regions:

- **conservation of reproducibility** across the shared ``canon_label`` positions
  — reproducible in all / a majority / some / none of the serotypes it is present
  in (design §3.3);
- **serotype-divergent positions** — signed and reproducible in some serotypes
  but absent or non-signed in others (candidate serotype-specific differences,
  §3.3);
- **direction concordance** — for shared positions that are signed in ≥ 2
  serotypes, whether the serotypes agree on increase vs decrease (§3.3);
- the **domain × serotype ρ matrix** (tidy long) over the NS3 domains + NS2B, with
  the conserved catalytic machinery (Catalytic Triad, Oxyanion Loop) flagged
  (§3.3); and
- a **per-serotype cross-serotype scorecard** (n_loci, %residue-gated, %signed,
  %mixed, ρ median, shared-position counts; §3.6).

This module is the single source of truth for:

- the S0 STRIDE-table and S1A conservation-table columns S5 depends on,
- the provisional gate constant ρ\* and the two-tier (licensed / exploratory)
  vocabulary,
- the conserved catalytic machinery labels,
- the conservation- and concordance-class vocabularies, and
- the four output-table schemas, their keys, and the artifact filenames.

Every design decision here follows the authoritative design document:

- **serotype is the unit of replication (n = 4)** — S5 aggregates to one value per
  serotype per region *first*, then compares across the four serotypes; it never
  treats residues as independent samples and prefers descriptive statistics and
  effect sizes over p-values at n = 4 (§5.2, table §2.3 pseudoreplication note);
- **the gate is uncalibrated** — every reproducibility statement is descriptive
  and relative to the provisional ρ\* = 0.5, clearly labelled provisional; S5
  makes no calibrated pass/fail claim (§0.1, §5.3);
- ρ and the variance components are **read**, never recomputed (§1.2);
- **domain scale is the licensed claim level at K = 3**; residue-scale products
  (position conservation, direction concordance, the scorecard's residue fields)
  are labelled exploratory / outside the operating range (§0.1, §5.3, §5.4).

S5 reads only the S0 STRIDE table and the S1A conservation table; it never
consumes the S2/S3/S4 reduction outputs, never re-reads the raw STRIDE files,
never produces figures, and never modifies any earlier stage.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S0 STRIDE-table columns S5 depends on (a subset of the S0 frozen schema)
# ===========================================================================
#: Columns S5 reads from the S0 STRIDE table. The residue-scale rows carry ρ and
#: the structural labels; the gated row of each locus carries the mechanism
#: direction; the domain-scale rows carry the region-constant ρ/β/variance
#: components used by the domain × serotype matrix.
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
    "h_chain",
    "h_domain",
    "is_gated_scale",
    "mech_direction",
)

# ===========================================================================
# S1A conservation-table columns S5 depends on (the shared-position index)
# ===========================================================================
#: Columns S5 reads from the S1A ``conservation_table`` — the authoritative
#: cross-serotype presence map keyed by ``canon_label``.
CONSERVATION_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "canon_label",
    "n_serotypes",
    "serotypes_present",
    "in_all_serotypes",
    "chain",
    "domain",
)

# ===========================================================================
# Scale grammar (mirrors the S0 grammar; re-declared so S5 is self-contained)
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
#: design is emphatic that this is *provisional* and uncalibrated (§0.1). S5 uses
#: it descriptively to decide whether a position's residue-scale ρ clears the
#: threshold ("reproducible" at the provisional gate); it makes no calibrated
#: claim. A locus's residue-scale ρ ≥ ρ\* is exactly equivalent to its being
#: gated at the residue scale (the finest scale where ρ ≥ ρ\*).
PROVISIONAL_RHO_STAR: Final[float] = 0.5

#: Number of decimals ρ/β/variance-derived quantities are rounded to (keeps
#: emitted rows deterministic and free of binary-float noise).
RHO_DECIMALS: Final[int] = 6

# ===========================================================================
# Two-tier output vocabulary (design §5.4)
# ===========================================================================
#: Tier A — licensed at K = 3: domain-scale and coarser.
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
#: Placeholder direction for a position that is not reproducible at the residue
#: scale (so it carries no residue-scale signed claim).
DIRECTION_NONE: Final[str] = "none"
DIRECTIONS: Final[tuple[str, ...]] = (
    DIRECTION_INCREASE,
    DIRECTION_DECREASE,
    DIRECTION_MIXED,
)

# ===========================================================================
# Conserved catalytic machinery (design §3.3)
# ===========================================================================
#: The Catalytic Triad residues, by canonical label. Biologically these should be
#: the most conserved across serotypes; the position tables flag them.
CATALYTIC_TRIAD_LABELS: Final[frozenset[str]] = frozenset(
    {"NS3:51", "NS3:75", "NS3:135"}
)
#: The catalytic-machinery domains flagged in the domain × serotype matrix.
CATALYTIC_DOMAINS: Final[frozenset[str]] = frozenset(
    {"Catalytic Triad", "Oxyanion Loop"}
)

# ===========================================================================
# Conservation-of-reproducibility vocabulary (design §3.3)
# ===========================================================================
#: A shared position is reproducible in *all* the serotypes it is present in …
CONSERVATION_ALL: Final[str] = "reproducible_all"
#: … in a strict majority of them …
CONSERVATION_MAJORITY: Final[str] = "reproducible_majority"
#: … in at least one but not a majority …
CONSERVATION_SOME: Final[str] = "reproducible_some"
#: … or in none of them.
CONSERVATION_NONE: Final[str] = "reproducible_none"
CONSERVATION_CLASSES: Final[tuple[str, ...]] = (
    CONSERVATION_ALL,
    CONSERVATION_MAJORITY,
    CONSERVATION_SOME,
    CONSERVATION_NONE,
)

# ===========================================================================
# Direction-concordance vocabulary (design §3.3)
# ===========================================================================
#: All signed serotypes push the same way (increase-only or decrease-only).
CONCORDANCE_AGREE: Final[str] = "agree"
#: Both directions occur but one is a strict majority.
CONCORDANCE_MAJORITY: Final[str] = "majority"
#: The signed serotypes split evenly (a tie between increase and decrease).
CONCORDANCE_CONFLICT: Final[str] = "conflict"
CONCORDANCE_CLASSES: Final[tuple[str, ...]] = (
    CONCORDANCE_AGREE,
    CONCORDANCE_MAJORITY,
    CONCORDANCE_CONFLICT,
)
#: A shared position must be signed-and-reproducible in at least this many
#: serotypes for direction concordance to be well-posed.
MIN_SEROTYPES_FOR_CONCORDANCE: Final[int] = 2

# ===========================================================================
# Output table 1 — position_conservation.parquet  (Tier B — exploratory)
# ===========================================================================
#: One row per shared ``canon_label`` (every position in the S1A conservation
#: index): conservation of reproducibility across the serotypes it is present in.
#: ``n_serotypes_reproducible`` counts serotypes whose residue-scale ρ ≥ ρ\*;
#: ``conservation_class`` labels all / majority / some / none; a divergent flag
#: marks positions signed-and-reproducible in some but not all serotypes.
#: Residue-scale ⇒ Tier B (exploratory). Unique key: (canon_label,).
POSITION_CONSERVATION_COLUMNS: Final[tuple[str, ...]] = (
    "canon_label",
    "chain",
    "domain",
    "n_serotypes_total",             # serotypes in the dataset (denominator context)
    "n_serotypes_present",           # from the S1A conservation table
    "serotypes_present",             # sorted list, from the S1A conservation table
    "in_all_serotypes",              # present in every serotype (S1A flag)
    "n_serotypes_reproducible",      # present AND residue ρ ≥ ρ*
    "n_serotypes_signed_reproducible",  # present AND reproducible AND signed
    "frac_reproducible",             # reproducible / present (NaN if present == 0)
    "conservation_class",            # reproducible_all | majority | some | none
    "is_serotype_divergent",         # signed+reproducible in some but not all serotypes
    "is_catalytic_triad",            # canon_label in the Catalytic Triad
    "rho_residue_min",               # residue-scale ρ across serotypes present
    "rho_residue_median",
    "rho_residue_max",
    "rho_star",
    "is_provisional_rho_star",
    "tier",                          # always exploratory
)
POSITION_CONSERVATION_KEY: Final[tuple[str, ...]] = ("canon_label",)

# ===========================================================================
# Output table 2 — direction_concordance.parquet  (Tier B — exploratory)
# ===========================================================================
#: One row per shared ``canon_label`` that is signed-and-reproducible in at least
#: ``MIN_SEROTYPES_FOR_CONCORDANCE`` serotypes: do those serotypes agree on the
#: sign of the effect? ``concordance_class`` is agree / majority / conflict.
#: Residue-scale ⇒ Tier B (exploratory). Unique key: (canon_label,).
DIRECTION_CONCORDANCE_COLUMNS: Final[tuple[str, ...]] = (
    "canon_label",
    "chain",
    "domain",
    "n_serotypes_signed",            # present AND reproducible AND signed
    "n_increase",
    "n_decrease",
    "majority_direction",            # increase | decrease | none (tie)
    "frac_majority",                 # max(n_increase, n_decrease) / n_serotypes_signed
    "concordance_class",             # agree | majority | conflict
    "is_catalytic_triad",
    "rho_star",
    "is_provisional_rho_star",
    "tier",                          # always exploratory
)
DIRECTION_CONCORDANCE_KEY: Final[tuple[str, ...]] = ("canon_label",)

# ===========================================================================
# Output table 3 — domain_serotype_matrix.parquet  (Tier A — licensed)
# ===========================================================================
#: One row per (serotype, chain, domain): the domain-scale ρ (and β, β_se, and
#: the variance components), read from the domain region's region-constant
#: domain-scale profile row. This is the tidy-long form of the design's
#: ρ(domain × serotype) heatmap (§3.3, F3); the conserved catalytic-machinery
#: domains are flagged. Domain-scale ⇒ Tier A (licensed). Unique key:
#: (serotype, chain, domain).
DOMAIN_SEROTYPE_MATRIX_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
    "region_id",
    "rho_domain",
    "beta_domain",
    "beta_se_domain",
    "tau2_domain",
    "sigma2_bar_domain",
    "is_catalytic_domain",           # domain in {Catalytic Triad, Oxyanion Loop}
    "tier",                          # always licensed
)
DOMAIN_SEROTYPE_MATRIX_KEY: Final[tuple[str, ...]] = (
    "serotype",
    "chain",
    "domain",
)

# ===========================================================================
# Output table 4 — cross_serotype_scorecard.parquet  (per serotype)
# ===========================================================================
#: One row per serotype: the cross-serotype reproducibility scorecard of the
#: design (§3.6) — n_loci, the residue-gated fraction, the signed / mixed
#: mechanism composition, the residue-scale ρ median/min/max, and how many of the
#: pan-serotype shared positions are reproducible in this serotype. Its core
#: quantities are residue-scale ⇒ Tier B (exploratory). Unique key: (serotype,).
CROSS_SEROTYPE_SCORECARD_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "n_loci",                        # distinct residue loci (canon_labels)
    "n_reproducible_residue",        # residue ρ ≥ ρ*
    "frac_reproducible_residue",     # n_reproducible_residue / n_loci
    "n_mechanisms",                  # gated mechanisms (one per gated region)
    "n_signed",                      # mechanisms with direction != mixed
    "n_mixed",                       # mechanisms with direction == mixed
    "frac_signed",                   # n_signed / n_mechanisms
    "frac_mixed",                    # n_mixed / n_mechanisms
    "rho_residue_median",
    "rho_residue_min",
    "rho_residue_max",
    "n_shared_positions",            # present here AND pan-serotype (in_all_serotypes)
    "n_shared_reproducible",         # of those, reproducible here
    "rho_star",
    "is_provisional_rho_star",
    "tier",                          # always exploratory
)
CROSS_SEROTYPE_SCORECARD_KEY: Final[tuple[str, ...]] = ("serotype",)

# ===========================================================================
# input + output artifact filenames
# ===========================================================================
IN_STRIDE_TABLE: Final[str] = "stride_table.parquet"
IN_CONSERVATION_TABLE: Final[str] = "conservation_table.parquet"

OUT_POSITION_CONSERVATION: Final[str] = "position_conservation.parquet"
OUT_DIRECTION_CONCORDANCE: Final[str] = "direction_concordance.parquet"
OUT_DOMAIN_SEROTYPE_MATRIX: Final[str] = "domain_serotype_matrix.parquet"
OUT_CROSS_SEROTYPE_SCORECARD: Final[str] = "cross_serotype_scorecard.parquet"
OUT_CONSERVATION_SUMMARY: Final[str] = "conservation_summary.json"
