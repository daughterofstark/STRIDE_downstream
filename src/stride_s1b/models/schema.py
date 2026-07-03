"""Frozen schema constants for Stage S1B outputs.

S1B is the biological *annotation* layer. It consumes **only** the four S1A
parquet tables and derives reusable, deterministic structural annotations. It
performs no statistics, no ranking, no clustering, no hypothesis generation, and
produces no figures — every column here is a categorical or boolean label
derived by a fixed rule from S1A facts.

This module is the single source of truth for:
- the S1A columns S1B depends on,
- the deterministic annotation category vocabularies, and
- the four output-table schemas, keys, and filenames.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S1A columns S1B depends on (a subset of the S1A frozen schema)
# ===========================================================================
CANONICAL_RESIDUES_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "motif",
    "secondary_structure",
    "hierarchy_path",
    "complex",
    "protein",
    "residue",
)
DOMAIN_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "complex",
    "protein",
    "chain",
    "domain",
    "n_residues",
    "canon_labels",
    "hierarchy_path",
)
REPLICATE_INVENTORY_REQUIRED: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "n_replicates",
    "available",
    "in_all_replicates",
)
CONSERVATION_TABLE_REQUIRED: Final[tuple[str, ...]] = (
    "canon_label",
    "n_serotypes",
    "serotypes_present",
    "serotypes_absent",
    "in_all_serotypes",
    "in_any_serotype",
)

# ===========================================================================
# Deterministic annotation vocabularies
# ===========================================================================
# STRIDE marks an unresolved hierarchy level with one of these sentinels. S1B
# reads them to decide whether a level is "assigned" — it never rewrites them.
UNRESOLVED_HIERARCHY_SENTINELS: Final[frozenset[str]] = frozenset(
    {"unassigned", "none", "unknown", ""}
)

#: Domain-assignment status of a residue/domain.
DOMAIN_STATUS_ASSIGNED: Final[str] = "assigned"
DOMAIN_STATUS_UNASSIGNED: Final[str] = "unassigned"

#: Secondary-structure resolution status of a residue.
SS_STATUS_RESOLVED: Final[str] = "resolved"
SS_STATUS_UNRESOLVED: Final[str] = "unresolved"

#: Conservation class of a residue (from the S1A conservation table).
#: - ``pan_serotype``: present in every serotype (the intersection)
#: - ``partial``:      present in more than one but not all serotypes
#: - ``serotype_unique``: present in exactly one serotype
CONSERVATION_PAN: Final[str] = "pan_serotype"
CONSERVATION_PARTIAL: Final[str] = "partial"
CONSERVATION_UNIQUE: Final[str] = "serotype_unique"

#: Replicate-availability class of a residue (from the S1A replicate inventory).
#: - ``all_replicates``: observed in every replicate of its serotype
#: - ``some_replicates``: observed in at least one but not all replicates
#: - ``no_replicates``: not observed in any replicate (or summaries-only run)
AVAIL_ALL: Final[str] = "all_replicates"
AVAIL_SOME: Final[str] = "some_replicates"
AVAIL_NONE: Final[str] = "no_replicates"

#: Ordered category lists (documentation + validation of the closed vocabulary).
CONSERVATION_CLASSES: Final[tuple[str, ...]] = (
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
)
AVAILABILITY_CLASSES: Final[tuple[str, ...]] = (
    AVAIL_ALL,
    AVAIL_SOME,
    AVAIL_NONE,
)
DOMAIN_STATUSES: Final[tuple[str, ...]] = (
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
)
SS_STATUSES: Final[tuple[str, ...]] = (SS_STATUS_RESOLVED, SS_STATUS_UNRESOLVED)


# ===========================================================================
# Output table 1 — residue_annotation.parquet
# ===========================================================================
#: One row per (serotype, canon_label): the biological annotation of a residue.
#: Unique key: (serotype, canon_label). Every field is a deterministic label.
RESIDUE_ANNOTATION_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "hierarchy_path",
    "domain_status",          # assigned | unassigned
    "secondary_structure_status",  # resolved | unresolved
    "conservation_class",     # pan_serotype | partial | serotype_unique
    "n_serotypes_present",    # from conservation table
    "availability_class",     # all_replicates | some_replicates | no_replicates
    "n_replicates",           # from replicate inventory
)
RESIDUE_ANNOTATION_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")


# ===========================================================================
# Output table 2 — domain_annotation.parquet
# ===========================================================================
#: One row per (serotype, chain, domain): structural annotation of a domain.
#: Unique key: (serotype, chain, domain).
DOMAIN_ANNOTATION_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "complex",
    "protein",
    "chain",
    "domain",
    "hierarchy_path",
    "domain_status",              # assigned | unassigned
    "n_residues",                 # carried from S1A domain table
    "n_pan_serotype_residues",    # members conserved across all serotypes
    "fully_conserved",            # bool: every member residue is pan_serotype
)
DOMAIN_ANNOTATION_KEY: Final[tuple[str, ...]] = ("serotype", "chain", "domain")


# ===========================================================================
# Output table 3 — hierarchy_annotation.parquet
# ===========================================================================
#: One row per (serotype, hierarchy_path): annotation of a residue's structural
#: path. Because a residue's hierarchy_path is unique to that residue, this table
#: annotates path resolution completeness per residue path.
#: Unique key: (serotype, hierarchy_path).
HIERARCHY_ANNOTATION_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "hierarchy_path",
    "canon_label",
    "complex",
    "protein",
    "chain",
    "domain",
    "motif",
    "secondary_structure",
    "n_levels_total",       # always 7 (the grammar depth)
    "n_levels_resolved",    # how many of the 7 levels are non-sentinel
    "fully_resolved",       # bool: every level resolved
)
HIERARCHY_ANNOTATION_KEY: Final[tuple[str, ...]] = ("serotype", "hierarchy_path")


# ===========================================================================
# Output table 4 — serotype_annotation.parquet
# ===========================================================================
#: One row per serotype: dataset-level structural annotation.
#: Unique key: (serotype,).
SEROTYPE_ANNOTATION_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "n_residues",
    "n_domains",
    "n_assigned_domain_residues",
    "n_unassigned_domain_residues",
    "n_pan_serotype_residues",
    "n_partial_residues",
    "n_serotype_unique_residues",
    "n_residues_all_replicates",
    "n_residues_some_replicates",
    "n_residues_no_replicates",
)
SEROTYPE_ANNOTATION_KEY: Final[tuple[str, ...]] = ("serotype",)


# ===========================================================================
# hierarchy grammar depth (levels per residue path)
# ===========================================================================
HIERARCHY_LEVEL_COLUMNS: Final[tuple[str, ...]] = (
    "complex",
    "protein",
    "chain",
    "domain",
    "motif",
    "secondary_structure",
    "residue",
)
N_HIERARCHY_LEVELS: Final[int] = len(HIERARCHY_LEVEL_COLUMNS)


# ===========================================================================
# input + output artifact filenames
# ===========================================================================
IN_CANONICAL_RESIDUES: Final[str] = "canonical_residues.parquet"
IN_DOMAIN_TABLE: Final[str] = "domain_table.parquet"
IN_REPLICATE_INVENTORY: Final[str] = "replicate_inventory.parquet"
IN_CONSERVATION_TABLE: Final[str] = "conservation_table.parquet"

OUT_RESIDUE_ANNOTATION: Final[str] = "residue_annotation.parquet"
OUT_DOMAIN_ANNOTATION: Final[str] = "domain_annotation.parquet"
OUT_HIERARCHY_ANNOTATION: Final[str] = "hierarchy_annotation.parquet"
OUT_SEROTYPE_ANNOTATION: Final[str] = "serotype_annotation.parquet"
OUT_ANNOTATION_SUMMARY: Final[str] = "annotation_summary.json"
