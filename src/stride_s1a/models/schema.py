"""Frozen schema constants for Stage S1A outputs.

S1A consumes **only** the S0 canonical tables (``replicate_table.parquet`` and
``stride_table.parquet``) and produces reusable biological-data-layer tables.
This module is the single source of truth for the shape of everything S1A emits,
and for the S0 columns it depends on.

S1A contains no statistics, no scoring, no interpretation — every column here is
a structural annotation or an availability count.
"""
from __future__ import annotations

from typing import Final

# ===========================================================================
# S0 columns S1A depends on (a subset of the S0 frozen schema)
# ===========================================================================
#: The residue-scale row of the STRIDE table carries the full hierarchy for a
#: residue. S1A reads these columns and no others from the STRIDE table.
STRIDE_HIERARCHY_COLUMNS: Final[tuple[str, ...]] = (
    "h_complex",
    "h_protein",
    "h_chain",
    "h_domain",
    "h_motif",
    "h_secondary_structure",
    "h_residue",
)
#: Columns required to exist on the STRIDE table for S1A to run.
STRIDE_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "scale_level",
    "locus",
    "region_id",
    *STRIDE_HIERARCHY_COLUMNS,
)
#: Columns required to exist on the replicate table for S1A to run.
REPLICATE_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "replicate",
    "replicate_index",
    "canon_label",
)

#: The scale level whose rows define a residue (identity grain).
RESIDUE_SCALE_LEVEL: Final[str] = "residue"

#: The four dengue serotypes (used only for ordering/reporting, never hardcoded
#: as a requirement — S1A discovers whatever serotypes the tables contain).
DENGUE_SEROTYPES: Final[tuple[str, ...]] = ("DENV1", "DENV2", "DENV3", "DENV4")


# ===========================================================================
# Output table 1 — canonical_residues.parquet
# ===========================================================================
#: One row per (serotype, canon_label): the canonical residue object.
#: Unique key: (serotype, canon_label).
CANONICAL_RESIDUES_COLUMNS: Final[tuple[str, ...]] = (
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
CANONICAL_RESIDUES_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")


# ===========================================================================
# Output table 2 — domain_table.parquet
# ===========================================================================
#: One row per (serotype, chain, domain): structural summary of a domain.
#: Unique key: (serotype, chain, domain). Structural counts only — no scoring.
DOMAIN_TABLE_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "complex",
    "protein",
    "chain",
    "domain",
    "n_residues",
    "canon_labels",       # sorted list of member canon_labels (identifiers)
    "hierarchy_path",     # complex/protein/chain/domain
)
DOMAIN_TABLE_KEY: Final[tuple[str, ...]] = ("serotype", "chain", "domain")


# ===========================================================================
# Output table 3 — replicate_inventory.parquet
# ===========================================================================
#: One row per (serotype, canon_label): replicate availability for a residue.
#: Unique key: (serotype, canon_label). Availability only — never averages.
REPLICATE_INVENTORY_COLUMNS: Final[tuple[str, ...]] = (
    "serotype",
    "canon_label",
    "n_replicates",       # how many replicates observed this residue
    "replicates",         # sorted list of replicate names (e.g. 1st_run)
    "replicate_indices",  # sorted list of 1-based replicate indices
    "available",          # bool: n_replicates > 0
    "in_all_replicates",  # bool: observed in every replicate of the serotype
)
REPLICATE_INVENTORY_KEY: Final[tuple[str, ...]] = ("serotype", "canon_label")


# ===========================================================================
# Output table 4 — conservation_table.parquet
# ===========================================================================
#: One row per canon_label (across serotypes): cross-serotype presence map.
#: Unique key: (canon_label,). Set membership only — no scoring, no alignment.
CONSERVATION_TABLE_COLUMNS: Final[tuple[str, ...]] = (
    "canon_label",
    "n_serotypes",           # in how many serotypes this residue is present
    "serotypes_present",     # sorted list of serotypes containing it
    "serotypes_absent",      # sorted list of serotypes lacking it
    "in_all_serotypes",      # bool: present in every serotype (intersection)
    "in_any_serotype",       # bool: always True by construction (union member)
    "chain",                 # structural annotation (consistent across serotypes)
    "domain",
)
CONSERVATION_TABLE_KEY: Final[tuple[str, ...]] = ("canon_label",)


# ===========================================================================
# output artifact filenames
# ===========================================================================
OUT_CANONICAL_RESIDUES: Final[str] = "canonical_residues.parquet"
OUT_DOMAIN_TABLE: Final[str] = "domain_table.parquet"
OUT_REPLICATE_INVENTORY: Final[str] = "replicate_inventory.parquet"
OUT_CONSERVATION_TABLE: Final[str] = "conservation_table.parquet"
OUT_DATASET_SUMMARY: Final[str] = "dataset_summary.json"

#: Sentinel values used by STRIDE for unresolved hierarchy levels. S1A preserves
#: them verbatim (it does not clean or reinterpret them — that would be
#: interpretation). Listed here only for documentation/reporting.
UNRESOLVED_HIERARCHY_SENTINELS: Final[frozenset[str]] = frozenset(
    {"unassigned", "none", "unknown"}
)
