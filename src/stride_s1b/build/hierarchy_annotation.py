"""Hierarchy annotation.

Builds one annotation row per (serotype, hierarchy_path) from the S1A canonical
residues. Because a residue's ``hierarchy_path`` is unique to that residue, this
table annotates the *resolution completeness* of each residue's structural path:
how many of the seven hierarchy levels are real assignments versus STRIDE
sentinels.

No statistics, no scoring, no ranking.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    HIERARCHY_ANNOTATION_COLUMNS,
    N_HIERARCHY_LEVELS,
)
from ._classify import is_resolved


def build_hierarchy_annotation(canonical_residues: pd.DataFrame) -> pd.DataFrame:
    """Build the hierarchy annotation table.

    One row per (serotype, hierarchy_path). Raises :class:`ConsistencyError` if
    that key is not unique (a hierarchy path must identify a single residue
    within a serotype).
    """
    if canonical_residues.empty:
        return pd.DataFrame(columns=list(HIERARCHY_ANNOTATION_COLUMNS))

    dup = canonical_residues.duplicated(
        ["serotype", "hierarchy_path"], keep=False
    )
    if dup.any():
        examples = (
            canonical_residues.loc[dup, ["serotype", "hierarchy_path"]]
            .drop_duplicates()
            .head(5)
            .to_dict("records")
        )
        raise ConsistencyError(
            f"hierarchy key (serotype, hierarchy_path) not unique; "
            f"examples: {examples}"
        )

    records = []
    for serotype, hierarchy_path, canon_label, complex_, protein, chain, domain, motif, ss, residue in zip(
        canonical_residues["serotype"],
        canonical_residues["hierarchy_path"],
        canonical_residues["canon_label"],
        canonical_residues["complex"],
        canonical_residues["protein"],
        canonical_residues["chain"],
        canonical_residues["domain"],
        canonical_residues["motif"],
        canonical_residues["secondary_structure"],
        canonical_residues["residue"],
        strict=True,
    ):
        n_resolved = sum(
            1
            for v in (complex_, protein, chain, domain, motif, ss, residue)
            if is_resolved(v)
        )
        records.append(
            {
                "serotype": serotype,
                "hierarchy_path": hierarchy_path,
                "canon_label": canon_label,
                "complex": complex_,
                "protein": protein,
                "chain": chain,
                "domain": domain,
                "motif": motif,
                "secondary_structure": ss,
                "n_levels_total": N_HIERARCHY_LEVELS,
                "n_levels_resolved": n_resolved,
                "fully_resolved": bool(n_resolved == N_HIERARCHY_LEVELS),
            }
        )

    out = pd.DataFrame.from_records(records)[list(HIERARCHY_ANNOTATION_COLUMNS)]
    return out.sort_values(["serotype", "hierarchy_path"]).reset_index(
        drop=True
    )
