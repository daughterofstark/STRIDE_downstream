"""Task 3 — domain summary tables.

Constructs one structural summary row per (serotype, chain, domain) from the
canonical residues. Structural summaries only: residue count, hierarchy
membership, and the member canonical identifiers. No biological scoring, no
ranking, no interpretation.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import DOMAIN_TABLE_COLUMNS


def build_domain_table(canonical_residues: pd.DataFrame) -> pd.DataFrame:
    """Build the domain summary table from the canonical residues table.

    One row per (serotype, chain, domain). ``canon_labels`` is the sorted list
    of member residue identifiers; ``n_residues`` is their count.
    """
    if canonical_residues.empty:
        return pd.DataFrame(columns=list(DOMAIN_TABLE_COLUMNS))

    records = []
    group_cols = ["serotype", "complex", "protein", "chain", "domain"]
    for keys, grp in canonical_residues.groupby(group_cols, sort=True):
        serotype, complex_, protein, chain, domain = keys
        labels = sorted(grp["canon_label"].tolist())
        records.append(
            {
                "serotype": serotype,
                "complex": complex_,
                "protein": protein,
                "chain": chain,
                "domain": domain,
                "n_residues": len(labels),
                "canon_labels": labels,
                "hierarchy_path": f"{complex_}/{protein}/{chain}/{domain}",
            }
        )

    out = pd.DataFrame.from_records(records)[list(DOMAIN_TABLE_COLUMNS)]
    return out.sort_values(["serotype", "chain", "domain"]).reset_index(
        drop=True
    )
