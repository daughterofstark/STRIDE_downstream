"""Task 1 — canonical residue identifiers.

Constructs one canonical residue object per (serotype, canon_label) from the
residue-scale rows of the S0 STRIDE table. Each object carries the residue's
structural annotations (chain, domain, motif, secondary structure) and its full
hierarchy path. No interpretation: the hierarchy values (including STRIDE's
``unassigned``/``none``/``unknown`` sentinels) are preserved verbatim.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    CANONICAL_RESIDUES_COLUMNS,
    RESIDUE_SCALE_LEVEL,
)


def build_canonical_residues(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the canonical residues table from the STRIDE table.

    One row per (serotype, canon_label). Raises :class:`ConsistencyError` if
    that key is not unique among residue-scale rows, or if any residue carries
    more than one distinct hierarchy path.
    """
    res = stride_table[
        stride_table["scale_level"] == RESIDUE_SCALE_LEVEL
    ].copy()

    if res.empty:
        return pd.DataFrame(columns=list(CANONICAL_RESIDUES_COLUMNS))

    # key uniqueness (each residue appears once per serotype at residue scale)
    dup = res.duplicated(["serotype", "canon_label"], keep=False)
    if dup.any():
        examples = (
            res.loc[dup, ["serotype", "canon_label"]]
            .drop_duplicates()
            .head(5)
            .to_dict("records")
        )
        raise ConsistencyError(
            f"canonical residue key (serotype, canon_label) not unique; "
            f"examples: {examples}"
        )

    out = pd.DataFrame(
        {
            "serotype": res["serotype"].to_numpy(),
            "canon_label": res["canon_label"].to_numpy(),
            "chain": res["h_chain"].to_numpy(),
            "domain": res["h_domain"].to_numpy(),
            "motif": res["h_motif"].to_numpy(),
            "secondary_structure": res["h_secondary_structure"].to_numpy(),
            "hierarchy_path": res["region_id"].to_numpy(),
            "complex": res["h_complex"].to_numpy(),
            "protein": res["h_protein"].to_numpy(),
            "residue": res["h_residue"].to_numpy(),
        }
    )

    # each (serotype, canon_label) must have exactly one hierarchy path
    path_counts = out.groupby(["serotype", "canon_label"])[
        "hierarchy_path"
    ].nunique()
    multi = path_counts[path_counts > 1]
    if len(multi):
        raise ConsistencyError(
            f"{len(multi)} residue(s) carry more than one hierarchy path, e.g. "
            f"{list(multi.index[:3])}"
        )

    out = out[list(CANONICAL_RESIDUES_COLUMNS)]
    return out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
