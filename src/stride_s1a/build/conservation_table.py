"""Task 2 — conserved residue mapping.

Builds a cross-serotype presence map keyed by ``canon_label``: which serotypes
contain each residue, which lack it, plus the union/intersection flags. This is
pure set membership — no sequence alignment, no scoring, no interpretation of
why a residue is present or absent.

The structural annotations (chain, domain) are carried through for convenience;
they are consistent across serotypes for a shared ``canon_label`` (S1A validates
this in :mod:`stride_s1a.validation`).
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import CONSERVATION_TABLE_COLUMNS


def build_conservation_table(canonical_residues: pd.DataFrame) -> pd.DataFrame:
    """Build the conservation table from the canonical residues table.

    One row per distinct ``canon_label`` across all serotypes (the union).
    """
    if canonical_residues.empty:
        return pd.DataFrame(columns=list(CONSERVATION_TABLE_COLUMNS))

    all_serotypes = sorted(canonical_residues["serotype"].unique().tolist())
    present_by_label: dict[str, set[str]] = {}
    annot_by_label: dict[str, tuple[str, str]] = {}

    for serotype, canon_label, chain, domain in zip(
        canonical_residues["serotype"],
        canonical_residues["canon_label"],
        canonical_residues["chain"],
        canonical_residues["domain"],
        strict=True,
    ):
        present_by_label.setdefault(canon_label, set()).add(serotype)
        # first-seen annotation; validation guarantees consistency across serotypes
        annot_by_label.setdefault(canon_label, (chain, domain))

    records = []
    for canon_label in sorted(present_by_label):
        present = sorted(present_by_label[canon_label])
        absent = sorted(set(all_serotypes) - set(present))
        chain, domain = annot_by_label[canon_label]
        records.append(
            {
                "canon_label": canon_label,
                "n_serotypes": len(present),
                "serotypes_present": present,
                "serotypes_absent": absent,
                "in_all_serotypes": len(present) == len(all_serotypes),
                "in_any_serotype": len(present) > 0,
                "chain": chain,
                "domain": domain,
            }
        )

    out = pd.DataFrame.from_records(records)[list(CONSERVATION_TABLE_COLUMNS)]
    return out.sort_values("canon_label").reset_index(drop=True)
