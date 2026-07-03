"""Residue annotation.

Builds one biological annotation row per (serotype, canon_label) by joining the
S1A canonical residues with the conservation table and replicate inventory, then
applying deterministic classification rules (domain assignment status, secondary
structure resolution, conservation class, replicate availability class).

No statistics, no scoring, no ranking — every derived column is a fixed-rule
label.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import RESIDUE_ANNOTATION_COLUMNS
from ._classify import (
    availability_class,
    conservation_class,
    domain_status,
    secondary_structure_status,
)


def build_residue_annotation(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
    n_serotypes_total: int,
) -> pd.DataFrame:
    """Build the residue annotation table.

    Parameters
    ----------
    canonical_residues
        S1A canonical residues (one row per (serotype, canon_label)).
    conservation_table
        S1A conservation table (one row per canon_label), providing
        ``n_serotypes`` per residue.
    replicate_inventory
        S1A replicate inventory (one row per (serotype, canon_label)), providing
        ``n_replicates``, ``available`` and ``in_all_replicates``.
    n_serotypes_total
        Total number of serotypes in the dataset (denominator for the
        conservation class).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label), keyed and column-ordered per
        :data:`~stride_s1b.models.schema.RESIDUE_ANNOTATION_COLUMNS`.

    Raises
    ------
    ConsistencyError
        If the canonical residue key is not unique.
    """
    if canonical_residues.empty:
        return pd.DataFrame(columns=list(RESIDUE_ANNOTATION_COLUMNS))

    dup = canonical_residues.duplicated(["serotype", "canon_label"], keep=False)
    if dup.any():
        examples = (
            canonical_residues.loc[dup, ["serotype", "canon_label"]]
            .drop_duplicates()
            .head(5)
            .to_dict("records")
        )
        raise ConsistencyError(
            f"canonical residue key (serotype, canon_label) not unique; "
            f"examples: {examples}"
        )

    # lookup: canon_label -> n_serotypes present (from conservation table)
    n_serotypes_by_label: dict[str, int] = {}
    if not conservation_table.empty:
        n_serotypes_by_label = dict(
            zip(
                conservation_table["canon_label"],
                (int(x) for x in conservation_table["n_serotypes"]),
                strict=True,
            )
        )

    # lookup: (serotype, canon_label) -> availability facts
    avail_by_key: dict[tuple[str, str], tuple[int, bool, bool]] = {}
    if not replicate_inventory.empty:
        for serotype, canon_label, n_rep, avail, in_all in zip(
            replicate_inventory["serotype"],
            replicate_inventory["canon_label"],
            replicate_inventory["n_replicates"],
            replicate_inventory["available"],
            replicate_inventory["in_all_replicates"],
            strict=True,
        ):
            avail_by_key[(serotype, canon_label)] = (
                int(n_rep),
                bool(avail),
                bool(in_all),
            )

    records = []
    for serotype, canon_label, chain, domain, ss, path in zip(
        canonical_residues["serotype"],
        canonical_residues["canon_label"],
        canonical_residues["chain"],
        canonical_residues["domain"],
        canonical_residues["secondary_structure"],
        canonical_residues["hierarchy_path"],
        strict=True,
    ):
        n_present = n_serotypes_by_label.get(canon_label, 1)
        n_rep, avail, in_all = avail_by_key.get(
            (serotype, canon_label), (0, False, False)
        )
        records.append(
            {
                "serotype": serotype,
                "canon_label": canon_label,
                "chain": chain,
                "domain": domain,
                "hierarchy_path": path,
                "domain_status": domain_status(domain),
                "secondary_structure_status": secondary_structure_status(ss),
                "conservation_class": conservation_class(
                    n_present, n_serotypes_total
                ),
                "n_serotypes_present": n_present,
                "availability_class": availability_class(n_rep, avail, in_all),
                "n_replicates": n_rep,
            }
        )

    out = pd.DataFrame.from_records(records)[list(RESIDUE_ANNOTATION_COLUMNS)]
    return out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
