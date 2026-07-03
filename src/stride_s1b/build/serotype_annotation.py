"""Serotype annotation.

Builds one dataset-level structural annotation row per serotype by tallying the
residue annotation and domain annotation. Every column is a deterministic count
of category members — no statistics, no scoring, no ranking.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    AVAIL_ALL,
    AVAIL_NONE,
    AVAIL_SOME,
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
    SEROTYPE_ANNOTATION_COLUMNS,
)


def build_serotype_annotation(
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Build the serotype annotation table.

    One row per serotype, tallying category counts from the residue annotation
    and the number of domains from the domain annotation.
    """
    if residue_annotation.empty:
        return pd.DataFrame(columns=list(SEROTYPE_ANNOTATION_COLUMNS))

    n_domains_by_serotype: dict[str, int] = {}
    if not domain_annotation.empty:
        n_domains_by_serotype = {
            str(s): int(n)
            for s, n in domain_annotation.groupby("serotype").size().items()
        }

    records = []
    for serotype, grp in residue_annotation.groupby("serotype", sort=True):
        dom = grp["domain_status"]
        cons = grp["conservation_class"]
        avail = grp["availability_class"]
        records.append(
            {
                "serotype": serotype,
                "n_residues": int(len(grp)),
                "n_domains": n_domains_by_serotype.get(str(serotype), 0),
                "n_assigned_domain_residues": int(
                    (dom == DOMAIN_STATUS_ASSIGNED).sum()
                ),
                "n_unassigned_domain_residues": int(
                    (dom == DOMAIN_STATUS_UNASSIGNED).sum()
                ),
                "n_pan_serotype_residues": int((cons == CONSERVATION_PAN).sum()),
                "n_partial_residues": int((cons == CONSERVATION_PARTIAL).sum()),
                "n_serotype_unique_residues": int(
                    (cons == CONSERVATION_UNIQUE).sum()
                ),
                "n_residues_all_replicates": int((avail == AVAIL_ALL).sum()),
                "n_residues_some_replicates": int((avail == AVAIL_SOME).sum()),
                "n_residues_no_replicates": int((avail == AVAIL_NONE).sum()),
            }
        )

    out = pd.DataFrame.from_records(records)[list(SEROTYPE_ANNOTATION_COLUMNS)]
    return out.sort_values("serotype").reset_index(drop=True)
