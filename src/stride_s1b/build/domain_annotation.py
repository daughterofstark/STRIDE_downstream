"""Domain annotation.

Builds one structural annotation row per (serotype, chain, domain) from the S1A
domain table plus the residue annotation (for conservation composition). Derived
columns are deterministic: domain assignment status, member residue count, count
of pan-serotype members, and whether the domain is fully conserved.

No statistics, no scoring, no ranking.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    CONSERVATION_PAN,
    DOMAIN_ANNOTATION_COLUMNS,
)
from ._classify import domain_status


def build_domain_annotation(
    domain_table: pd.DataFrame,
    residue_annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Build the domain annotation table.

    Parameters
    ----------
    domain_table
        S1A domain table (one row per (serotype, chain, domain)).
    residue_annotation
        The S1B residue annotation, used to count pan-serotype members per
        domain.

    Returns
    -------
    DataFrame
        One row per (serotype, chain, domain), keyed and column-ordered per
        :data:`~stride_s1b.models.schema.DOMAIN_ANNOTATION_COLUMNS`.
    """
    if domain_table.empty:
        return pd.DataFrame(columns=list(DOMAIN_ANNOTATION_COLUMNS))

    # count pan-serotype residues per (serotype, chain, domain) from residue annot
    pan_by_domain: dict[tuple[str, str, str], int] = {}
    if not residue_annotation.empty:
        pan = residue_annotation[
            residue_annotation["conservation_class"] == CONSERVATION_PAN
        ]
        grouped = pan.groupby(["serotype", "chain", "domain"]).size()
        pan_by_domain = {
            (str(s), str(c), str(d)): int(n)
            for (s, c, d), n in grouped.items()
        }

    records = []
    for serotype, complex_, protein, chain, domain, n_residues, path in zip(
        domain_table["serotype"],
        domain_table["complex"],
        domain_table["protein"],
        domain_table["chain"],
        domain_table["domain"],
        domain_table["n_residues"],
        domain_table["hierarchy_path"],
        strict=True,
    ):
        n_res = int(n_residues)
        n_pan = pan_by_domain.get((serotype, chain, domain), 0)
        records.append(
            {
                "serotype": serotype,
                "complex": complex_,
                "protein": protein,
                "chain": chain,
                "domain": domain,
                "hierarchy_path": path,
                "domain_status": domain_status(domain),
                "n_residues": n_res,
                "n_pan_serotype_residues": n_pan,
                "fully_conserved": bool(n_res > 0 and n_pan == n_res),
            }
        )

    out = pd.DataFrame.from_records(records)[list(DOMAIN_ANNOTATION_COLUMNS)]
    return out.sort_values(["serotype", "chain", "domain"]).reset_index(
        drop=True
    )
