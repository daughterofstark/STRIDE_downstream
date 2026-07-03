"""Task 4 — replicate summaries (availability inventory).

Summarises, per (serotype, canon_label), which replicates observed the residue.
Availability only: replicate count, replicate identifiers, and availability
flags. This module deliberately does **not** average, aggregate effect values,
or compute anything statistical — it only records presence.
"""
from __future__ import annotations

from typing import NamedTuple

import pandas as pd

from ..models.schema import REPLICATE_INVENTORY_COLUMNS


class _Observed(NamedTuple):
    replicates: list[str]
    replicate_indices: list[int]


def build_replicate_inventory(
    replicate_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
) -> pd.DataFrame:
    """Build the replicate availability inventory.

    One row per (serotype, canon_label) present in the canonical residues table.
    Residues with no replicate observations (e.g. a summaries-only S0 run, or a
    residue absent from every replicate) get ``n_replicates == 0`` and
    ``available == False`` — they are recorded, never dropped.
    """
    if canonical_residues.empty:
        return pd.DataFrame(columns=list(REPLICATE_INVENTORY_COLUMNS))

    # total replicates available per serotype (denominator for in_all_replicates)
    total_by_serotype: dict[str, int] = {}
    observed: dict[tuple[str, str], _Observed] = {}
    if not replicate_table.empty:
        total_by_serotype = (
            replicate_table.groupby("serotype")["replicate"].nunique().to_dict()
        )
        grouped = replicate_table.groupby(["serotype", "canon_label"])
        for (serotype, canon_label), grp in grouped:
            observed[(serotype, canon_label)] = _Observed(
                replicates=sorted(grp["replicate"].unique().tolist()),
                replicate_indices=sorted(
                    int(i) for i in grp["replicate_index"].unique()
                ),
            )

    records = []
    for serotype, canon_label in zip(
        canonical_residues["serotype"],
        canonical_residues["canon_label"],
        strict=True,
    ):
        obs = observed.get((serotype, canon_label))
        reps = obs.replicates if obs is not None else []
        idxs = obs.replicate_indices if obs is not None else []
        n = len(reps)
        total = total_by_serotype.get(serotype, 0)
        records.append(
            {
                "serotype": serotype,
                "canon_label": canon_label,
                "n_replicates": n,
                "replicates": reps,
                "replicate_indices": idxs,
                "available": n > 0,
                "in_all_replicates": bool(total > 0 and n == total),
            }
        )

    out = pd.DataFrame.from_records(records)[list(REPLICATE_INVENTORY_COLUMNS)]
    return out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
