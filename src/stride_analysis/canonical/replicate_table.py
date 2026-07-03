"""Replicate canonical table construction.

Builds the Level-1 canonical table: one row per (serotype, replicate, residue),
preserving the original per-run quantities. This table is kept **separate** from
the STRIDE summary table — the two data levels are never collapsed.

Unique key: ``(serotype, replicate, canon_label)``.
"""
from __future__ import annotations

import pandas as pd

from ..models import ReplicateInput
from ..models.errors import ConsistencyError


def build_replicate_rows(
    df: pd.DataFrame, rep: ReplicateInput
) -> pd.DataFrame:
    """Attach identity + provenance columns to one replicate's table.

    The original correlation columns are preserved verbatim; identity and
    provenance columns are prepended. ``canon_label`` mirrors the STRIDE
    summary join key (the replicate's ``label`` column).
    """
    out = df.copy()
    # canonical join key: the replicate 'label' == the summary canon_label
    out.insert(0, "serotype", rep.serotype)
    out.insert(1, "replicate", rep.run_dir)
    out.insert(2, "replicate_index", rep.replicate_index)
    out.insert(3, "canon_label", out["label"])
    out["source_path"] = str(rep.correlations_path)
    out["run_dir"] = rep.run_dir
    return out


def assemble_replicate_table(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate per-replicate frames and assert the uniqueness key.

    Union-aligns columns across replicates (milestone drift means some
    replicates may carry extra columns); missing cells become NA. Raises
    :class:`ConsistencyError` if ``(serotype, replicate, canon_label)`` is not
    unique.
    """
    if not frames:
        return pd.DataFrame()
    master = pd.concat(frames, ignore_index=True, sort=False)
    key = ["serotype", "replicate", "canon_label"]
    dup = master.duplicated(subset=key, keep=False)
    if dup.any():
        examples = master.loc[dup, key].drop_duplicates().head(5).to_dict("records")
        raise ConsistencyError(
            f"replicate key {key} not unique; duplicate examples: {examples}"
        )
    lead = ["serotype", "replicate", "replicate_index", "canon_label"]
    rest = [c for c in master.columns if c not in lead]
    master = master[lead + rest]
    return master.sort_values(
        ["serotype", "replicate_index", "canon_label"]
    ).reset_index(drop=True)
