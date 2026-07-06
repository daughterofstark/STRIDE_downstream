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


def replicate_canon_label(df: pd.DataFrame) -> pd.Series:
    """Construct the canonical join key for a replicate table.

    The key must be **identical** to the Level-2 profile's ``canon_label``, which
    STRIDE writes as ``"{chain}:{canon_resid}"`` (e.g. ``"NS3:51"``,
    ``"NS2B:-24"``). The real STRIDE correlations table carries the explicit
    ``chain`` and ``canon_resid`` structural-metadata columns, so the key is
    rebuilt from those — its raw ``label`` column there is a residue-name label
    (e.g. ``"ALA-24"``), *not* the canonical join key.

    When no ``chain`` column is present (e.g. the synthetic example, whose
    ``label`` column already *is* the canonical ``"{chain}:{canon_resid}"`` key),
    the existing ``label`` column is used unchanged — so this preserves the prior
    behaviour on any dataset that lacks explicit ``chain`` metadata.
    """
    label = df["label"].astype(str)
    if "chain" not in df.columns:
        return label
    resid = pd.to_numeric(df["canon_resid"], errors="coerce").astype("Int64")
    chain = df["chain"]
    built = chain.astype(str).str.cat(resid.astype(str), sep=":")
    # only trust the rebuilt key where both components are present; otherwise
    # fall back to the label so partial/absent metadata never yields "nan:…".
    usable = chain.notna() & resid.notna()
    return built.where(usable, label)


def build_replicate_rows(
    df: pd.DataFrame, rep: ReplicateInput
) -> pd.DataFrame:
    """Attach identity + provenance columns to one replicate's table.

    The original correlation columns are preserved verbatim; identity and
    provenance columns are prepended. ``canon_label`` mirrors the STRIDE summary
    join key ``"{chain}:{canon_resid}"`` (see :func:`replicate_canon_label`).
    """
    out = df.copy()
    # canonical join key: identical to the Level-2 profile canon_label
    out.insert(0, "serotype", rep.serotype)
    out.insert(1, "replicate", rep.run_dir)
    out.insert(2, "replicate_index", rep.replicate_index)
    out.insert(3, "canon_label", replicate_canon_label(df))
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
