r"""Internal frame-extraction helpers shared by the S6 builders.

All helpers are pure: they read the loaded replicate inputs and return new frames
without mutating their arguments. ``per_run_effects`` reshapes the Level-1
``replicate_table`` into one finite per-run θ per ``(serotype, canon_label, run)``;
``serotype_effect_matrix`` extracts, for one serotype, the complete-case
position × run matrix used by the concordance computation.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import MIN_REPLICATES_FOR_CONCORDANCE


def chain_from_canon(canon_label: object) -> str:
    """Return the chain prefix of a ``canon_label`` (e.g. ``NS3:51`` → ``NS3``)."""
    text = str(canon_label)
    return text.split(":", 1)[0] if ":" in text else text


def per_run_effects(replicate_table: pd.DataFrame | None) -> pd.DataFrame:
    """Tidy the Level-1 observations to one finite per-run θ per position-run.

    Returns an empty frame with the expected columns when the replicate table is
    absent or carries no finite ``r``. Rows are the finite per-run effects, keyed
    by ``(serotype, canon_label, replicate_index)``; if a run repeats a position
    (it should not) the first finite value in stable order is kept.
    """
    cols = ["serotype", "canon_label", "replicate_index", "r"]
    if replicate_table is None or replicate_table.empty:
        return pd.DataFrame(columns=cols)
    df = replicate_table.loc[:, cols].copy()
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    df = df[df["r"].notna()]
    df["serotype"] = df["serotype"].astype(str)
    df["canon_label"] = df["canon_label"].astype(str)
    df["replicate_index"] = df["replicate_index"].astype(int)
    df = df.sort_values(
        ["serotype", "canon_label", "replicate_index"], kind="stable"
    )
    df = df.drop_duplicates(
        ["serotype", "canon_label", "replicate_index"], keep="first"
    )
    return df.reset_index(drop=True)


def runs_with_effects(effects: pd.DataFrame, serotype: str) -> list[int]:
    """The sorted distinct run indices carrying a finite θ for ``serotype``."""
    if effects.empty:
        return []
    sub = effects[effects["serotype"] == serotype]
    return sorted(int(x) for x in sub["replicate_index"].unique())


def per_replicate_effects_available(effects: pd.DataFrame) -> bool:
    """True iff *any* serotype carries finite θ in ≥ the minimum number of runs."""
    if effects.empty:
        return False
    per_serotype = effects.groupby("serotype")["replicate_index"].nunique()
    return bool((per_serotype >= MIN_REPLICATES_FOR_CONCORDANCE).any())


def serotype_effect_matrix(
    effects: pd.DataFrame, serotype: str
) -> tuple[list[str], list[int], list[list[float]]]:
    """Complete-case position × run θ matrix for one serotype.

    Returns ``(positions, runs, matrix)`` where ``matrix[i][j]`` is the θ of
    ``positions[i]`` in ``runs[j]``; only positions with a finite θ in *every* run
    are kept (complete cases), and both axes are sorted for determinism.
    """
    sub = effects[effects["serotype"] == serotype]
    runs = sorted(int(x) for x in sub["replicate_index"].unique())
    if not runs:
        return [], [], []
    wide = sub.pivot_table(
        index="canon_label",
        columns="replicate_index",
        values="r",
        aggfunc="first",
    )
    wide = wide.reindex(columns=runs)
    complete = wide.dropna(axis=0, how="any")
    complete = complete.sort_index(kind="stable")
    positions = [str(p) for p in complete.index.tolist()]
    matrix = [[float(v) for v in row] for row in complete.to_numpy()]
    return positions, runs, matrix
