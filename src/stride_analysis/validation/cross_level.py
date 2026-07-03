"""Cross-level consistency between replicate observations and STRIDE summaries.

These checks are advisory-strength structural checks that both data levels
describe the *same serotype*: the residue labels a replicate reports should be a
subset of the residue-scale canon_labels the profile knows about. This does not
recompute anything STRIDE did — it only asserts the two levels are talking about
the same residues, so later stages can safely relate them.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError


def check_replicate_summary_alignment(
    replicate_labels: set[str],
    profile_df: pd.DataFrame,
    serotype: str,
    *,
    strict: bool = False,
) -> str:
    """Check that replicate residue labels align with profile residue labels.

    Returns a short human-readable detail string describing the overlap.

    With ``strict=True``, raises :class:`ConsistencyError` if any replicate
    label is absent from the profile's residue labels. With ``strict=False``
    (default) it only reports the overlap — replicate and profile residue sets
    can legitimately differ slightly (e.g. terminal residues), and S0 does not
    force them equal.
    """
    prof_labels = set(
        profile_df.loc[profile_df["scale_level"] == "residue", "canon_label"]
    )
    missing = replicate_labels - prof_labels
    overlap = replicate_labels & prof_labels
    detail = (
        f"{len(overlap)} shared, {len(missing)} replicate-only, "
        f"{len(prof_labels - replicate_labels)} profile-only"
    )
    if strict and missing:
        raise ConsistencyError(
            f"[{serotype}] {len(missing)} replicate residue label(s) absent "
            f"from the profile, e.g. {sorted(missing)[:5]}"
        )
    return detail
