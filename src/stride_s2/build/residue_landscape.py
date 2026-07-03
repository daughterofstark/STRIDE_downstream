"""Residue-scale reproducibility landscape (Tier B — exploratory).

One row per (serotype, canon_label): the residue-scale ρ, unsigned reproducible
magnitude, coherence and variance components, annotated with the residue's
structural labels (chain/domain/conservation, from S1B), plus the finest scale
the locus re-gates at under the provisional ρ*.

Every row is Tier B: residue-scale claims are exploratory, uncalibrated, and
outside the licensed operating range at K=3 (design §0.1, §5.3–5.4). The tier
label makes that explicit on every row so downstream consumers cannot mistake a
residue-scale number for a licensed claim.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    PROVISIONAL_RHO_STAR,
    RESIDUE_LANDSCAPE_COLUMNS,
    RHO_STAR_DECIMALS,
    TIER_EXPLORATORY,
)
from ._frames import residue_slice, rho_by_scale_map
from ._screens import is_residue_scale, regate_scale


def build_residue_landscape(
    stride_table: pd.DataFrame,
    residue_annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Build the residue-scale reproducibility landscape.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    residue_annotation
        The S1B residue annotation (structural labels per residue).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label), column-ordered per
        :data:`~stride_s2.models.schema.RESIDUE_LANDSCAPE_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(RESIDUE_LANDSCAPE_COLUMNS))

    res = residue_slice(stride_table)
    rho_maps = rho_by_scale_map(stride_table)
    provisional = round(PROVISIONAL_RHO_STAR, RHO_STAR_DECIMALS)

    # structural labels keyed by (serotype, canon_label)
    labels: dict[tuple[str, str], tuple[str, str, str, str]] = {}
    if not residue_annotation.empty:
        for serotype, canon_label, chain, domain, dstatus, cclass in zip(
            residue_annotation["serotype"],
            residue_annotation["canon_label"],
            residue_annotation["chain"],
            residue_annotation["domain"],
            residue_annotation["domain_status"],
            residue_annotation["conservation_class"],
            strict=True,
        ):
            labels[(str(serotype), str(canon_label))] = (
                str(chain),
                str(domain),
                str(dstatus),
                str(cclass),
            )

    records = []
    for row in res.itertuples(index=False):
        key = (str(row.serotype), str(row.canon_label))
        chain, domain, dstatus, cclass = labels.get(
            key, (str(row.h_chain), str(row.h_domain), "", "")
        )
        gated_level, gated_idx = regate_scale(rho_maps[key], provisional)
        records.append(
            {
                "serotype": row.serotype,
                "canon_label": row.canon_label,
                "chain": chain,
                "domain": domain,
                "domain_status": dstatus,
                "conservation_class": cclass,
                "rho_residue": _f(row.rho),
                "beta_residue": _f(row.beta),
                "coherence_residue": _f(row.coherence),
                "tau2_residue": _f(row.tau2),
                "sigma2_bar_residue": _f(row.sigma2_bar),
                "gated_scale_level_provisional": gated_level,
                "gated_scale_index_provisional": gated_idx,
                "gates_at_residue_provisional": is_residue_scale(gated_level),
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame.from_records(records)[list(RESIDUE_LANDSCAPE_COLUMNS)]
    out = out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)

    if out.duplicated(["serotype", "canon_label"]).any():
        raise ConsistencyError(
            "residue landscape has duplicate (serotype, canon_label) rows"
        )
    return out


def _f(value: object) -> float:
    """Coerce a numeric cell to float, preserving NaN."""
    if value is None:
        return float("nan")
    return float(value)  # type: ignore[arg-type]
