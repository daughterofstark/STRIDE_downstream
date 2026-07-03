"""Domain-scale reproducibility (Tier A — licensed).

One row per (serotype, chain, domain): the domain-scale ρ, unsigned reproducible
magnitude, coherence and variance components read from the domain region's
domain-scale profile row, plus a directional-coherence label and the member
residue count (from the S1B domain annotation).

Domain scale is the **licensed** claim level at K=3 (design §0.1, §5.3): every
row is Tier A. A domain region's ρ/β/coherence at the domain scale are
region-level constants shared by all member loci; the builder verifies that
invariant and collapses the shared rows to one per domain.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    COHERENCE_THRESHOLD,
    DOMAIN_REPRODUCIBILITY_COLUMNS,
    TIER_LICENSED,
)
from ._frames import domain_slice


def build_domain_reproducibility(
    stride_table: pd.DataFrame,
    domain_annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Build the domain-scale reproducibility table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    domain_annotation
        The S1B domain annotation (member residue counts + domain status).

    Returns
    -------
    DataFrame
        One row per (serotype, chain, domain), column-ordered per
        :data:`~stride_s2.models.schema.DOMAIN_REPRODUCIBILITY_COLUMNS`.

    Raises
    ------
    ConsistencyError
        If a domain region's domain-scale ρ/β/coherence is not constant across
        the member loci (the region-level invariant the design relies on).
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(DOMAIN_REPRODUCIBILITY_COLUMNS))

    dom = domain_slice(stride_table)
    if dom.empty:
        return pd.DataFrame(columns=list(DOMAIN_REPRODUCIBILITY_COLUMNS))

    # member residue count + domain status keyed by (serotype, chain, domain)
    n_res: dict[tuple[str, str, str], int] = {}
    dstatus_map: dict[tuple[str, str, str], str] = {}
    if not domain_annotation.empty:
        for serotype, chain, domain, dstatus, nres in zip(
            domain_annotation["serotype"],
            domain_annotation["chain"],
            domain_annotation["domain"],
            domain_annotation["domain_status"],
            domain_annotation["n_residues"],
            strict=True,
        ):
            k = (str(serotype), str(chain), str(domain))
            n_res[k] = int(nres)
            dstatus_map[k] = str(dstatus)

    records = []
    group_cols = ["serotype", "h_chain", "h_domain"]
    for keys, grp in dom.groupby(group_cols, sort=True):
        serotype, chain, domain = (str(k) for k in keys)
        _assert_region_constant(grp, serotype, chain, domain, "rho")
        _assert_region_constant(grp, serotype, chain, domain, "coherence")

        first = grp.iloc[0]
        coherence = _f(first["coherence"])
        key = (serotype, chain, domain)
        records.append(
            {
                "serotype": serotype,
                "chain": chain,
                "domain": domain,
                "domain_status": dstatus_map.get(key, ""),
                "region_id": str(first["region_id"]),
                "n_residues": n_res.get(key, int(grp["canon_label"].nunique())),
                "rho_domain": _f(first["rho"]),
                "beta_domain": _f(first["beta"]),
                "coherence_domain": coherence,
                "tau2_domain": _f(first["tau2"]),
                "sigma2_bar_domain": _f(first["sigma2_bar"]),
                "is_coherent": bool(coherence >= COHERENCE_THRESHOLD),
                "tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(DOMAIN_REPRODUCIBILITY_COLUMNS)
    ]
    out = out.sort_values(["serotype", "chain", "domain"]).reset_index(
        drop=True
    )
    if out.duplicated(["serotype", "chain", "domain"]).any():
        raise ConsistencyError(
            "domain reproducibility has duplicate (serotype, chain, domain) rows"
        )
    return out


def _assert_region_constant(
    grp: pd.DataFrame, serotype: str, chain: str, domain: str, column: str
) -> None:
    """Raise if ``column`` is not a single value across the domain's rows."""
    values = grp[column].dropna().unique()
    if len(values) > 1:
        raise ConsistencyError(
            f"domain ({serotype}, {chain}, {domain}) has non-constant "
            f"domain-scale {column} across member loci: {sorted(values)[:5]}"
        )


def _f(value: object) -> float:
    if value is None:
        return float("nan")
    return float(value)  # type: ignore[arg-type]
