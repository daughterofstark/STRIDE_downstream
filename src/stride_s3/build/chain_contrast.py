"""Chain-level contrast: NS2B (cofactor) vs NS3 (protease) and any other chain.

One row per (serotype, chain): aggregate reproducibility and signed-direction
composition over the chain's member loci (design §3.4). Reproducibility is
summarised at both residue scale (Tier B, exploratory) and chain scale (Tier A,
licensed); direction counts come from the gated mechanisms.

Aggregates are taken over the chain's loci with **no region-constant
assumption** — the mean/median summaries are well-defined whether or not every
locus in a chain happens to share an identical coarse-scale ρ.

Pure: no IO, no mutation of inputs. ρ is read from the profile, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    CHAIN_CONTRAST_COLUMNS,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    DIRECTION_MIXED,
    RHO_DECIMALS,
    TIER_LICENSED,
)
from ._frames import (
    gated_slice,
    locus_scale_frame,
    residue_slice,
    scale_slice,
)

_CHAIN_SCALE_LEVEL = "chain"


def build_chain_contrast(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-chain contrast table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, chain), column-ordered per
        :data:`~stride_s3.models.schema.CHAIN_CONTRAST_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(CHAIN_CONTRAST_COLUMNS))

    locus_scale = locus_scale_frame(stride_table)
    residue = residue_slice(locus_scale)
    chain_rows = scale_slice(locus_scale, _CHAIN_SCALE_LEVEL)
    gated = gated_slice(locus_scale)

    # residue-scale aggregates per (serotype, chain)
    res_stats = _residue_stats(residue)
    # chain-scale aggregates per (serotype, chain)
    chain_stats = _chain_stats(chain_rows)
    # gated direction counts per (serotype, chain)
    dir_counts = _direction_counts(gated)

    keys = sorted(set(res_stats) | set(chain_stats) | set(dir_counts))
    records = []
    for serotype, chain in keys:
        n_loci, r_mean, r_median, r_min, r_max, b_mean = res_stats.get(
            (serotype, chain),
            (0, float("nan"), float("nan"), float("nan"), float("nan"), float("nan")),
        )
        c_mean, c_median = chain_stats.get(
            (serotype, chain), (float("nan"), float("nan"))
        )
        n_inc, n_dec, n_mix = dir_counts.get((serotype, chain), (0, 0, 0))
        records.append(
            {
                "serotype": serotype,
                "chain": chain,
                "n_loci": int(n_loci),
                "rho_residue_mean": _round(r_mean),
                "rho_residue_median": _round(r_median),
                "rho_residue_min": _round(r_min),
                "rho_residue_max": _round(r_max),
                "rho_chain_mean": _round(c_mean),
                "rho_chain_median": _round(c_median),
                "beta_residue_mean": _round(b_mean),
                "n_mechanisms": int(n_inc + n_dec + n_mix),
                "n_increase": int(n_inc),
                "n_decrease": int(n_dec),
                "n_mixed": int(n_mix),
                "n_signed": int(n_inc + n_dec),
                "tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[list(CHAIN_CONTRAST_COLUMNS)]
    out = out.sort_values(["serotype", "chain"]).reset_index(drop=True)
    if out.duplicated(["serotype", "chain"]).any():
        raise ConsistencyError(
            "chain contrast has duplicate (serotype, chain) rows"
        )
    return out


def _residue_stats(
    residue: pd.DataFrame,
) -> dict[tuple[str, str], tuple[int, float, float, float, float, float]]:
    stats: dict[tuple[str, str], tuple[int, float, float, float, float, float]] = {}
    if residue.empty:
        return stats
    for (serotype, chain), grp in residue.groupby(
        ["serotype", "h_chain"], sort=True
    ):
        rho = grp["rho"].dropna()
        beta = grp["beta"].dropna()
        stats[(str(serotype), str(chain))] = (
            int(grp["canon_label"].nunique()),
            float(rho.mean()) if not rho.empty else float("nan"),
            float(rho.median()) if not rho.empty else float("nan"),
            float(rho.min()) if not rho.empty else float("nan"),
            float(rho.max()) if not rho.empty else float("nan"),
            float(beta.mean()) if not beta.empty else float("nan"),
        )
    return stats


def _chain_stats(
    chain_rows: pd.DataFrame,
) -> dict[tuple[str, str], tuple[float, float]]:
    stats: dict[tuple[str, str], tuple[float, float]] = {}
    if chain_rows.empty:
        return stats
    for (serotype, chain), grp in chain_rows.groupby(
        ["serotype", "h_chain"], sort=True
    ):
        rho = grp["rho"].dropna()
        stats[(str(serotype), str(chain))] = (
            float(rho.mean()) if not rho.empty else float("nan"),
            float(rho.median()) if not rho.empty else float("nan"),
        )
    return stats


def _direction_counts(
    gated: pd.DataFrame,
) -> dict[tuple[str, str], tuple[int, int, int]]:
    counts: dict[tuple[str, str], tuple[int, int, int]] = {}
    if gated.empty:
        return counts
    for (serotype, chain), grp in gated.groupby(
        ["serotype", "h_chain"], sort=True
    ):
        directions = grp["mech_direction"].astype("string")
        n_inc = int((directions == DIRECTION_INCREASE).sum())
        n_dec = int((directions == DIRECTION_DECREASE).sum())
        n_mix = int((directions == DIRECTION_MIXED).sum())
        counts[(str(serotype), str(chain))] = (n_inc, n_dec, n_mix)
    return counts


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
