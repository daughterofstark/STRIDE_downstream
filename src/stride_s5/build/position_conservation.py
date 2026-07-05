"""Conservation of reproducibility across shared positions (Tier B — exploratory).

One row per shared ``canon_label`` (every position in the S1A conservation
index): across the serotypes the position is present in, in how many is it
reproducible at the provisional gate (residue-scale ρ ≥ ρ\\*), and how conserved
is that reproducibility (all / majority / some / none). The row also flags
**serotype-divergent** positions — signed and reproducible in some but not all
serotypes — and the **Catalytic Triad** residues, per the design's cross-serotype
conservation and conserved-catalytic-machinery products (§3.3).

Serotype is the unit of replication (n = 4): the per-serotype reproducibility is
aggregated *first* (:func:`~stride_s5.build._frames.position_frame`), then counted
across serotypes — never treating residues as independent samples (§5.2). ρ is
read, never recomputed; the gate is uncalibrated, so every reproducibility
statement is descriptive and stamped with the provisional ρ\\*.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    POSITION_CONSERVATION_COLUMNS,
    PROVISIONAL_RHO_STAR,
    RHO_DECIMALS,
    TIER_EXPLORATORY,
)
from ._classify import conservation_class, is_catalytic_triad
from ._frames import position_frame


def build_position_conservation(
    stride_table: pd.DataFrame,
    conservation_table: pd.DataFrame,
    rho_star: float = PROVISIONAL_RHO_STAR,
) -> pd.DataFrame:
    """Build the per-position cross-serotype conservation table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    conservation_table
        The S1A ``conservation_table`` (the shared-position index).
    rho_star
        The provisional gate threshold (default :data:`PROVISIONAL_RHO_STAR`).

    Returns
    -------
    DataFrame
        One row per ``canon_label``, column-ordered per
        :data:`~stride_s5.models.schema.POSITION_CONSERVATION_COLUMNS`.
    """
    if conservation_table.empty:
        return pd.DataFrame(columns=list(POSITION_CONSERVATION_COLUMNS))

    pf = position_frame(stride_table, rho_star)
    n_total = int(pf["serotype"].nunique()) if not pf.empty else 0

    # index the per-serotype position view by canon_label for fast lookup
    repro_by_label: dict[str, set[str]] = {}
    signed_by_label: dict[str, set[str]] = {}
    rho_by_label: dict[str, list[float]] = {}
    for row in pf.itertuples(index=False):
        canon = str(row.canon_label)
        rho_by_label.setdefault(canon, []).append(float(row.rho_residue))
        if bool(row.reproducible):
            repro_by_label.setdefault(canon, set()).add(str(row.serotype))
            if bool(row.is_signed):
                signed_by_label.setdefault(canon, set()).add(str(row.serotype))

    records = []
    for row in conservation_table.itertuples(index=False):
        canon = str(row.canon_label)
        n_present = int(row.n_serotypes)
        serotypes_present = _as_str_list(row.serotypes_present)
        n_repro = len(repro_by_label.get(canon, set()))
        n_signed_repro = len(signed_by_label.get(canon, set()))
        frac = (n_repro / n_present) if n_present > 0 else float("nan")
        rhos = rho_by_label.get(canon, [])
        divergent = 0 < n_signed_repro < n_total
        records.append(
            {
                "canon_label": canon,
                "chain": str(row.chain),
                "domain": str(row.domain),
                "n_serotypes_total": n_total,
                "n_serotypes_present": n_present,
                "serotypes_present": serotypes_present,
                "in_all_serotypes": bool(row.in_all_serotypes),
                "n_serotypes_reproducible": n_repro,
                "n_serotypes_signed_reproducible": n_signed_repro,
                "frac_reproducible": _round(frac),
                "conservation_class": conservation_class(n_repro, n_present),
                "is_serotype_divergent": bool(divergent),
                "is_catalytic_triad": is_catalytic_triad(canon),
                "rho_residue_min": _round(min(rhos)) if rhos else float("nan"),
                "rho_residue_median": _round(_median(rhos)),
                "rho_residue_max": _round(max(rhos)) if rhos else float("nan"),
                "rho_star": rho_star,
                "is_provisional_rho_star": True,
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(POSITION_CONSERVATION_COLUMNS)
    ]
    out = out.sort_values(["canon_label"]).reset_index(drop=True)
    if out.duplicated(["canon_label"]).any():
        raise ConsistencyError(
            "position conservation has duplicate canon_label rows"
        )
    return out


def _as_str_list(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Iterable):
        return [str(v) for v in value]
    return []


def _median(values: list[float]) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
