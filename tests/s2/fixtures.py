"""Synthetic fixtures for S2 tests.

Builds tiny, S0/S1B-shaped tables directly (no dependence on the real DENV data
and no need to run S0/S1A/S1B). The scenario has two serotypes and a deliberately
varied set of loci so every reduction branch is exercised:

Loci (canon_label → behaviour, ρ ascending across the 7 scales):
- ``NS3:51``  — Catalytic Triad, gates at **residue** at ρ*=0.5 (ρ starts high,
  0.80); a **signed** ``increase`` mechanism whose β CI excludes 0.
- ``NS3:200`` — C-Terminal Tail, gates at **domain** at ρ*=0.5 (ρ climbs only
  when aggregated); a **mixed** mechanism (no signed β).
- ``NS3:99``  — Gly45 Turn, low ρ everywhere so it re-gates coarse / can
  fall **unresolved** at high ρ*; a **signed** ``decrease`` mechanism whose β CI
  **touches** 0 (so it never passes the screen).

The three loci therefore cover: residue vs domain resolution, signed vs mixed,
CI-excludes-0 vs CI-touches-0, and the unresolved sentinel at high ρ*.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# the seven scales, finest → coarsest
_SCALES = [
    (0, "residue"),
    (1, "secondary_structure"),
    (2, "motif"),
    (3, "domain"),
    (4, "chain"),
    (5, "protein"),
    (6, "complex"),
]

# canon_label -> (chain, domain, rho-by-scale-index, direction, beta_signed,
#                 ci_lower, ci_upper, coherence)
_LOCI = {
    "NS3:51": (
        "NS3",
        "Catalytic Triad",
        {0: 0.80, 1: 0.82, 2: 0.85, 3: 0.88, 4: 0.90, 5: 0.92, 6: 0.95},
        "increase",
        0.12,
        0.02,
        0.22,
        0.90,
    ),
    "NS3:200": (
        "NS3",
        "C-Terminal Tail",
        {0: 0.30, 1: 0.35, 2: 0.40, 3: 0.70, 4: 0.75, 5: 0.80, 6: 0.85},
        "mixed",
        None,
        None,
        None,
        0.50,
    ),
    "NS3:99": (
        "NS3",
        "Gly45 Turn",
        {0: 0.20, 1: 0.25, 2: 0.30, 3: 0.55, 4: 0.60, 5: 0.65, 6: 0.70},
        "decrease",
        -0.08,
        -0.15,
        0.00,  # CI touches 0 → never excludes zero
        0.80,
    ),
}
_COMPLEX = "CPLX"
_PROTEIN = "prot"
_SEROTYPES = ("DENVA", "DENVB")

# the gated scale index for each locus at the provisional ρ*=0.5 (finest ρ≥0.5)
_GATED_INDEX = {"NS3:51": 0, "NS3:200": 3, "NS3:99": 3}


def _region_at(chain: str, domain: str, canon: str, scale_index: int) -> str:
    full = [_COMPLEX, _PROTEIN, chain, domain, "none", "unknown", canon]
    depth = 7 - scale_index
    return "/".join(full[:depth])


def make_stride_table() -> pd.DataFrame:
    """A valid S0-shaped STRIDE table (all 7 scales × loci × serotypes)."""
    rows = []
    for serotype in _SEROTYPES:
        for canon, (
            chain,
            domain,
            rho_by_idx,
            direction,
            beta_signed,
            ci_lo,
            ci_hi,
            coherence,
        ) in _LOCI.items():
            gated_idx = _GATED_INDEX[canon]
            for si, sl in _SCALES:
                is_gated = si == gated_idx
                rows.append(
                    {
                        "serotype": serotype,
                        "canon_label": canon,
                        "scale_level": sl,
                        "scale_index": si,
                        "region_id": _region_at(chain, domain, canon, si),
                        "rho": rho_by_idx[si],
                        "gated": is_gated,
                        "beta": 1.0 + si * 0.1,
                        "beta_se": 0.1,
                        "tau2": 0.2,
                        "sigma2_bar": 0.3,
                        "coherence": coherence,
                        "h_chain": chain,
                        "h_domain": domain,
                        "is_gated_scale": is_gated,
                        "mech_direction": direction if is_gated else None,
                        "mech_beta_signed": beta_signed if is_gated else None,
                        "mech_beta_ci_lower": ci_lo if is_gated else None,
                        "mech_beta_ci_upper": ci_hi if is_gated else None,
                        "mech_reproducible_magnitude_energy": (
                            1.0 + si * 0.1 if is_gated else None
                        ),
                    }
                )
    return pd.DataFrame(rows)


def make_residue_annotation() -> pd.DataFrame:
    """A valid S1B-shaped residue annotation for the scenario loci."""
    rows = []
    for serotype in _SEROTYPES:
        for canon, (chain, domain, *_rest) in _LOCI.items():
            dstatus = "assigned"
            rows.append(
                {
                    "serotype": serotype,
                    "canon_label": canon,
                    "chain": chain,
                    "domain": domain,
                    "domain_status": dstatus,
                    "conservation_class": "pan_serotype",
                }
            )
    return pd.DataFrame(rows)


def make_domain_annotation() -> pd.DataFrame:
    """A valid S1B-shaped domain annotation derived from the scenario loci."""
    ra = make_residue_annotation()
    records = []
    for (serotype, chain, domain), grp in ra.groupby(
        ["serotype", "chain", "domain"], sort=True
    ):
        records.append(
            {
                "serotype": serotype,
                "chain": chain,
                "domain": domain,
                "domain_status": "assigned",
                "n_residues": int(len(grp)),
            }
        )
    return pd.DataFrame.from_records(records)


def write_inputs(
    tmpdir: str | Path,
    stride_table: pd.DataFrame | None = None,
    residue_annotation: pd.DataFrame | None = None,
    domain_annotation: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write the three S2 input parquet tables into ``tmpdir``; return paths."""
    d = Path(tmpdir)
    d.mkdir(parents=True, exist_ok=True)
    st = stride_table if stride_table is not None else make_stride_table()
    ra = (
        residue_annotation
        if residue_annotation is not None
        else make_residue_annotation()
    )
    da = (
        domain_annotation
        if domain_annotation is not None
        else make_domain_annotation()
    )
    paths = {
        "stride_table": d / "stride_table.parquet",
        "residue_annotation": d / "residue_annotation.parquet",
        "domain_annotation": d / "domain_annotation.parquet",
    }
    st.to_parquet(paths["stride_table"], index=False)
    ra.to_parquet(paths["residue_annotation"], index=False)
    da.to_parquet(paths["domain_annotation"], index=False)
    return paths


# scenario constants for test expectations
SEROTYPES = _SEROTYPES
N_LOCI_PER_SEROTYPE = len(_LOCI)
