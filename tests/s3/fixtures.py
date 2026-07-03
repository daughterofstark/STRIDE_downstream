"""Synthetic fixtures for S3 tests.

Builds a tiny, S0-shaped STRIDE table directly (no dependence on the real DENV
data and no need to run S0). The scenario has two serotypes and a deliberately
varied set of loci across **two chains** so every hierarchy-reduction branch is
exercised:

Loci (canon_label → chain / domain → behaviour, ρ ascending across the 7 scales
unless noted):

- ``NS3:51``  — NS3 / Catalytic Triad. High ρ throughout (0.80→0.95), gates at
  **residue** at the provisional ρ*=0.5; small domain−residue gap; monotone;
  a **signed** ``increase`` mechanism.
- ``NS3:200`` — NS3 / C-Terminal Tail. Low-then-high ρ (0.30→0.85), gates at
  **domain**; large domain−residue gap → a **distributed** effect; monotone;
  a **mixed** mechanism.
- ``NS3:99``  — NS3 / Gly45 Turn. **Non-monotone**: ρ dips at the motif scale
  before rising, so the upward-closure audit flags one violation; gates at
  **domain**; a **signed** ``decrease`` mechanism.
- ``NS2B:-1`` — NS2B / Cofactor Interface. Moderate ρ (0.55→0.88), gates at
  **residue**; monotone; a **signed** ``decrease`` mechanism. Present so the
  chain contrast has a second chain (NS2B, negative numbering) to compare.

The four loci therefore cover: residue vs domain resolution, distributed vs
locally-resolved effects, monotone vs non-monotone curves, signed vs mixed
direction, and a genuine two-chain (NS2B vs NS3) contrast.
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

# canon_label -> (chain, domain, rho-by-scale-index, gated_index, direction, beta)
_LOCI = {
    "NS3:51": (
        "NS3",
        "Catalytic Triad",
        {0: 0.80, 1: 0.82, 2: 0.85, 3: 0.88, 4: 0.90, 5: 0.92, 6: 0.95},
        0,
        "increase",
        1.0,
    ),
    "NS3:200": (
        "NS3",
        "C-Terminal Tail",
        {0: 0.30, 1: 0.35, 2: 0.40, 3: 0.70, 4: 0.75, 5: 0.80, 6: 0.85},
        3,
        "mixed",
        1.2,
    ),
    "NS3:99": (
        "NS3",
        "Gly45 Turn",
        # non-monotone: ρ dips at motif (index 2) from 0.45 → 0.40
        {0: 0.42, 1: 0.45, 2: 0.40, 3: 0.72, 4: 0.78, 5: 0.82, 6: 0.90},
        3,
        "decrease",
        1.1,
    ),
    "NS2B:-1": (
        "NS2B",
        "Cofactor Interface",
        {0: 0.55, 1: 0.58, 2: 0.62, 3: 0.70, 4: 0.75, 5: 0.80, 6: 0.88},
        0,
        "decrease",
        0.9,
    ),
}
_COMPLEX = "CPLX"
_PROTEIN = "prot"
_SEROTYPES = ("DENVA", "DENVB")


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
            gated_idx,
            direction,
            beta,
        ) in _LOCI.items():
            for si, sl in _SCALES:
                is_gated = si == gated_idx
                rows.append(
                    {
                        "serotype": serotype,
                        "canon_label": canon,
                        "scale_level": sl,
                        "scale_index": si,
                        "locus": _region_at(chain, domain, canon, 0),
                        "region_id": _region_at(chain, domain, canon, si),
                        "rho": rho_by_idx[si],
                        "gated": is_gated,
                        "beta": beta + si * 0.1,
                        "coherence": 0.9 if direction != "mixed" else 0.5,
                        "h_chain": chain,
                        "h_domain": domain,
                        "is_gated_scale": is_gated,
                        "mech_direction": direction if is_gated else None,
                    }
                )
    return pd.DataFrame(rows)


def write_inputs(
    tmpdir: str | Path,
    stride_table: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write the S3 input parquet table into ``tmpdir``; return the path."""
    d = Path(tmpdir)
    d.mkdir(parents=True, exist_ok=True)
    st = stride_table if stride_table is not None else make_stride_table()
    paths = {"stride_table": d / "stride_table.parquet"}
    st.to_parquet(paths["stride_table"], index=False)
    return paths


# scenario constants for test expectations
SEROTYPES = _SEROTYPES
N_LOCI_PER_SEROTYPE = len(_LOCI)
CHAINS = ("NS2B", "NS3")
