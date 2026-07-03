"""Synthetic fixtures for S4 tests.

Builds a tiny, S0-shaped STRIDE table directly (no dependence on the real DENV
data and no need to run S0). The scenario has two serotypes and a set of loci
across two chains and several domains, chosen so every uncertainty-layer branch
is exercised:

Loci (canon_label → chain / domain, with the variance components and gated
mechanism payload that drive each S4 product):

- ``NS3:51``  — NS3 / Catalytic Triad. Gated at **domain**; a **signed**
  ``increase`` mechanism whose CI ``[0.02, 0.22]`` excludes 0 (z = 2.4).
  Domain variance is **replicate-dominated** (τ²=0.42, σ̄²=0.18 → frac_tau2=0.7).
- ``NS3:200`` — NS3 / C-Terminal Tail. Gated at **domain**; a **signed**
  ``decrease`` mechanism whose CI ``[-0.30, -0.02]`` excludes 0. Domain variance
  is **sampling-dominated** (τ²=0.10, σ̄²=0.40 → frac_tau2=0.2).
- ``NS3:99``  — NS3 / Catalytic Triad (same domain as NS3:51). Gated at
  **residue**; a **signed** ``increase`` mechanism whose CI ``[-0.01, 0.21]``
  **touches 0** (not significant). Gives the Catalytic Triad domain a second
  mechanism so the β_se-weighted mean aggregates >1 signed effect.
- ``NS2B:-1`` — NS2B / Cofactor Interface. Gated at **domain**; a **mixed**
  mechanism (unsigned; no CI). Domain variance is **balanced** (τ²=0.25,
  σ̄²=0.25 → frac_tau2=0.5). Present so a second chain and a mixed/zero-variance
  path are covered.

The four loci therefore cover: replicate- vs sampling-dominated vs balanced
variance regimes; signed increase/decrease vs mixed; CI-excludes-0 vs
CI-touches-0; residue- vs domain-gated tiers; multiple signed mechanisms in one
domain (for the weighted mean); and a genuine two-chain layout.
"""
from __future__ import annotations

from dataclasses import dataclass
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

_COMPLEX = "CPLX"
_PROTEIN = "prot"
_SEROTYPES = ("DENVA", "DENVB")


@dataclass(frozen=True)
class _LocusSpec:
    """A synthetic locus: hierarchy labels, gate scale, mechanism, variances."""

    chain: str
    domain: str
    gated_index: int
    #: mechanism direction on the gated row; ``None`` → mixed (unsigned)
    direction: str | None
    #: (beta_signed, beta_se, ci_lower, ci_upper) on the gated row, or ``None``
    mech: tuple[float, float, float, float] | None
    #: (tau2, sigma2_bar) region-constant values on the domain-scale rows
    dom_var: tuple[float, float]
    #: (tau2, sigma2_bar) values on the residue-scale row
    res_var: tuple[float, float]


# canon_label -> spec
_LOCI: dict[str, _LocusSpec] = {
    "NS3:51": _LocusSpec(
        chain="NS3",
        domain="Catalytic Triad",
        gated_index=3,
        direction="increase",
        mech=(0.12, 0.05, 0.02, 0.22),
        dom_var=(0.42, 0.18),
        res_var=(0.30, 0.30),
    ),
    "NS3:99": _LocusSpec(
        chain="NS3",
        domain="Catalytic Triad",
        gated_index=0,
        direction="increase",
        mech=(0.10, 0.055, -0.01, 0.21),
        dom_var=(0.42, 0.18),  # same domain region as NS3:51
        res_var=(0.20, 0.20),
    ),
    "NS3:200": _LocusSpec(
        chain="NS3",
        domain="C-Terminal Tail",
        gated_index=3,
        direction="decrease",
        mech=(-0.16, 0.05, -0.30, -0.02),
        dom_var=(0.10, 0.40),
        res_var=(0.05, 0.45),
    ),
    "NS2B:-1": _LocusSpec(
        chain="NS2B",
        domain="Cofactor Interface",
        gated_index=3,
        direction=None,  # mixed
        mech=None,
        dom_var=(0.25, 0.25),
        res_var=(0.25, 0.25),
    ),
}


def _region_at(chain: str, domain: str, canon: str, scale_index: int) -> str:
    full = [_COMPLEX, _PROTEIN, chain, domain, "none", "unknown", canon]
    depth = 7 - scale_index
    return "/".join(full[:depth])


def make_stride_table() -> pd.DataFrame:
    """A valid S0-shaped STRIDE table (all 7 scales × loci × serotypes)."""
    rows = []
    for serotype in _SEROTYPES:
        for canon, spec in _LOCI.items():
            chain = spec.chain
            domain = spec.domain
            gated_idx = spec.gated_index
            direction = spec.direction
            mech = spec.mech
            dom_tau2, dom_sigma2 = spec.dom_var
            res_tau2, res_sigma2 = spec.res_var
            for si, sl in _SCALES:
                is_gated = si == gated_idx
                # variance components: residue row uses res_var; domain+ use dom_var
                if si == 0:
                    tau2, sigma2 = res_tau2, res_sigma2
                else:
                    tau2, sigma2 = dom_tau2, dom_sigma2
                row: dict[str, object] = {
                    "serotype": serotype,
                    "canon_label": canon,
                    "scale_level": sl,
                    "scale_index": si,
                    "region_id": _region_at(chain, domain, canon, si),
                    "rho": 0.3 + 0.08 * si,
                    "gated": is_gated,
                    "beta": 1.0 + si * 0.1,
                    "beta_se": 0.2,
                    "tau2": tau2,
                    "sigma2_bar": sigma2,
                    "coherence": 0.9 if direction else 0.4,
                    "h_chain": chain,
                    "h_domain": domain,
                    "is_gated_scale": is_gated,
                    "mech_direction": None,
                    "mech_beta_signed": None,
                    "mech_beta_ci_lower": None,
                    "mech_beta_ci_upper": None,
                    "mech_beta_se": None,
                }
                if is_gated:
                    row["mech_direction"] = direction if direction else "mixed"
                    if mech is not None:
                        bs, bse, lo, hi = mech
                        row["mech_beta_signed"] = bs
                        row["mech_beta_se"] = bse
                        row["mech_beta_ci_lower"] = lo
                        row["mech_beta_ci_upper"] = hi
                rows.append(row)
    return pd.DataFrame(rows)


def write_inputs(
    tmpdir: str | Path,
    stride_table: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write the S4 input parquet table into ``tmpdir``; return the path."""
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
