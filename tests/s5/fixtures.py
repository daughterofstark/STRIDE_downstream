"""Synthetic fixtures for S5 tests.

Builds a tiny, S0-shaped STRIDE table and an S1A-shaped conservation table
directly (no dependence on the real DENV data and no need to run S0/S1A). The
scenario has the **four** dengue serotypes and a set of shared canonical
positions across two chains and four domains, chosen so every cross-serotype
branch is exercised.

Positions (canon_label → chain / domain, with per-serotype residue-scale ρ and
the residue-gated mechanism direction; a serotype omitted from a locus is
*absent* there):

- ``NS3:51``  — NS3 / Catalytic Triad (catalytic triad). Reproducible + signed
  ``increase`` in **all four** serotypes → conservation ``reproducible_all``,
  concordance ``agree`` (increase), not divergent.
- ``NS3:75``  — NS3 / Catalytic Triad (catalytic triad). Reproducible + signed in
  all four but split 2 ``increase`` / 2 ``decrease`` → concordance ``conflict``.
- ``NS3:135`` — NS3 / Catalytic Triad (catalytic triad). Reproducible + signed in
  all four, 3 ``increase`` / 1 ``decrease`` → concordance ``majority``.
- ``NS3:200`` — NS3 / C-Terminal Tail (non-catalytic). ρ = 0.3 in every serotype →
  reproducible in **none** → conservation ``reproducible_none``; absent from the
  concordance table (no signed serotype).
- ``NS2B:-1`` — NS2B / Cofactor Interface (non-catalytic). ``increase`` in DENV1,
  reproducible-but-``mixed`` in DENV2, not reproducible in DENV3, **absent** in
  DENV4 → present in 3, reproducible in 2, signed in 1 → conservation
  ``reproducible_majority``, **divergent**; excluded from concordance (one signed
  serotype only).
- ``NS3:250`` — NS3 / Oxyanion Loop (catalytic domain). Signed ``decrease`` in
  DENV1/DENV2, not reproducible in DENV3/DENV4 → conservation ``reproducible_some``,
  **divergent**, concordance ``agree`` (decrease).

This covers: all four conservation classes; agree / majority / conflict
concordance and the ``< 2`` signed exclusion; divergent True/False; the Catalytic
Triad residue flag and the catalytic-domain flag (Catalytic Triad + Oxyanion
Loop); a reproducible-but-mixed position; an absent-in-one-serotype position; and
a genuine two-chain / four-domain layout.
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
_DOMAIN_INDEX = 3
_RESIDUE_INDEX = 0

_COMPLEX = "CPLX"
_PROTEIN = "prot"
SEROTYPES = ("DENV1", "DENV2", "DENV3", "DENV4")

# canon_label -> (chain, domain, {serotype: (rho_residue, direction | None)})
# direction is the residue-gated mechanism direction; used only when ρ ≥ 0.5.
# A serotype absent from the inner dict means the position is absent there.
_LOCI: dict[str, tuple[str, str, dict[str, tuple[float, str | None]]]] = {
    "NS3:51": (
        "NS3",
        "Catalytic Triad",
        {
            "DENV1": (0.80, "increase"),
            "DENV2": (0.82, "increase"),
            "DENV3": (0.78, "increase"),
            "DENV4": (0.85, "increase"),
        },
    ),
    "NS3:75": (
        "NS3",
        "Catalytic Triad",
        {
            "DENV1": (0.80, "increase"),
            "DENV2": (0.76, "increase"),
            "DENV3": (0.70, "decrease"),
            "DENV4": (0.72, "decrease"),
        },
    ),
    "NS3:135": (
        "NS3",
        "Catalytic Triad",
        {
            "DENV1": (0.80, "increase"),
            "DENV2": (0.70, "increase"),
            "DENV3": (0.62, "increase"),
            "DENV4": (0.60, "decrease"),
        },
    ),
    "NS3:200": (
        "NS3",
        "C-Terminal Tail",
        {
            "DENV1": (0.30, None),
            "DENV2": (0.30, None),
            "DENV3": (0.30, None),
            "DENV4": (0.30, None),
        },
    ),
    "NS2B:-1": (
        "NS2B",
        "Cofactor Interface",
        {
            "DENV1": (0.80, "increase"),
            "DENV2": (0.60, "mixed"),
            "DENV3": (0.30, None),
            # absent in DENV4
        },
    ),
    "NS3:250": (
        "NS3",
        "Oxyanion Loop",
        {
            "DENV1": (0.70, "decrease"),
            "DENV2": (0.60, "decrease"),
            "DENV3": (0.30, None),
            "DENV4": (0.30, None),
        },
    ),
}

# region-constant domain-scale values per domain (identical across member loci)
_DOMAIN_VALUES: dict[str, tuple[float, float, float, float, float]] = {
    # domain -> (rho_domain, beta, beta_se, tau2, sigma2_bar)
    "Catalytic Triad": (0.88, 1.30, 0.10, 0.20, 0.30),
    "C-Terminal Tail": (0.55, 1.10, 0.12, 0.25, 0.35),
    "Cofactor Interface": (0.60, 1.00, 0.15, 0.25, 0.25),
    "Oxyanion Loop": (0.66, 1.20, 0.11, 0.22, 0.28),
}

_RHO_STAR = 0.5


def _region_at(chain: str, domain: str, canon: str, scale_index: int) -> str:
    full = [_COMPLEX, _PROTEIN, chain, domain, "none", "unknown", canon]
    depth = 7 - scale_index
    return "/".join(full[:depth])


def make_stride_table() -> pd.DataFrame:
    """A valid S0-shaped STRIDE table (all 7 scales × present loci × serotypes)."""
    rows = []
    for canon, (chain, domain, per_serotype) in _LOCI.items():
        dom_rho, dom_beta, dom_bse, dom_tau2, dom_sig2 = _DOMAIN_VALUES[domain]
        for serotype, (rho_residue, direction) in per_serotype.items():
            reproducible = rho_residue >= _RHO_STAR
            gated_index = _RESIDUE_INDEX if reproducible else _DOMAIN_INDEX
            for si, sl in _SCALES:
                is_gated = si == gated_index
                if si == _RESIDUE_INDEX:
                    rho = rho_residue
                    beta, bse, tau2, sig2 = 1.0, 0.2, 0.15, 0.15
                elif si == _DOMAIN_INDEX:
                    rho = dom_rho
                    beta, bse, tau2, sig2 = dom_beta, dom_bse, dom_tau2, dom_sig2
                else:
                    rho = 0.4 + 0.05 * si
                    beta, bse, tau2, sig2 = 1.0 + si * 0.1, 0.2, 0.2, 0.2
                mech_direction: str | None = None
                if is_gated:
                    if reproducible:
                        mech_direction = direction if direction else "mixed"
                    else:
                        # coarse-gated: irrelevant to residue-scale signed claims
                        mech_direction = "mixed"
                rows.append(
                    {
                        "serotype": serotype,
                        "canon_label": canon,
                        "scale_level": sl,
                        "scale_index": si,
                        "region_id": _region_at(chain, domain, canon, si),
                        "rho": rho,
                        "gated": is_gated,
                        "beta": beta,
                        "beta_se": bse,
                        "tau2": tau2,
                        "sigma2_bar": sig2,
                        "h_chain": chain,
                        "h_domain": domain,
                        "is_gated_scale": is_gated,
                        "mech_direction": mech_direction,
                    }
                )
    return pd.DataFrame(rows)


def make_conservation_table() -> pd.DataFrame:
    """An S1A-shaped conservation table consistent with :func:`make_stride_table`."""
    all_serotypes = set(SEROTYPES)
    records = []
    for canon, (chain, domain, per_serotype) in _LOCI.items():
        present = sorted(per_serotype.keys())
        absent = sorted(all_serotypes - set(present))
        records.append(
            {
                "canon_label": canon,
                "n_serotypes": len(present),
                "serotypes_present": present,
                "serotypes_absent": absent,
                "in_all_serotypes": len(present) == len(SEROTYPES),
                "in_any_serotype": True,
                "chain": chain,
                "domain": domain,
            }
        )
    # sort by canon_label for determinism, mirroring S1A
    records.sort(key=lambda r: str(r["canon_label"]))
    return pd.DataFrame(records)


def write_inputs(
    tmpdir: str | Path,
    stride_table: pd.DataFrame | None = None,
    conservation_table: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write the S5 input parquet tables into ``tmpdir``; return the paths."""
    d = Path(tmpdir)
    d.mkdir(parents=True, exist_ok=True)
    st = stride_table if stride_table is not None else make_stride_table()
    ct = (
        conservation_table
        if conservation_table is not None
        else make_conservation_table()
    )
    paths = {
        "stride_table": d / "stride_table.parquet",
        "conservation_table": d / "conservation_table.parquet",
    }
    st.to_parquet(paths["stride_table"], index=False)
    ct.to_parquet(paths["conservation_table"], index=False)
    return paths


# scenario constants for test expectations
N_SEROTYPES = len(SEROTYPES)
N_POSITIONS = len(_LOCI)
CHAINS = ("NS2B", "NS3")
CATALYTIC_DOMAIN_CELLS = 8  # Catalytic Triad×4 + Oxyanion Loop×4
