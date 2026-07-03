"""Synthetic fixtures for S1B tests.

Builds tiny, S1A-shaped tables directly (no dependence on the real DENV data and
no need to run S0/S1A). The scenario has three serotypes and a deliberately
varied set of annotations so every category is exercised:

Residues (canon_label → structural facts):
- ``NS3:51``  — domain "Catalytic Triad" (assigned), SS "helix" (resolved),
  present in all three serotypes (pan_serotype).
- ``NS3:72``  — domain "Catalytic Triad" (assigned), SS "unknown" (unresolved),
  present in DENVA and DENVB (partial).
- ``NS3:99``  — domain "unassigned" (unassigned), SS "sheet" (resolved),
  present in DENVA only (serotype_unique).

Replicate availability:
- most residues observed in all replicates;
- ``NS3:72`` in DENVB observed in only one replicate (some_replicates).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# residue -> (chain, domain, motif, secondary_structure)
_RES_ANNOT = {
    "NS3:51": ("NS3", "Catalytic Triad", "none", "helix"),
    "NS3:72": ("NS3", "Catalytic Triad", "none", "unknown"),
    "NS3:99": ("NS3", "unassigned", "none", "sheet"),
}
_COMPLEX = "CPLX"
_PROTEIN = "protease"

# which serotypes contain which residues (the conservation pattern)
_PRESENCE = {
    "DENVA": ["NS3:51", "NS3:72", "NS3:99"],
    "DENVB": ["NS3:51", "NS3:72"],
    "DENVC": ["NS3:51"],
}
_N_SEROTYPES = len(_PRESENCE)


def _hierarchy_path(res: str) -> str:
    chain, domain, motif, ss = _RES_ANNOT[res]
    return "/".join([_COMPLEX, _PROTEIN, chain, domain, motif, ss, res])


def make_canonical_residues() -> pd.DataFrame:
    """A valid S1A-shaped canonical residues table."""
    rows = []
    for serotype, residues in _PRESENCE.items():
        for res in residues:
            chain, domain, motif, ss = _RES_ANNOT[res]
            rows.append(
                {
                    "serotype": serotype,
                    "canon_label": res,
                    "chain": chain,
                    "domain": domain,
                    "motif": motif,
                    "secondary_structure": ss,
                    "hierarchy_path": _hierarchy_path(res),
                    "complex": _COMPLEX,
                    "protein": _PROTEIN,
                    "residue": res,
                }
            )
    return pd.DataFrame(rows)


def make_domain_table() -> pd.DataFrame:
    """A valid S1A-shaped domain table derived from the canonical residues."""
    cr = make_canonical_residues()
    records = []
    group_cols = ["serotype", "complex", "protein", "chain", "domain"]
    for keys, grp in cr.groupby(group_cols, sort=True):
        serotype, complex_, protein, chain, domain = keys
        labels = sorted(grp["canon_label"].tolist())
        records.append(
            {
                "serotype": serotype,
                "complex": complex_,
                "protein": protein,
                "chain": chain,
                "domain": domain,
                "n_residues": len(labels),
                "canon_labels": labels,
                "hierarchy_path": f"{complex_}/{protein}/{chain}/{domain}",
            }
        )
    return pd.DataFrame.from_records(records)


def make_conservation_table() -> pd.DataFrame:
    """A valid S1A-shaped conservation table."""
    all_serotypes = sorted(_PRESENCE)
    present_by_label: dict[str, set[str]] = {}
    for serotype, residues in _PRESENCE.items():
        for res in residues:
            present_by_label.setdefault(res, set()).add(serotype)
    records = []
    for res in sorted(present_by_label):
        present = sorted(present_by_label[res])
        absent = sorted(set(all_serotypes) - set(present))
        chain, domain, _, _ = _RES_ANNOT[res]
        records.append(
            {
                "canon_label": res,
                "n_serotypes": len(present),
                "serotypes_present": present,
                "serotypes_absent": absent,
                "in_all_serotypes": len(present) == _N_SEROTYPES,
                "in_any_serotype": True,
                "chain": chain,
                "domain": domain,
            }
        )
    return pd.DataFrame.from_records(records)


def make_replicate_inventory(n_replicates: int = 3) -> pd.DataFrame:
    """A valid S1A-shaped replicate inventory.

    Every residue is observed in all replicates except ``NS3:72`` in DENVB,
    which is observed in only one replicate (``some_replicates``).
    """
    run_names = ["1st_run", "2nd_run", "3rd_run"][:n_replicates]
    rows = []
    for serotype, residues in _PRESENCE.items():
        for res in residues:
            if serotype == "DENVB" and res == "NS3:72":
                reps = run_names[:1]
            else:
                reps = list(run_names)
            n = len(reps)
            rows.append(
                {
                    "serotype": serotype,
                    "canon_label": res,
                    "n_replicates": n,
                    "replicates": reps,
                    "replicate_indices": list(range(1, n + 1)),
                    "available": n > 0,
                    "in_all_replicates": n == n_replicates,
                }
            )
    return pd.DataFrame.from_records(rows)


def make_empty_replicate_inventory() -> pd.DataFrame:
    """A replicate inventory where every residue has zero replicates."""
    cr = make_canonical_residues()
    rows = []
    for serotype, res in zip(cr["serotype"], cr["canon_label"], strict=True):
        rows.append(
            {
                "serotype": serotype,
                "canon_label": res,
                "n_replicates": 0,
                "replicates": [],
                "replicate_indices": [],
                "available": False,
                "in_all_replicates": False,
            }
        )
    return pd.DataFrame.from_records(rows)


def write_s1a_tables(
    tmpdir: str | Path,
    canonical_residues: pd.DataFrame | None = None,
    domain_table: pd.DataFrame | None = None,
    replicate_inventory: pd.DataFrame | None = None,
    conservation_table: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write the four S1A parquet tables into ``tmpdir``; return their paths."""
    d = Path(tmpdir)
    d.mkdir(parents=True, exist_ok=True)
    cr = canonical_residues if canonical_residues is not None else make_canonical_residues()
    dt = domain_table if domain_table is not None else make_domain_table()
    ri = (
        replicate_inventory
        if replicate_inventory is not None
        else make_replicate_inventory()
    )
    ct = conservation_table if conservation_table is not None else make_conservation_table()
    paths = {
        "canonical_residues": d / "canonical_residues.parquet",
        "domain_table": d / "domain_table.parquet",
        "replicate_inventory": d / "replicate_inventory.parquet",
        "conservation_table": d / "conservation_table.parquet",
    }
    cr.to_parquet(paths["canonical_residues"], index=False)
    dt.to_parquet(paths["domain_table"], index=False)
    ri.to_parquet(paths["replicate_inventory"], index=False)
    ct.to_parquet(paths["conservation_table"], index=False)
    return paths


# number of serotypes in the scenario (for conservation-class expectations)
N_SEROTYPES_TOTAL = _N_SEROTYPES
