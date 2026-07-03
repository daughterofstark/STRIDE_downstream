"""Synthetic fixtures for S1A tests.

Builds tiny, S0-shaped ``stride_table`` and ``replicate_table`` DataFrames
directly (no dependence on the real DENV data, and no need to run S0). The
canonical scenario has three serotypes and a deliberately non-trivial
conservation pattern:

- residue ``NS3:51``  — present in all three serotypes (conserved)
- residue ``NS3:72``  — present in DENVA and DENVB only (partial)
- residue ``NS3:99``  — present in DENVA only (serotype-unique)

so union = 3 residues, intersection = 1. Two domains exercise the domain table.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# hierarchy grammar: complex/protein/chain/domain/motif/secondary_structure/residue
_SCALES = [
    (0, "residue"),
    (1, "secondary_structure"),
    (2, "motif"),
    (3, "domain"),
    (4, "chain"),
    (5, "protein"),
    (6, "complex"),
]

# residue -> hierarchy path segments (7 deep)
_RES_PATHS = {
    "NS3:51": ["CPLX", "protease", "NS3", "Catalytic Triad", "none", "helix", "NS3:51"],
    "NS3:72": ["CPLX", "protease", "NS3", "Catalytic Triad", "none", "sheet", "NS3:72"],
    "NS3:99": ["CPLX", "protease", "NS3", "C-Terminal Tail", "none", "unknown", "NS3:99"],
}

# which serotypes contain which residues (the conservation pattern)
_PRESENCE = {
    "DENVA": ["NS3:51", "NS3:72", "NS3:99"],
    "DENVB": ["NS3:51", "NS3:72"],
    "DENVC": ["NS3:51"],
}

# S0 STRIDE-table columns S1A does not read — filled with harmless placeholders
# so the fixture is shaped like a real stride_table.
_PASSTHROUGH_DEFAULTS: dict[str, object] = {
    "rho": 0.5, "gated": False, "beta": 0.0, "beta_se": 0.1, "tau2": 0.0,
    "sigma2_bar": 0.0, "a_signed": 0.0, "coherence": 0.5, "method": "bayesian",
    "status": "gate_uncertain", "is_gated_scale": False, "mech_label": None,
    "mech_direction": None, "mech_beta_signed": None, "mech_beta_ci_lower": None,
    "mech_beta_ci_upper": None, "mech_beta_se": None, "mech_coherence": None,
    "mech_reproducible_magnitude_energy": None, "mech_rho_star": None,
    "mech_calibrated": None, "mech_gate_uncertain": None, "mech_status": None,
    "mech_region_id": None, "mech_n_loci": None, "profile_source": "x",
    "mechanism_source": "y", "gate_rho_star": 0.5, "gate_alpha": 0.05,
    "gate_coherence_threshold": 0.6, "mechanism_calibrated": False,
    "mechanism_schema_version": "m5",
}


def _region_at(path: list[str], scale_index: int) -> str:
    depth = 7 - scale_index
    return "/".join(path[:depth])


def make_stride_table() -> pd.DataFrame:
    """A valid S0-shaped stride_table with the conservation scenario above."""
    rows = []
    for serotype, residues in _PRESENCE.items():
        for res in residues:
            path = _RES_PATHS[res]
            locus = "/".join(path)
            for si, sl in _SCALES:
                row = {
                    "serotype": serotype,
                    "canon_label": res,
                    "scale_level": sl,
                    "scale_index": si,
                    "locus": locus,
                    "region_id": _region_at(path, si),
                    "h_complex": path[0],
                    "h_protein": path[1],
                    "h_chain": path[2],
                    "h_domain": path[3],
                    "h_motif": path[4],
                    "h_secondary_structure": path[5],
                    "h_residue": path[6],
                    **_PASSTHROUGH_DEFAULTS,
                }
                rows.append(row)
    return pd.DataFrame(rows)


def make_replicate_table(n_replicates: int = 3) -> pd.DataFrame:
    """A valid S0-shaped replicate_table.

    Every residue present in a serotype is observed in every replicate, except
    ``NS3:72`` in ``DENVB`` which is observed in only the first replicate (to
    exercise partial availability / ``in_all_replicates == False``).
    """
    run_names = ["1st_run", "2nd_run", "3rd_run"][:n_replicates]
    rows = []
    for serotype, residues in _PRESENCE.items():
        for res in residues:
            for idx, run in enumerate(run_names, start=1):
                if serotype == "DENVB" and res == "NS3:72" and idx > 1:
                    continue  # partial availability
                rows.append(
                    {
                        "serotype": serotype,
                        "replicate": run,
                        "replicate_index": idx,
                        "canon_label": res,
                        "file_resid": 51,
                        "canon_resid": 51,
                        "name": "HIS",
                        "r": 0.4,
                        "abs_r": 0.4,
                    }
                )
    return pd.DataFrame(rows)


def write_tables(
    tmpdir: str | Path,
    stride: pd.DataFrame | None = None,
    replicate: pd.DataFrame | None = None,
) -> tuple[Path, Path]:
    """Write stride/replicate parquet into ``tmpdir``; return the two paths."""
    d = Path(tmpdir)
    d.mkdir(parents=True, exist_ok=True)
    s = stride if stride is not None else make_stride_table()
    r = replicate if replicate is not None else make_replicate_table()
    sp = d / "stride_table.parquet"
    rp = d / "replicate_table.parquet"
    s.to_parquet(sp, index=False)
    r.to_parquet(rp, index=False)
    return sp, rp
