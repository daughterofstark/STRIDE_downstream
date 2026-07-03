"""Synthetic STRIDE-like dataset generator.

Produces *tiny, schema-valid* Level-1 and Level-2 files in the nested run-dir
layout the framework discovers. Used both to materialise the committed example
dataset and as the basis for test fixtures. Contains no real data and no
biology — the numbers are arbitrary but internally consistent.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from stride_analysis.models.schema import ANALYSIS_OUTPUT_DIRNAME

# two residue "loci" per serotype: A gates at residue, B gates at domain.
_PATHS = {
    "A": ["CPLX", "prot", "NS3", "Catalytic Triad", "none", "unknown", "NS3:51"],
    "B": ["CPLX", "prot", "NS3", "C-Terminal Tail", "none", "unknown", "NS3:200"],
}
_SCALES = [
    (0, "residue"), (1, "secondary_structure"), (2, "motif"), (3, "domain"),
    (4, "chain"), (5, "protein"), (6, "complex"),
]


def _region_at(path: list[str], scale_index: int) -> str:
    depth = 7 - scale_index
    return "/".join(path[:depth])


def make_correlations_df(serotype: str, replicate_index: int) -> pd.DataFrame:
    """A minimal valid Level-1 replicate table (2 residues)."""
    rows = []
    for key, path in _PATHS.items():
        canon = path[-1]
        # per-replicate jitter so replicates differ but stay consistent
        r = (0.4 if key == "A" else -0.2) + 0.01 * replicate_index
        rows.append(
            {
                "file_resid": 51 if key == "A" else 200,
                "canon_resid": 51 if key == "A" else 200,
                "name": "HIS" if key == "A" else "LYS",
                "label": canon,
                "r": r,
                "abs_r": abs(r),
                # a couple of known-optional columns to exercise that path
                "domain_label": "NS3pro",
                "rmsf": 0.8 + 0.05 * replicate_index,
                "theta_se": 0.05,
            }
        )
    return pd.DataFrame(rows)


def make_profile_df(serotype: str) -> pd.DataFrame:
    rho_by_locus = {
        "A": {0: 0.80, 1: 0.82, 2: 0.85, 3: 0.88, 4: 0.90, 5: 0.92, 6: 0.95},
        "B": {0: 0.30, 1: 0.35, 2: 0.40, 3: 0.70, 4: 0.75, 5: 0.80, 6: 0.85},
    }
    gated_scale = {"A": 0, "B": 3}
    rows = []
    for key, path in _PATHS.items():
        locus = "/".join(path)
        canon = path[-1]
        for si, sl in _SCALES:
            rows.append(
                {
                    "protein": serotype,
                    "locus": locus,
                    "canon_label": canon,
                    "scale_index": si,
                    "scale_level": sl,
                    "region_id": _region_at(path, si),
                    "rho": rho_by_locus[key][si],
                    "gated": (si == gated_scale[key]),
                    "beta": 1.0 + si * 0.1,
                    "beta_se": 0.1,
                    "tau2": 0.2,
                    "sigma2_bar": 0.3,
                    "a_signed": (1.0 if key == "A" else -1.0) * (1.0 + si),
                    "coherence": 0.9 if key == "A" else 0.5,
                    "method": "bayesian",
                    "status": "gate_uncertain",
                }
            )
    return pd.DataFrame(rows)


def make_mechanism_dict(serotype: str) -> dict:
    a_path, b_path = _PATHS["A"], _PATHS["B"]
    a_locus, b_locus = "/".join(a_path), "/".join(b_path)
    mech_a = {
        "region_id": a_locus, "label": a_path[-1], "scale_level": "residue",
        "scale_index": 0, "n_loci": 1, "loci": [a_locus], "rho": 0.80,
        "rho_star": 0.5, "calibrated": False, "direction": "increase",
        "beta_signed": 0.12, "beta_ci_lower": 0.02, "beta_ci_upper": 0.22,
        "beta_se": 0.05, "coherence": 0.9,
        "reproducible_magnitude_energy": 1.0, "method": "bayesian",
        "gate_uncertain": True, "status": "gate_uncertain",
    }
    mech_b = {
        "region_id": _region_at(b_path, 3), "label": b_path[3],
        "scale_level": "domain", "scale_index": 3, "n_loci": 1,
        "loci": [b_locus], "rho": 0.70, "rho_star": 0.5, "calibrated": False,
        "direction": "mixed", "beta_signed": None, "beta_ci_lower": None,
        "beta_ci_upper": None, "beta_se": None, "coherence": 0.5,
        "reproducible_magnitude_energy": 1.3, "method": "bayesian",
        "gate_uncertain": True, "status": "gate_uncertain",
    }
    return {
        "schema_version": "m5", "calibrated": False,
        "uncalibrated_note": "provisional threshold; do not make gate claims.",
        "gate": {"rho_star": 0.5, "alpha": 0.05, "coherence_threshold": 0.6},
        "summary": {"n_loci": 2, "n_mechanisms": 2, "n_unresolved": 0,
                    "n_gate_uncertain": 2},
        "mechanisms": [mech_a, mech_b], "unresolved_loci": [],
    }


def write_dataset(
    root: Path,
    serotypes: list[str],
    run_dirs: list[str],
    *,
    with_summaries: bool = True,
) -> Path:
    """Materialise a nested dataset under ``root`` and return it."""
    root = Path(root)
    for idx, run in enumerate(run_dirs, start=1):
        for sero in serotypes:
            analysis = root / run / sero / ANALYSIS_OUTPUT_DIRNAME
            analysis.mkdir(parents=True, exist_ok=True)
            make_correlations_df(sero, idx).to_csv(
                analysis / f"{sero}_correlations_v5.csv", index=False
            )
    if with_summaries:
        summaries = root / "summaries"
        summaries.mkdir(parents=True, exist_ok=True)
        for sero in serotypes:
            make_profile_df(sero).to_csv(
                summaries / f"{sero}_profile.csv", index=False
            )
            (summaries / f"{sero}_mechanism.json").write_text(
                json.dumps(make_mechanism_dict(sero), indent=2)
            )
    return root
