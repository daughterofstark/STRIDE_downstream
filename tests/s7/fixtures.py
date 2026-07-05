"""Synthetic fixtures for S7 tests.

Builds tiny, S2–S6-shaped output frames directly — no dependence on the real DENV
data, on running S0–S6, or on any earlier stage's tests. The default scenario has
the four dengue serotypes over six canonical positions and five domains (two of
them catalytic: ``Catalytic Triad`` and ``Oxyanion Loop``), rich enough that every
figure/table has multiple rows, the catalytic filters are non-empty, and the T4
ranking is non-degenerate.

``write_inputs`` lays the twelve tables out across ``s2/ s3/ s4/ s5/ s6``
subdirectories exactly as the prior-stage CLIs would. ``write_empty_inputs`` writes
the same twelve tables **present but empty** (correct columns, zero rows) to
exercise S7's empty-input handling.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from stride_s7.io.loaders import INPUT_REQUIRED_COLUMNS, STAGE_INPUTS
from stride_s7.models.schema import (
    IN_DOMAIN_EFFECT_SUMMARY,
    IN_DOMAIN_REPRODUCIBILITY,
    IN_DOMAIN_SEROTYPE_MATRIX,
    IN_POSITION_CONSERVATION,
    IN_REPLICATE_BLOCKED_ANALYSES,
    IN_REPLICATE_REGIME,
    IN_RESIDUE_LANDSCAPE,
    IN_RESOLUTION_CENSUS,
    IN_SCALE_CURVE,
    IN_SEROTYPE_SUMMARY,
    IN_SIGNIFICANCE_SCREEN,
    IN_VARIANCE_BUDGET,
    PROVISIONAL_RHO_STAR,
)

SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]

# (canon_label, chain, domain, is_catalytic_domain)
POSITIONS = [
    ("NS3:51", "NS3", "Catalytic Triad", True),
    ("NS3:75", "NS3", "Catalytic Triad", True),
    ("NS3:135", "NS3", "Catalytic Triad", True),
    ("NS3:38", "NS3", "Oxyanion Loop", True),
    ("NS3:200", "NS3", "Helicase P-loop", False),
    ("NS2B:1", "NS2B", "Cofactor", False),
]
# (domain, chain, is_catalytic_domain)
DOMAINS = [
    ("Catalytic Triad", "NS3", True),
    ("Oxyanion Loop", "NS3", True),
    ("Helicase P-loop", "NS3", False),
    ("Cofactor", "NS2B", False),
    ("RNA Groove", "NS3", False),
]
SCALE_LEVELS = [
    "residue",
    "secondary_structure",
    "motif",
    "domain",
    "chain",
    "protein",
    "complex",
]


def _rho(seed: int) -> float:
    return round(0.5 + (seed % 5) * 0.09, 6)


def make_residue_landscape() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for pi, (canon, chain, domain, _cat) in enumerate(POSITIONS):
            rows.append(
                {
                    "serotype": s,
                    "canon_label": canon,
                    "chain": chain,
                    "domain": domain,
                    "domain_status": "named",
                    "conservation_class": "reproducible_all",
                    "rho_residue": _rho(si + pi),
                    "coherence_residue": round(0.6 + 0.05 * pi, 6),
                    "tier": "exploratory",
                }
            )
    return pd.DataFrame(rows)


def make_resolution_census() -> pd.DataFrame:
    rows = []
    levels = [("residue", 0, 4), ("secondary_structure", 1, 1), ("chain", 4, 1)]
    for s in SEROTYPES:
        for prov in (True, False):
            rho_star = PROVISIONAL_RHO_STAR if prov else 0.7
            for level, idx, base in levels:
                rows.append(
                    {
                        "serotype": s,
                        "rho_star": rho_star,
                        "is_provisional_rho_star": prov,
                        "gated_scale_level": level,
                        "gated_scale_index": idx,
                        "tier": "licensed" if idx >= 3 else "exploratory",
                        "n_loci": base + (1 if prov else 0),
                    }
                )
    return pd.DataFrame(rows)


def make_serotype_summary() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for prov in (True, False):
            rho_star = PROVISIONAL_RHO_STAR if prov else 0.7
            rows.append(
                {
                    "serotype": s,
                    "rho_star": rho_star,
                    "is_provisional_rho_star": prov,
                    "n_loci": 200 + si,
                    "n_gated_residue": 120 + si,
                    "n_gated_domain_or_coarser": 8,
                    "n_unresolved": 0,
                    "frac_gated_residue": round(0.6 + 0.01 * si, 6),
                    "n_mechanisms": 129 + si,
                    "n_signed": 30 + si,
                    "n_mixed": 60 + si,
                    "n_signed_significant": 20 + si,
                    "frac_mixed": round(0.45 + 0.01 * si, 6),
                    "rho_residue_median": round(0.70 + 0.01 * si, 6),
                    "rho_residue_q1": round(0.60 + 0.01 * si, 6),
                    "rho_residue_q3": round(0.80 + 0.01 * si, 6),
                    "rho_residue_min": 0.5,
                    "rho_residue_max": round(0.95 + 0.005 * si, 6),
                }
            )
    return pd.DataFrame(rows)


def make_domain_reproducibility() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for di, (domain, chain, _cat) in enumerate(DOMAINS):
            rows.append(
                {
                    "serotype": s,
                    "chain": chain,
                    "domain": domain,
                    "domain_status": "named",
                    "region_id": f"NS2B-NS3/protease/{chain}/{domain}",
                    "n_residues": 3 + di,
                    "rho_domain": _rho(si + di + 1),
                    "beta_domain": round(1.0 + 0.1 * di, 6),
                    "coherence_domain": round(0.55 + 0.06 * di, 6),
                    "tau2_domain": round(0.02 + 0.01 * di, 6),
                    "sigma2_bar_domain": round(0.05 + 0.01 * di, 6),
                    "is_coherent": bool(di % 2 == 0),
                    "tier": "licensed",
                }
            )
    return pd.DataFrame(rows)


def make_scale_curve() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for canon, chain, domain, _cat in POSITIONS:
            for idx, level in enumerate(SCALE_LEVELS):
                rows.append(
                    {
                        "serotype": s,
                        "canon_label": canon,
                        "chain": chain,
                        "domain": domain,
                        "scale_index": idx,
                        "scale_level": level,
                        "rho": round(min(0.99, 0.5 + 0.06 * idx + 0.01 * si), 6),
                        "rho_prev": None,
                        "rho_step_gain": 0.0,
                        "rho_cumulative_gain": 0.0,
                        "is_gated_scale": idx == 0,
                        "tier": "exploratory" if idx == 0 else "licensed",
                    }
                )
    return pd.DataFrame(rows)


def make_significance_screen() -> pd.DataFrame:
    rows = []
    directions = ["increase", "decrease", "mixed"]
    for si, s in enumerate(SEROTYPES):
        for pi, (canon, chain, domain, _cat) in enumerate(POSITIONS):
            direction = directions[(si + pi) % 3]
            is_signed = direction != "mixed"
            beta = round(0.4 - 0.1 * pi + 0.05 * si, 6) if is_signed else float("nan")
            rows.append(
                {
                    "serotype": s,
                    "canon_label": canon,
                    "chain": chain,
                    "domain": domain,
                    "gated_scale_level": "residue",
                    "tier": "exploratory",
                    "direction": direction,
                    "is_signed": is_signed,
                    "beta_signed": beta,
                    "beta_se": 0.05 if is_signed else float("nan"),
                    "beta_ci_lower": round(beta - 0.1, 6) if is_signed else float("nan"),
                    "beta_ci_upper": round(beta + 0.1, 6) if is_signed else float("nan"),
                    "ci_excludes_zero": is_signed and pi < 3,
                    "z_score": 2.0 if is_signed else float("nan"),
                    "p_value": 0.01 if is_signed else float("nan"),
                    "p_value_bh": 0.02 if is_signed else float("nan"),
                    "significant_raw": is_signed and pi < 3,
                    "significant_fdr": is_signed and pi < 2,
                }
            )
    return pd.DataFrame(rows)


def make_variance_budget() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for di, (domain, chain, _cat) in enumerate(DOMAINS):
            tau2 = round(0.02 + 0.01 * di, 6)
            sigma2 = round(0.08 - 0.005 * di, 6)
            total = tau2 + sigma2
            rows.append(
                {
                    "serotype": s,
                    "chain": chain,
                    "domain": domain,
                    "region_id": f"NS2B-NS3/protease/{chain}/{domain}",
                    "rho_domain": _rho(si + di + 1),
                    "beta_domain": round(1.0 + 0.1 * di, 6),
                    "beta_se_domain": 0.05,
                    "tau2": tau2,
                    "sigma2_bar": sigma2,
                    "total_unreproduced": round(total, 6),
                    "frac_tau2": round(tau2 / total, 6),
                    "frac_sigma2": round(sigma2 / total, 6),
                    "tau2_sigma2_ratio": round(tau2 / sigma2, 6),
                    "variance_regime": "sampling_dominated" if di % 2 else "replicate_dominated",
                    "tier": "licensed",
                }
            )
    return pd.DataFrame(rows)


def make_domain_effect_summary() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for di, (domain, chain, _cat) in enumerate(DOMAINS):
            rows.append(
                {
                    "serotype": s,
                    "chain": chain,
                    "domain": domain,
                    "n_mechanisms": 3 + di,
                    "n_signed": 2 + (di % 2),
                    "n_mixed": 1,
                    "n_ci_excludes_zero": 1 + (di % 2),
                    "frac_ci_excludes_zero": round(0.5 + 0.05 * di, 6),
                    "n_significant_fdr": di % 2,
                    "beta_weighted_mean": round(0.3 - 0.05 * di + 0.02 * si, 6),
                    "beta_weighted_se": round(0.04 + 0.005 * di, 6),
                    "beta_unweighted_mean": round(0.28 - 0.05 * di, 6),
                    "tier": "licensed",
                }
            )
    return pd.DataFrame(rows)


def make_domain_serotype_matrix() -> pd.DataFrame:
    rows = []
    for si, s in enumerate(SEROTYPES):
        for di, (domain, chain, cat) in enumerate(DOMAINS):
            rows.append(
                {
                    "serotype": s,
                    "chain": chain,
                    "domain": domain,
                    "region_id": f"NS2B-NS3/protease/{chain}/{domain}",
                    "rho_domain": _rho(si + di + 1),
                    "beta_domain": round(1.0 + 0.1 * di, 6),
                    "beta_se_domain": 0.05,
                    "tau2_domain": round(0.02 + 0.01 * di, 6),
                    "sigma2_bar_domain": round(0.05 + 0.01 * di, 6),
                    "is_catalytic_domain": cat,
                    "tier": "licensed",
                }
            )
    return pd.DataFrame(rows)


def make_position_conservation() -> pd.DataFrame:
    classes = [
        "reproducible_all",
        "reproducible_majority",
        "reproducible_some",
        "reproducible_none",
    ]
    rows = []
    for pi, (canon, chain, domain, _cat) in enumerate(POSITIONS):
        n_repro = 4 - (pi % 4)
        rows.append(
            {
                "canon_label": canon,
                "chain": chain,
                "domain": domain,
                "n_serotypes_total": 4,
                "n_serotypes_present": 4,
                "serotypes_present": "DENV1;DENV2;DENV3;DENV4",
                "in_all_serotypes": True,
                "n_serotypes_reproducible": n_repro,
                "n_serotypes_signed_reproducible": max(0, n_repro - 1),
                "frac_reproducible": round(n_repro / 4.0, 6),
                "conservation_class": classes[pi % 4],
                "is_serotype_divergent": bool(pi % 2),
                "is_catalytic_triad": domain == "Catalytic Triad",
                "rho_residue_min": 0.5,
                "rho_residue_median": round(0.7 + 0.02 * pi, 6),
                "rho_residue_max": round(0.9 + 0.005 * pi, 6),
                "rho_star": PROVISIONAL_RHO_STAR,
                "is_provisional_rho_star": True,
                "tier": "exploratory",
            }
        )
    return pd.DataFrame(rows)


def make_replicate_regime() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        rows.append(
            {
                "serotype": s,
                "n_replicates": 3,
                "n_positions": 6,
                "n_positions_in_all_replicates": 6,
                "frac_complete": 1.0,
                "residue_claims_licensed": False,
                "per_replicate_effects_available": False,
                "n_replicates_with_effects": 0,
            }
        )
    return pd.DataFrame(rows)


def make_replicate_blocked_analyses() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "analysis_id": "per_run_rank_concordance",
                "description": "per-run rank concordance across replicates",
                "status": "blocked",
                "available": False,
                "reason": "per-run correlation CSVs not supplied",
                "required_input": "replicate_table",
                "design_ref": "3.1",
            },
            {
                "analysis_id": "leave_one_replicate_out",
                "description": "LOO-replicate stability",
                "status": "blocked",
                "available": False,
                "reason": "needs a STRIDE re-run",
                "required_input": "per-run CSVs",
                "design_ref": "4.2",
            },
        ]
    )


_BUILDERS = {
    IN_RESIDUE_LANDSCAPE: make_residue_landscape,
    IN_RESOLUTION_CENSUS: make_resolution_census,
    IN_SEROTYPE_SUMMARY: make_serotype_summary,
    IN_DOMAIN_REPRODUCIBILITY: make_domain_reproducibility,
    IN_SCALE_CURVE: make_scale_curve,
    IN_SIGNIFICANCE_SCREEN: make_significance_screen,
    IN_VARIANCE_BUDGET: make_variance_budget,
    IN_DOMAIN_EFFECT_SUMMARY: make_domain_effect_summary,
    IN_DOMAIN_SEROTYPE_MATRIX: make_domain_serotype_matrix,
    IN_POSITION_CONSERVATION: make_position_conservation,
    IN_REPLICATE_REGIME: make_replicate_regime,
    IN_REPLICATE_BLOCKED_ANALYSES: make_replicate_blocked_analyses,
}


def make_all_inputs() -> dict[str, pd.DataFrame]:
    """Every input table keyed by its filename."""
    return {name: builder() for name, builder in _BUILDERS.items()}


def make_empty_inputs() -> dict[str, pd.DataFrame]:
    """Every input table present but empty, with the columns S7 requires."""
    empty: dict[str, pd.DataFrame] = {}
    for name, full in make_all_inputs().items():
        empty[name] = full.iloc[0:0].copy()
    return empty


def _write(frames: dict[str, pd.DataFrame], root: Path) -> dict[str, Path]:
    stage_dirs: dict[str, Path] = {}
    for stage in set(STAGE_INPUTS.values()):
        d = root / stage
        d.mkdir(parents=True, exist_ok=True)
        stage_dirs[stage] = d
    for name in INPUT_REQUIRED_COLUMNS:
        stage = STAGE_INPUTS[name]
        frames[name].to_parquet(stage_dirs[stage] / name, index=False)
    return stage_dirs


def write_inputs(root: Path) -> dict[str, Path]:
    """Write the full synthetic inputs under ``root`` and return the stage dirs."""
    return _write(make_all_inputs(), root)


def write_empty_inputs(root: Path) -> dict[str, Path]:
    """Write present-but-empty inputs under ``root`` and return the stage dirs."""
    return _write(make_empty_inputs(), root)
