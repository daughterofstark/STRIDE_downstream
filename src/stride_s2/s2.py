"""Stage-S2 orchestration.

Thin composition of the reusable subpackages: load the S0 STRIDE table + the two
S1B annotation tables → build the five reduction tables over a ρ* band →
validate the reduction layer → (optionally) write artifacts. No figures, no
cross-serotype tests (that is S5), no calibrated pass/fail claims (the gate is
uncalibrated, §0.1).

Public entry points:

- :func:`build_s2` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s2`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_domain_reproducibility,
    build_residue_landscape,
    build_resolution_census,
    build_serotype_summary,
    build_signed_screen,
)
from .build._screens import normalise_band
from .io import (
    file_digest,
    load_domain_annotation,
    load_residue_annotation,
    load_stride_table,
    write_reduction_summary,
    write_tables,
)
from .models import S2Report
from .models.errors import ConfigError
from .models.schema import (
    DEFAULT_RHO_STAR_BAND,
    PROVISIONAL_RHO_STAR,
    RHO_STAR_DECIMALS,
)
from .validation import (
    validate_census_totals,
    validate_regating_monotonicity,
    validate_serotype_summary_consistency,
    validate_tiers,
    validate_unique_keys,
)


@dataclass
class S2Tables:
    """The five S2 reduction tables, returned together."""

    resolution_census: pd.DataFrame
    residue_landscape: pd.DataFrame
    domain_reproducibility: pd.DataFrame
    signed_screen: pd.DataFrame
    serotype_summary: pd.DataFrame


def _resolve_band(
    rho_star_band: Iterable[float] | None,
) -> tuple[float, ...]:
    """Normalise the requested ρ* band, or fall back to the frozen default."""
    raw = DEFAULT_RHO_STAR_BAND if rho_star_band is None else tuple(rho_star_band)
    band = normalise_band(raw, RHO_STAR_DECIMALS)
    if not band:
        raise ConfigError("rho_star_band is empty; supply at least one ρ* value")
    for value in band:
        if not (0.0 <= value <= 1.0):
            raise ConfigError(
                f"ρ* value {value} is outside [0, 1]; ρ is a reproducibility "
                "coefficient defined on [0, 1]"
            )
    return band


def build_s2(
    stride_table_path: str | Path,
    residue_annotation_path: str | Path,
    domain_annotation_path: str | Path,
    *,
    rho_star_band: Iterable[float] | None = None,
) -> tuple[S2Tables, S2Report]:
    """Load inputs, build the five reduction tables, and validate them.

    Parameters
    ----------
    stride_table_path
        Path to the S0 ``stride_table.parquet``.
    residue_annotation_path, domain_annotation_path
        Paths to the two S1B annotation parquet tables.
    rho_star_band
        The ρ* band to sweep. ``None`` (default) uses the frozen
        :data:`~stride_s2.models.schema.DEFAULT_RHO_STAR_BAND` (0.50–0.90 in
        steps of 0.05). Values are rounded, de-duplicated and sorted; every
        value must lie in ``[0, 1]``.

    Returns
    -------
    (tables, report)
        The five reduction tables and the run report (with provenance header).

    Raises
    ------
    S2Error
        Any subclass on the first input, config, or consistency problem.
    """
    report = S2Report()
    band = _resolve_band(rho_star_band)
    report.rho_star_band = list(band)

    stride_table = load_stride_table(stride_table_path)
    residue_annotation = load_residue_annotation(residue_annotation_path)
    domain_annotation = load_domain_annotation(domain_annotation_path)

    report.serotypes = sorted(
        stride_table["serotype"].astype(str).unique().tolist()
    )
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": PROVISIONAL_RHO_STAR,
        "rho_star_band": list(band),
        "n_replicates_note": (
            "K per serotype is a Level-1 fact (replicate table); at K=3 only "
            "domain-scale claims are licensed, residue-scale is exploratory"
        ),
        "inputs": {
            "stride_table": {
                "path": str(stride_table_path),
                "sha256": file_digest(stride_table_path),
            },
            "residue_annotation": {
                "path": str(residue_annotation_path),
                "sha256": file_digest(residue_annotation_path),
            },
            "domain_annotation": {
                "path": str(domain_annotation_path),
                "sha256": file_digest(domain_annotation_path),
            },
        },
    }
    report.add(
        "inputs loaded",
        "global",
        True,
        f"stride_table={len(stride_table)}, "
        f"residue_annotation={len(residue_annotation)}, "
        f"domain_annotation={len(domain_annotation)}, "
        f"ρ* band={list(band)}",
    )

    # -- build --------------------------------------------------------------
    resolution_census = build_resolution_census(
        stride_table, band, RHO_STAR_DECIMALS
    )
    residue_landscape = build_residue_landscape(stride_table, residue_annotation)
    domain_reproducibility = build_domain_reproducibility(
        stride_table, domain_annotation
    )
    signed_screen = build_signed_screen(stride_table, band, RHO_STAR_DECIMALS)
    serotype_summary = build_serotype_summary(
        stride_table,
        resolution_census,
        residue_landscape,
        signed_screen,
        band,
        RHO_STAR_DECIMALS,
    )

    # -- validate (structural / arithmetic only) ----------------------------
    validate_unique_keys(
        resolution_census,
        residue_landscape,
        domain_reproducibility,
        signed_screen,
        serotype_summary,
        report,
    )
    validate_census_totals(resolution_census, residue_landscape, report)
    validate_regating_monotonicity(resolution_census, report)
    validate_serotype_summary_consistency(
        serotype_summary, resolution_census, report
    )
    validate_tiers(residue_landscape, domain_reproducibility, report)

    # -- report facts -------------------------------------------------------
    report.n_resolution_census = int(len(resolution_census))
    report.n_residue_landscape = int(len(residue_landscape))
    report.n_domain_reproducibility = int(len(domain_reproducibility))
    report.n_signed_screen = int(len(signed_screen))
    report.n_serotype_summary = int(len(serotype_summary))
    report.facts = _facts_at_provisional(serotype_summary)
    report.add(
        "s2 reduction tables built",
        "global",
        True,
        f"{report.n_resolution_census} census, "
        f"{report.n_residue_landscape} residue, "
        f"{report.n_domain_reproducibility} domain, "
        f"{report.n_signed_screen} signed-screen, "
        f"{report.n_serotype_summary} serotype rows",
    )

    tables = S2Tables(
        resolution_census=resolution_census,
        residue_landscape=residue_landscape,
        domain_reproducibility=domain_reproducibility,
        signed_screen=signed_screen,
        serotype_summary=serotype_summary,
    )
    return tables, report


def run_s2(
    stride_table_path: str | Path,
    residue_annotation_path: str | Path,
    domain_annotation_path: str | Path,
    output_dir: str | Path,
    *,
    rho_star_band: Iterable[float] | None = None,
) -> tuple[S2Tables, S2Report]:
    """Full S2: build + validate, then write the six artifacts to ``output_dir``.

    Writes ``resolution_census.parquet``, ``residue_landscape.parquet``,
    ``domain_reproducibility.parquet``, ``signed_screen.parquet``,
    ``serotype_summary.parquet`` and ``reduction_summary.json``. Returns
    ``(tables, report)``.
    """
    tables, report = build_s2(
        stride_table_path,
        residue_annotation_path,
        domain_annotation_path,
        rho_star_band=rho_star_band,
    )
    write_tables(
        tables.resolution_census,
        tables.residue_landscape,
        tables.domain_reproducibility,
        tables.signed_screen,
        tables.serotype_summary,
        output_dir,
    )
    write_reduction_summary(report, output_dir)
    return tables, report


def _facts_at_provisional(serotype_summary: pd.DataFrame) -> dict[str, object]:
    """Compact per-serotype facts at the provisional ρ* (for the summary JSON)."""
    if serotype_summary.empty:
        return {}
    prov = serotype_summary[
        serotype_summary["is_provisional_rho_star"].astype(bool)
    ]
    if prov.empty:
        return {}
    per_serotype = {}
    for row in prov.itertuples(index=False):
        per_serotype[str(row.serotype)] = {
            "n_loci": int(row.n_loci),
            "n_gated_residue": int(row.n_gated_residue),
            "n_mechanisms": int(row.n_mechanisms),
            "n_mixed": int(row.n_mixed),
            "n_signed_significant": int(row.n_signed_significant),
            "rho_residue_median": _round_or_none(row.rho_residue_median),
        }
    return {"at_provisional_rho_star": per_serotype}


def _round_or_none(value: float) -> float | None:
    if value is None or pd.isna(value):
        return None
    return round(float(value), 4)
