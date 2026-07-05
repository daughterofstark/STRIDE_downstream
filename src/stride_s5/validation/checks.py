"""Structural validation of the S5 cross-serotype layer.

Structural / arithmetic checks only — no statistical assertions and no biological
claims (the gate is uncalibrated, and serotype is the unit of replication at
n = 4, §0.1, §5.2). Each check appends a
:class:`~stride_s5.models.ValidationCheck` to the report and raises
:class:`~stride_s5.models.errors.ConsistencyError` on failure, so later stages can
trust the layer.

Checks:

- every output table's key is unique;
- the position-conservation table is *arithmetically consistent* — the
  reproducible / signed counts are bounded (``0 ≤ signed ≤ reproducible ≤
  present``), the reproducible fraction matches the counts, the
  ``conservation_class`` label matches the recomputed counts, the divergent flag
  matches its definition, and the Catalytic-Triad flag matches the label set;
- the direction-concordance table *partitions* each position's signed serotypes
  (``n_increase + n_decrease == n_serotypes_signed`` with at least the minimum
  required), and the ``concordance_class`` / ``majority_direction`` /
  ``frac_majority`` match the recomputed counts;
- the domain × serotype matrix carries ρ in ``[0, 1]``, is Tier A (licensed), and
  its catalytic-domain flag matches the label set;
- the cross-serotype scorecard *partitions* each serotype's mechanisms
  (``n_signed + n_mixed == n_mechanisms``) and its residue / shared counts are
  bounded, with every fraction in ``[0, 1]``.
"""
from __future__ import annotations

import math

import pandas as pd

from ..build._classify import (
    concordance_class,
    conservation_class,
    is_catalytic_domain,
    is_catalytic_triad,
    majority_direction,
)
from ..models import S5Report
from ..models.errors import ConsistencyError
from ..models.schema import (
    MIN_SEROTYPES_FOR_CONCORDANCE,
    RHO_DECIMALS,
    TIER_LICENSED,
)

_TOL = 10.0 ** (-(RHO_DECIMALS - 1))


def validate_unique_keys(
    position_conservation: pd.DataFrame,
    direction_concordance: pd.DataFrame,
    domain_serotype_matrix: pd.DataFrame,
    cross_serotype_scorecard: pd.DataFrame,
    report: S5Report,
) -> None:
    """Every output table's declared key is unique."""
    _assert_unique(
        position_conservation, ["canon_label"], "position_conservation"
    )
    _assert_unique(
        direction_concordance, ["canon_label"], "direction_concordance"
    )
    _assert_unique(
        domain_serotype_matrix,
        ["serotype", "chain", "domain"],
        "domain_serotype_matrix",
    )
    _assert_unique(
        cross_serotype_scorecard, ["serotype"], "cross_serotype_scorecard"
    )
    report.add(
        "every output table key is unique", "global", True, "4 tables checked"
    )


def validate_position_conservation(
    position_conservation: pd.DataFrame, report: S5Report
) -> None:
    """Counts are bounded and the conservation / divergent / catalytic labels match."""
    if position_conservation.empty:
        report.add(
            "position conservation self-consistent", "global", True, "empty"
        )
        return
    for row in position_conservation.itertuples(index=False):
        key = f"({row.canon_label})"
        n_present = int(row.n_serotypes_present)
        n_repro = int(row.n_serotypes_reproducible)
        n_signed = int(row.n_serotypes_signed_reproducible)
        n_total = int(row.n_serotypes_total)
        if not (0 <= n_signed <= n_repro <= n_present):
            raise ConsistencyError(
                f"position_conservation {key}: counts not ordered "
                f"0 <= signed({n_signed}) <= reproducible({n_repro}) <= "
                f"present({n_present})"
            )
        if n_present > n_total:
            raise ConsistencyError(
                f"position_conservation {key}: n_serotypes_present {n_present} "
                f"exceeds n_serotypes_total {n_total}"
            )
        frac = float(row.frac_reproducible)
        if n_present > 0:
            expected_frac = n_repro / n_present
            if math.isnan(frac) or abs(frac - expected_frac) > _TOL:
                raise ConsistencyError(
                    f"position_conservation {key}: frac_reproducible {frac} "
                    f"!= n_reproducible/n_present ({expected_frac})"
                )
        expected_class = conservation_class(n_repro, n_present)
        if row.conservation_class != expected_class:
            raise ConsistencyError(
                f"position_conservation {key}: conservation_class "
                f"{row.conservation_class!r} disagrees with counts "
                f"(expected {expected_class!r})"
            )
        expected_divergent = 0 < n_signed < n_total
        if bool(row.is_serotype_divergent) != expected_divergent:
            raise ConsistencyError(
                f"position_conservation {key}: is_serotype_divergent "
                f"{bool(row.is_serotype_divergent)} disagrees with "
                f"0 < signed({n_signed}) < total({n_total})"
            )
        if bool(row.is_catalytic_triad) != is_catalytic_triad(
            str(row.canon_label)
        ):
            raise ConsistencyError(
                f"position_conservation {key}: is_catalytic_triad flag "
                f"disagrees with the Catalytic-Triad label set"
            )
    report.add(
        "position conservation self-consistent (counts/class/divergent/catalytic)",
        "global",
        True,
        f"{len(position_conservation)} positions checked",
    )


def validate_direction_concordance(
    direction_concordance: pd.DataFrame, report: S5Report
) -> None:
    """Direction counts partition the signed serotypes; the labels match."""
    if direction_concordance.empty:
        report.add(
            "direction concordance self-consistent", "global", True, "empty"
        )
        return
    for row in direction_concordance.itertuples(index=False):
        key = f"({row.canon_label})"
        n_signed = int(row.n_serotypes_signed)
        n_inc = int(row.n_increase)
        n_dec = int(row.n_decrease)
        if n_inc + n_dec != n_signed:
            raise ConsistencyError(
                f"direction_concordance {key}: n_increase + n_decrease "
                f"({n_inc} + {n_dec}) != n_serotypes_signed ({n_signed})"
            )
        if n_signed < MIN_SEROTYPES_FOR_CONCORDANCE:
            raise ConsistencyError(
                f"direction_concordance {key}: n_serotypes_signed {n_signed} "
                f"below the minimum {MIN_SEROTYPES_FOR_CONCORDANCE}"
            )
        expected_class = concordance_class(n_inc, n_dec)
        if row.concordance_class != expected_class:
            raise ConsistencyError(
                f"direction_concordance {key}: concordance_class "
                f"{row.concordance_class!r} disagrees with counts "
                f"(expected {expected_class!r})"
            )
        expected_dir = majority_direction(n_inc, n_dec)
        if row.majority_direction != expected_dir:
            raise ConsistencyError(
                f"direction_concordance {key}: majority_direction "
                f"{row.majority_direction!r} disagrees with counts "
                f"(expected {expected_dir!r})"
            )
        frac = float(row.frac_majority)
        if not (-_TOL <= frac <= 1.0 + _TOL):
            raise ConsistencyError(
                f"direction_concordance {key}: frac_majority {frac} "
                f"outside [0, 1]"
            )
    report.add(
        "direction concordance self-consistent (counts/class/majority)",
        "global",
        True,
        f"{len(direction_concordance)} positions checked",
    )


def validate_domain_serotype_matrix(
    domain_serotype_matrix: pd.DataFrame, report: S5Report
) -> None:
    """ρ is in [0, 1]; rows are licensed; the catalytic-domain flag matches."""
    if domain_serotype_matrix.empty:
        report.add(
            "domain serotype matrix self-consistent", "global", True, "empty"
        )
        return
    for row in domain_serotype_matrix.itertuples(index=False):
        key = f"({row.serotype}, {row.chain}, {row.domain})"
        rho = float(row.rho_domain)
        if math.isnan(rho) or not (-_TOL <= rho <= 1.0 + _TOL):
            raise ConsistencyError(
                f"domain_serotype_matrix {key}: rho_domain {rho} outside [0, 1]"
            )
        if row.tier != TIER_LICENSED:
            raise ConsistencyError(
                f"domain_serotype_matrix {key}: tier {row.tier!r} is not "
                f"{TIER_LICENSED!r} (domain scale is licensed)"
            )
        if bool(row.is_catalytic_domain) != is_catalytic_domain(
            str(row.domain)
        ):
            raise ConsistencyError(
                f"domain_serotype_matrix {key}: is_catalytic_domain flag "
                f"disagrees with the catalytic-domain label set"
            )
    report.add(
        "domain serotype matrix self-consistent (rho range/tier/catalytic)",
        "global",
        True,
        f"{len(domain_serotype_matrix)} (serotype, chain, domain) rows checked",
    )


def validate_cross_serotype_scorecard(
    cross_serotype_scorecard: pd.DataFrame, report: S5Report
) -> None:
    """Mechanisms partition into signed/mixed; residue/shared counts are bounded."""
    if cross_serotype_scorecard.empty:
        report.add(
            "cross serotype scorecard self-consistent", "global", True, "empty"
        )
        return
    for row in cross_serotype_scorecard.itertuples(index=False):
        key = f"({row.serotype})"
        if int(row.n_signed) + int(row.n_mixed) != int(row.n_mechanisms):
            raise ConsistencyError(
                f"cross_serotype_scorecard {key}: n_signed + n_mixed "
                f"({row.n_signed} + {row.n_mixed}) != n_mechanisms "
                f"({row.n_mechanisms})"
            )
        if int(row.n_reproducible_residue) > int(row.n_loci):
            raise ConsistencyError(
                f"cross_serotype_scorecard {key}: n_reproducible_residue "
                f"{row.n_reproducible_residue} > n_loci {row.n_loci}"
            )
        if int(row.n_shared_reproducible) > int(row.n_shared_positions):
            raise ConsistencyError(
                f"cross_serotype_scorecard {key}: n_shared_reproducible "
                f"{row.n_shared_reproducible} > n_shared_positions "
                f"{row.n_shared_positions}"
            )
        for label, val in (
            ("frac_reproducible_residue", float(row.frac_reproducible_residue)),
            ("frac_signed", float(row.frac_signed)),
            ("frac_mixed", float(row.frac_mixed)),
        ):
            if not math.isnan(val) and not (-_TOL <= val <= 1.0 + _TOL):
                raise ConsistencyError(
                    f"cross_serotype_scorecard {key}: {label}={val} "
                    f"outside [0, 1]"
                )
    report.add(
        "cross serotype scorecard self-consistent (partition/bounds/fractions)",
        "global",
        True,
        f"{len(cross_serotype_scorecard)} serotype rows checked",
    )


def _assert_unique(df: pd.DataFrame, key: list[str], what: str) -> None:
    if df.empty:
        return
    dup = df.duplicated(key, keep=False)
    if dup.any():
        examples = df.loc[dup, key].drop_duplicates().head(3).to_dict("records")
        raise ConsistencyError(
            f"{what} key {tuple(key)} not unique; examples: {examples}"
        )
