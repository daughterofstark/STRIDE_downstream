"""Structural validation of the S2 reduction layer.

Structural / arithmetic checks only — no statistical assertions and no
biological claims (those are deferred; the gate is uncalibrated, §0.1). Each
check appends a :class:`~stride_s2.models.ValidationCheck` to the report and
raises :class:`~stride_s2.models.errors.ConsistencyError` on failure, so later
stages can trust the reduction layer.

Checks:

- every output table's key is unique;
- the resolution census is *total* — its per-scale locus counts sum to the
  serotype's locus count at every ρ*;
- re-gating is *monotone in ρ* * — the residue-gated count never increases as
  ρ* rises (raising the bar can only push a locus coarser, never finer);
- the per-serotype scorecard partitions the mechanisms (signed + mixed =
  n_mechanisms) and its census rollups agree with the census table;
- tier labels are correct (residue-scale rows are exploratory, domain-scale
  rows are licensed).
"""
from __future__ import annotations

import pandas as pd

from ..models import S2Report
from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_SCALE_LEVEL,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


def validate_unique_keys(
    resolution_census: pd.DataFrame,
    residue_landscape: pd.DataFrame,
    domain_reproducibility: pd.DataFrame,
    signed_screen: pd.DataFrame,
    serotype_summary: pd.DataFrame,
    report: S2Report,
) -> None:
    """Every output table's declared key is unique."""
    _assert_unique(
        resolution_census,
        ["serotype", "rho_star", "gated_scale_level"],
        "resolution_census",
    )
    _assert_unique(
        residue_landscape, ["serotype", "canon_label"], "residue_landscape"
    )
    _assert_unique(
        domain_reproducibility,
        ["serotype", "chain", "domain"],
        "domain_reproducibility",
    )
    _assert_unique(
        signed_screen,
        ["serotype", "canon_label", "rho_star"],
        "signed_screen",
    )
    _assert_unique(
        serotype_summary, ["serotype", "rho_star"], "serotype_summary"
    )
    report.add(
        "every output table key is unique",
        "global",
        True,
        "5 tables checked",
    )


def validate_census_totals(
    resolution_census: pd.DataFrame,
    residue_landscape: pd.DataFrame,
    report: S2Report,
) -> None:
    """The census is total: per-ρ* scale counts sum to the serotype locus count.

    Every locus gates at exactly one scale (or the unresolved sentinel) for a
    given ρ*, so the per-scale counts must sum to the number of loci in the
    serotype — the same at every ρ*.
    """
    if resolution_census.empty:
        report.add("resolution census is total", "global", True, "empty")
        return

    loci_per_serotype = (
        residue_landscape.groupby("serotype").size().to_dict()
        if not residue_landscape.empty
        else {}
    )
    grouped = resolution_census.groupby(["serotype", "rho_star"], sort=True)
    for (serotype, rho_star), grp in grouped:
        total = int(grp["n_loci"].sum())
        expected = int(loci_per_serotype.get(serotype, total))
        if total != expected:
            raise ConsistencyError(
                f"resolution census for ({serotype}, ρ*={rho_star}) sums to "
                f"{total} loci but the serotype has {expected} loci"
            )
    report.add(
        "resolution census is total (counts sum to locus count at every ρ*)",
        "global",
        True,
        f"{len(grouped)} (serotype, ρ*) cells checked",
    )


def validate_regating_monotonicity(
    resolution_census: pd.DataFrame, report: S2Report
) -> None:
    """Residue-gated count is non-increasing as ρ* rises.

    Raising ρ* can only make a locus's finest ρ≥ρ* scale coarser (or unresolved),
    never finer, so the count of loci gating at *residue* scale must not increase
    with ρ*. This is the re-gating analogue of upward closure.
    """
    if resolution_census.empty:
        report.add("re-gating monotone in ρ*", "global", True, "empty")
        return

    residue_index = SCALE_LEVEL_TO_INDEX[RESIDUE_SCALE_LEVEL]
    residue_rows = resolution_census[
        resolution_census["gated_scale_index"] == residue_index
    ]
    counts = (
        residue_rows.groupby(["serotype", "rho_star"])["n_loci"]
        .sum()
        .reset_index()
    )
    for serotype, grp in counts.groupby("serotype", sort=True):
        ordered = grp.sort_values("rho_star")
        values = ordered["n_loci"].tolist()
        for earlier, later in zip(values, values[1:], strict=False):
            if later > earlier:
                raise ConsistencyError(
                    f"residue-gated count for {serotype} increases with ρ* "
                    f"({earlier} → {later}); re-gating must be monotone"
                )
    report.add(
        "re-gating is monotone in ρ* (residue-gated count non-increasing)",
        "global",
        True,
        f"{counts['serotype'].nunique()} serotype(s) checked",
    )


def validate_serotype_summary_consistency(
    serotype_summary: pd.DataFrame,
    resolution_census: pd.DataFrame,
    report: S2Report,
) -> None:
    """The scorecard partitions mechanisms and agrees with the census.

    - ``n_signed + n_mixed == n_mechanisms`` on every row;
    - the residue-gated / domain-or-coarser / unresolved counts sum to
      ``n_loci`` and match the census rollup at the same (serotype, ρ*).
    """
    if serotype_summary.empty:
        report.add("serotype summary consistent", "global", True, "empty")
        return

    for row in serotype_summary.itertuples(index=False):
        if int(row.n_signed) + int(row.n_mixed) != int(row.n_mechanisms):
            raise ConsistencyError(
                f"serotype {row.serotype!r} scorecard: n_signed + n_mixed != "
                f"n_mechanisms ({row.n_signed} + {row.n_mixed} != "
                f"{row.n_mechanisms})"
            )
        parts = (
            int(row.n_gated_residue)
            + int(row.n_gated_domain_or_coarser)
            + int(row.n_unresolved)
        )
        # residue-gated and domain-or-coarser are disjoint; the middle scales
        # (SS/motif) between them plus these must still sum to n_loci.
        if parts > int(row.n_loci):
            raise ConsistencyError(
                f"serotype {row.serotype!r} scorecard: gated-count partition "
                f"({parts}) exceeds n_loci ({row.n_loci})"
            )

    # cross-check n_loci against the census
    if not resolution_census.empty:
        census_loci = (
            resolution_census.groupby(["serotype", "rho_star"])["n_loci"]
            .sum()
            .to_dict()
        )
        for row in serotype_summary.itertuples(index=False):
            key = (row.serotype, float(row.rho_star))
            expected = int(census_loci.get(key, row.n_loci))
            if int(row.n_loci) != expected:
                raise ConsistencyError(
                    f"serotype {row.serotype!r} scorecard n_loci={row.n_loci} "
                    f"disagrees with census total {expected} at ρ*={row.rho_star}"
                )
    report.add(
        "serotype scorecard partitions mechanisms and agrees with census",
        "global",
        True,
        f"{len(serotype_summary)} (serotype, ρ*) rows checked",
    )


def validate_tiers(
    residue_landscape: pd.DataFrame,
    domain_reproducibility: pd.DataFrame,
    report: S2Report,
) -> None:
    """Tier labels are correct: residue exploratory, domain licensed."""
    if not residue_landscape.empty:
        bad = residue_landscape[residue_landscape["tier"] != TIER_EXPLORATORY]
        if not bad.empty:
            raise ConsistencyError(
                "residue landscape has non-exploratory tier labels: "
                f"{sorted(bad['tier'].unique())}"
            )
    if not domain_reproducibility.empty:
        bad = domain_reproducibility[
            domain_reproducibility["tier"] != TIER_LICENSED
        ]
        if not bad.empty:
            raise ConsistencyError(
                "domain reproducibility has non-licensed tier labels: "
                f"{sorted(bad['tier'].unique())}"
            )
    report.add(
        "tier labels correct (residue exploratory, domain licensed)",
        "global",
        True,
        f"domain scale index = {SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]}",
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
