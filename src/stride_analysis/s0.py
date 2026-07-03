"""Stage-S0 orchestration.

Thin composition of the reusable subpackages (``io`` → ``validation`` →
``canonical``). There is no monolithic engine: this module wires the pieces and
collects a report. Later stages import the same subpackages directly.

Public entry points:

- :func:`build_tables`  — discover → validate → build; returns the two tables +
  report, no file writes.
- :func:`run_s0`        — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .canonical import (
    assemble_replicate_table,
    assemble_stride_table,
    build_replicate_rows,
    build_stride_rows,
)
from .io import (
    discover_dataset,
    load_correlations,
    load_mechanism,
    load_profile,
)
from .models import Dataset, Report
from .models.schema import (
    OUT_REPLICATE_CSV,
    OUT_REPLICATE_PARQUET,
    OUT_SCHEMA_REPORT,
    OUT_STRIDE_CSV,
    OUT_STRIDE_PARQUET,
    OUT_VALIDATION_REPORT,
)
from .reporting import write_schema_report, write_validation_report
from .validation import (
    check_profile_mechanism_consistency,
    check_replicate_summary_alignment,
    validate_correlations_schema,
    validate_profile_schema,
)


def build_tables(
    data_root: str | Path,
    *,
    require_replicates: bool = True,
    require_summaries: bool = True,
    enforce_equal_replicate_counts: bool = True,
    strict_cross_level: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, Report]:
    """Discover, validate, and build both canonical tables in memory.

    Parameters
    ----------
    data_root
        Dataset root to discover. Any layout the discovery logic recognises
        works; no path is hardcoded. See ``data/README.md``.
    require_replicates, require_summaries
        If ``True`` (default), every serotype must provide that data level;
        set ``False`` to allow partial datasets (e.g. summaries only).
    enforce_equal_replicate_counts
        If ``True`` (default), all serotypes must have the same replicate
        count; a mismatch raises :class:`~stride_analysis.DiscoveryError`.
    strict_cross_level
        If ``True``, promote the replicate↔summary residue-alignment check
        from advisory to hard (raises on any replicate residue label absent
        from the profile). Default ``False``.

    Returns
    -------
    (replicate_table, stride_table, report)
        The Level-1 replicate table, the Level-2 STRIDE table (either may be
        empty for a partial dataset), and the run :class:`~stride_analysis.Report`.

    Raises
    ------
    StrideAnalysisError
        Any subclass on the first discovery, schema, hierarchy, or consistency
        problem encountered.
    """
    report = Report()
    dataset: Dataset = discover_dataset(
        data_root,
        require_replicates=require_replicates,
        require_summaries=require_summaries,
        enforce_equal_replicate_counts=enforce_equal_replicate_counts,
    )
    report.add(
        "discovery",
        "global",
        True,
        f"{len(dataset.names)} serotype(s): {list(dataset.names)}",
    )

    replicate_frames: list[pd.DataFrame] = []
    stride_frames: list[pd.DataFrame] = []

    for sero in dataset.serotypes:
        facts: dict[str, Any] = {
            "serotype": sero.serotype,
            "n_replicates": sero.n_replicates,
            "replicate_rows": 0,
        }

        # -- Level 1: replicate observations ---------------------------------
        rep_labels_union: set[str] = set()
        sero_rep_rows = 0
        for rep in sero.replicates:
            df = load_correlations(rep.correlations_path)
            validate_correlations_schema(df, rep.serotype, rep.run_dir)
            rep_labels_union |= set(df["label"].astype(str))
            rows = build_replicate_rows(df, rep)
            sero_rep_rows += len(rows)
            replicate_frames.append(rows)
        facts["replicate_rows"] = sero_rep_rows
        if sero.replicates:
            report.add(
                "replicate schema",
                sero.serotype,
                True,
                f"{sero.n_replicates} replicate(s), {sero_rep_rows} rows",
            )

        # -- Level 2: STRIDE summaries ---------------------------------------
        if sero.summary is not None:
            prof = load_profile(sero.summary.profile_path)
            mech = load_mechanism(sero.summary.mechanism_path)
            validate_profile_schema(prof, sero.serotype)
            check_profile_mechanism_consistency(prof, mech, sero.serotype)
            report.add(
                "profile↔mechanism consistency",
                sero.serotype,
                True,
                "gated rho match, exact loci partition, no orphans",
            )
            stride_frames.append(
                build_stride_rows(
                    prof,
                    mech,
                    sero.serotype,
                    sero.summary.profile_path,
                    sero.summary.mechanism_path,
                )
            )
            facts["profile_loci"] = int(prof["locus"].nunique())
            facts["n_mechanisms"] = len(mech.mechanisms)
            facts["calibrated"] = mech.calibrated

            # -- cross-level alignment (advisory unless strict) --------------
            if rep_labels_union:
                detail = check_replicate_summary_alignment(
                    rep_labels_union, prof, sero.serotype, strict=strict_cross_level
                )
                report.add(
                    "replicate↔summary residue alignment",
                    sero.serotype,
                    True,
                    detail,
                )

        report.serotype_facts.append(facts)

    replicate_table = assemble_replicate_table(replicate_frames)
    stride_table = assemble_stride_table(stride_frames)
    report.add(
        "replicate key uniqueness (serotype, replicate, canon_label)",
        "global",
        True,
        f"{len(replicate_table)} rows",
    )
    report.add(
        "stride key uniqueness (serotype, canon_label, scale_level)",
        "global",
        True,
        f"{len(stride_table)} rows",
    )
    report.replicate_rows = int(len(replicate_table))
    report.stride_rows = int(len(stride_table))
    return replicate_table, stride_table, report


def run_s0(
    data_root: str | Path,
    output_dir: str | Path,
    *,
    require_replicates: bool = True,
    require_summaries: bool = True,
    enforce_equal_replicate_counts: bool = True,
    strict_cross_level: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, Report]:
    """Full S0: build both tables and write all artifacts to ``output_dir``.

    Parameters
    ----------
    data_root
        Dataset root to discover (see :func:`build_tables`).
    output_dir
        Directory to write the artifacts into; created if absent.
    require_replicates, require_summaries, enforce_equal_replicate_counts, strict_cross_level
        Forwarded to :func:`build_tables`; see there for semantics.

    Returns
    -------
    (replicate_table, stride_table, report)
        The two canonical tables and the run report. Empty tables are not
        written to disk.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    replicate_table, stride_table, report = build_tables(
        data_root,
        require_replicates=require_replicates,
        require_summaries=require_summaries,
        enforce_equal_replicate_counts=enforce_equal_replicate_counts,
        strict_cross_level=strict_cross_level,
    )

    if not replicate_table.empty:
        replicate_table.to_parquet(out / OUT_REPLICATE_PARQUET, index=False)
        replicate_table.to_csv(out / OUT_REPLICATE_CSV, index=False)
    if not stride_table.empty:
        stride_table.to_parquet(out / OUT_STRIDE_PARQUET, index=False)
        stride_table.to_csv(out / OUT_STRIDE_CSV, index=False)
    write_schema_report(report, out / OUT_SCHEMA_REPORT)
    write_validation_report(report, out / OUT_VALIDATION_REPORT)
    return replicate_table, stride_table, report
