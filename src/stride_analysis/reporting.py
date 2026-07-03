"""Report serialisation for S0.

Writes ``schema_report.json`` (machine-readable) and ``validation_report.md``
(human-readable). Structural facts only — no biology, no statistics.
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from .models import Report
from .models.schema import REPLICATE_TABLE_IDENTITY, STRIDE_TABLE_COLUMNS


def write_schema_report(report: Report, path: str | Path) -> None:
    payload = {
        "stage": "S0",
        "all_checks_passed": report.all_passed,
        "replicate_table_rows": report.replicate_rows,
        "stride_table_rows": report.stride_rows,
        "replicate_table_identity_columns": list(REPLICATE_TABLE_IDENTITY),
        "stride_table_columns": list(STRIDE_TABLE_COLUMNS),
        "serotypes": report.serotype_facts,
        "checks": [asdict(c) for c in report.checks],
    }
    Path(path).write_text(json.dumps(payload, indent=2))


def write_validation_report(report: Report, path: str | Path) -> None:
    lines: list[str] = ["# S0 Validation Report", ""]
    lines.append(f"**Overall status: {'PASSED' if report.all_passed else 'FAILED'}**")
    lines.append("")
    lines.append(
        "Stage S0 ingests STRIDE Level-1 replicate observations and Level-2 "
        "summaries and builds two separate canonical tables. This report is "
        "structural only (no biology, no statistics, no figures)."
    )
    lines.append("")

    lines.append("## Datasets ingested")
    lines.append("")
    lines.append(
        "| Serotype | Replicates | Replicate rows | Profile loci | "
        "Mechanisms | Calibrated |"
    )
    lines.append("|---|---|---|---|---|---|")
    for f in report.serotype_facts:
        lines.append(
            f"| {f['serotype']} | {f['n_replicates']} | "
            f"{f['replicate_rows']} | {f.get('profile_loci', '—')} | "
            f"{f.get('n_mechanisms', '—')} | {f.get('calibrated', '—')} |"
        )
    lines.append("")

    lines.append("## Validation & consistency checks")
    lines.append("")
    lines.append("| Check | Scope | Result | Detail |")
    lines.append("|---|---|---|---|")
    for c in report.checks:
        mark = "PASS" if c.passed else "FAIL"
        detail = c.detail.replace("|", "\\|") if c.detail else ""
        lines.append(f"| {c.name} | {c.scope} | {mark} | {detail} |")
    lines.append("")

    lines.append("## Canonical tables")
    lines.append("")
    lines.append(
        f"- **Replicate table**: {report.replicate_rows} rows — key "
        "`(serotype, replicate, canon_label)`."
    )
    lines.append(
        f"- **STRIDE table**: {report.stride_rows} rows — key "
        "`(serotype, canon_label, scale_level)`."
    )
    lines.append("")
    lines.append(
        "The two tables are kept separate; later stages consume whichever "
        "they need. Replicate observations are never collapsed into the "
        "summaries."
    )
    lines.append("")
    Path(path).write_text("\n".join(lines))
