"""Report carriers for S7 (plain dataclasses serialised to report_summary.json)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ValidationCheck:
    name: str
    scope: str
    passed: bool
    detail: str = ""


@dataclass
class ArtifactRecord:
    """One generated artifact (figure or table) and its on-disk files."""

    artifact_id: str  # e.g. "F1" or "T1"
    kind: str  # "figure" | "table"
    title: str
    slug: str
    sources: list[str]  # source table filenames it was assembled from
    files: list[str]  # output filenames written (relative to output dir)
    n_rows: int  # rows in the prepared / assembled table
    tier: str = ""  # dominant tier label, when the source rows carry one


@dataclass
class S7Report:
    """Aggregate S7 report: generated artifacts, provenance, and checks.

    ``provenance`` carries the design's required per-run header (§5.4): the input
    file digests, the provisional ρ\\*, ``calibrated = false``, and the K = 3 note.
    ``limitations`` records the replicate-layer status assembled from S6 (which no
    design figure/table draws on) so the blocked per-run analyses are documented
    rather than silently dropped.
    """

    checks: list[ValidationCheck] = field(default_factory=list)
    figures: list[ArtifactRecord] = field(default_factory=list)
    tables: list[ArtifactRecord] = field(default_factory=list)
    serotypes: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    facts: dict[str, Any] = field(default_factory=dict)
    limitations: list[dict[str, Any]] = field(default_factory=list)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
