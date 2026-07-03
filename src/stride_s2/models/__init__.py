"""Report carriers for S2 (plain dataclasses serialised to reduction_summary.json)."""
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
class S2Report:
    """Aggregate S2 report: reduction facts, provenance, and every check outcome.

    ``provenance`` carries the design's required per-run header (§5.4): input
    file digests, the ρ* band used, ``calibrated=false``, and K where known.
    """

    checks: list[ValidationCheck] = field(default_factory=list)
    n_resolution_census: int = 0
    n_residue_landscape: int = 0
    n_domain_reproducibility: int = 0
    n_signed_screen: int = 0
    n_serotype_summary: int = 0
    serotypes: list[str] = field(default_factory=list)
    rho_star_band: list[float] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    facts: dict[str, Any] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
