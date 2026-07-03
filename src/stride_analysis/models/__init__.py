"""Structured data models.

- Level 2 mechanism JSON is validated with pydantic (nested, nullable, enum).
- Dataset-layout descriptors and reports are plain dataclasses (pure carriers).

Level 1 replicate tables and Level 2 profile tables are validated as pandas
DataFrames in :mod:`stride_analysis.validation` — columnar data where pandas is
the right tool and per-row models would add cost without safety.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from .schema import VALID_DIRECTIONS


# ---------------------------------------------------------------------------
# Level 2 — mechanism JSON models
# ---------------------------------------------------------------------------
class GateParams(BaseModel):
    model_config = ConfigDict(extra="forbid")
    rho_star: float
    alpha: float
    coherence_threshold: float


class MechanismSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")
    n_loci: int
    n_mechanisms: int
    n_unresolved: int
    n_gate_uncertain: int


class Mechanism(BaseModel):
    """One element of the mechanisms array (18 fields).

    Signed beta fields are Optional because STRIDE emits ``null`` for mixed
    regions; the "null iff mixed" rule is enforced in a model validator.
    """

    model_config = ConfigDict(extra="forbid")

    region_id: str
    label: str
    scale_level: str
    scale_index: int
    n_loci: int
    loci: list[str]
    rho: float
    rho_star: float
    calibrated: bool
    direction: str
    beta_signed: float | None
    beta_ci_lower: float | None
    beta_ci_upper: float | None
    beta_se: float | None
    coherence: float
    reproducible_magnitude_energy: float
    method: str
    gate_uncertain: bool
    status: str

    @field_validator("direction")
    @classmethod
    def _valid_direction(cls, v: str) -> str:
        if v not in VALID_DIRECTIONS:
            raise ValueError(
                f"direction={v!r} not in {sorted(VALID_DIRECTIONS)}"
            )
        return v

    @field_validator("rho", "coherence")
    @classmethod
    def _unit_interval(cls, v: float, info: Any) -> float:
        if not (0.0 <= v <= 1.0):
            raise ValueError(f"{info.field_name}={v} outside [0, 1]")
        return v

    @field_validator("n_loci")
    @classmethod
    def _positive(cls, v: int) -> int:
        if v < 1:
            raise ValueError(f"n_loci={v} must be >= 1")
        return v

    @model_validator(mode="after")
    def _mixed_nullability(self) -> Mechanism:
        beta = {
            "beta_signed": self.beta_signed,
            "beta_ci_lower": self.beta_ci_lower,
            "beta_ci_upper": self.beta_ci_upper,
            "beta_se": self.beta_se,
        }
        all_null = all(v is None for v in beta.values())
        any_null = any(v is None for v in beta.values())
        if self.direction == "mixed" and not all_null:
            present = [k for k, v in beta.items() if v is not None]
            raise ValueError(
                f"region {self.region_id!r}: direction 'mixed' but signed "
                f"beta fields present: {present}"
            )
        if self.direction != "mixed" and any_null:
            missing = [k for k, v in beta.items() if v is None]
            raise ValueError(
                f"region {self.region_id!r}: direction {self.direction!r} but "
                f"signed beta fields null: {missing}"
            )
        if self.n_loci != len(self.loci):
            raise ValueError(
                f"region {self.region_id!r}: n_loci={self.n_loci} != "
                f"len(loci)={len(self.loci)}"
            )
        return self


class MechanismFile(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: str
    calibrated: bool
    uncalibrated_note: str
    gate: GateParams
    summary: MechanismSummary
    mechanisms: list[Mechanism]
    unresolved_loci: list[Any]

    @model_validator(mode="after")
    def _summary_matches(self) -> MechanismFile:
        s = self.summary
        if s.n_mechanisms != len(self.mechanisms):
            raise ValueError(
                f"summary.n_mechanisms={s.n_mechanisms} != "
                f"len(mechanisms)={len(self.mechanisms)}"
            )
        if s.n_unresolved != len(self.unresolved_loci):
            raise ValueError(
                f"summary.n_unresolved={s.n_unresolved} != "
                f"len(unresolved_loci)={len(self.unresolved_loci)}"
            )
        gu = sum(1 for m in self.mechanisms if m.gate_uncertain)
        if s.n_gate_uncertain != gu:
            raise ValueError(
                f"summary.n_gate_uncertain={s.n_gate_uncertain} != count={gu}"
            )
        return self


# ---------------------------------------------------------------------------
# dataset-layout descriptors (plain dataclasses)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ReplicateInput:
    """One discovered Level-1 replicate table for a serotype."""

    serotype: str
    run_dir: str
    replicate_index: int  # 1-based, defined by sorted run-dir order
    correlations_path: Path


@dataclass(frozen=True)
class SummaryInput:
    """One discovered Level-2 profile/mechanism pair for a serotype."""

    serotype: str
    profile_path: Path
    mechanism_path: Path


@dataclass(frozen=True)
class SerotypeDataset:
    """All discovered inputs for a single serotype (both levels)."""

    serotype: str
    replicates: tuple[ReplicateInput, ...]
    summary: SummaryInput | None

    @property
    def n_replicates(self) -> int:
        return len(self.replicates)


@dataclass(frozen=True)
class Dataset:
    """A fully discovered dataset spanning all serotypes."""

    serotypes: tuple[SerotypeDataset, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(s.serotype for s in self.serotypes)


# ---------------------------------------------------------------------------
# report carriers
# ---------------------------------------------------------------------------
@dataclass
class ValidationCheck:
    name: str
    scope: str
    passed: bool
    detail: str = ""


@dataclass
class Report:
    """Aggregate ingest report: dataset facts + every check outcome."""

    checks: list[ValidationCheck] = field(default_factory=list)
    replicate_rows: int = 0
    stride_rows: int = 0
    serotype_facts: list[dict[str, Any]] = field(default_factory=list)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
