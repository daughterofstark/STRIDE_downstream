"""STRIDE summary canonical table construction.

Builds the Level-2 canonical table: one row per (serotype, canon_label,
scale_level), merging profile + explicit hierarchy + mechanism payload (on gated
rows only) + provenance. Kept **separate** from the replicate table.

Unique key: ``(serotype, canon_label, scale_level)``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..models import MechanismFile
from ..models.errors import ConsistencyError
from ..models.schema import HIERARCHY_COLUMNS, STRIDE_TABLE_COLUMNS
from ..validation.hierarchy import parse_residue_path


def _hierarchy_frame(df: pd.DataFrame) -> pd.DataFrame:
    residue = df[df["scale_level"] == "residue"][["locus", "region_id"]]
    records = []
    for locus, region_id in zip(
        residue["locus"], residue["region_id"], strict=True
    ):
        parsed = parse_residue_path(region_id)
        parsed["locus"] = locus
        records.append(parsed)
    return pd.DataFrame.from_records(records)[["locus", *HIERARCHY_COLUMNS]]


def _mechanism_frame(mech: MechanismFile) -> pd.DataFrame:
    records = []
    for m in mech.mechanisms:
        payload = {
            "mech_label": m.label,
            "mech_direction": m.direction,
            "mech_beta_signed": m.beta_signed,
            "mech_beta_ci_lower": m.beta_ci_lower,
            "mech_beta_ci_upper": m.beta_ci_upper,
            "mech_beta_se": m.beta_se,
            "mech_coherence": m.coherence,
            "mech_reproducible_magnitude_energy": m.reproducible_magnitude_energy,
            "mech_rho_star": m.rho_star,
            "mech_calibrated": m.calibrated,
            "mech_gate_uncertain": m.gate_uncertain,
            "mech_status": m.status,
            "mech_region_id": m.region_id,
            "mech_n_loci": m.n_loci,
        }
        for locus in m.loci:
            records.append({"locus": locus, **payload})
    return pd.DataFrame.from_records(records)


def build_stride_rows(
    df: pd.DataFrame,
    mech: MechanismFile,
    serotype: str,
    profile_source: str | Path,
    mechanism_source: str | Path,
) -> pd.DataFrame:
    """Build the STRIDE canonical rows for one serotype.

    Assumes ``df``/``mech`` already passed validation + consistency.
    """
    out = df.copy().rename(columns={"protein": "serotype"})

    hier = _hierarchy_frame(df)
    out = out.merge(hier, on="locus", how="left", validate="many_to_one")

    out["is_gated_scale"] = out["gated"]

    mech_df = _mechanism_frame(mech)
    out = out.merge(mech_df, on="locus", how="left", validate="many_to_one")
    # mechanism payload belongs to the gated scale only; null it elsewhere in a
    # dtype-safe way (object dtype honestly represents "present only on ℓ̂*").
    mech_cols = [c for c in mech_df.columns if c != "locus"]
    non_gated = ~out["is_gated_scale"]
    for col in mech_cols:
        out[col] = out[col].astype(object)
        out.loc[non_gated, col] = None

    out["profile_source"] = str(profile_source)
    out["mechanism_source"] = str(mechanism_source)
    out["gate_rho_star"] = mech.gate.rho_star
    out["gate_alpha"] = mech.gate.alpha
    out["gate_coherence_threshold"] = mech.gate.coherence_threshold
    out["mechanism_calibrated"] = mech.calibrated
    out["mechanism_schema_version"] = mech.schema_version

    missing = [c for c in STRIDE_TABLE_COLUMNS if c not in out.columns]
    if missing:  # pragma: no cover - internal guard
        raise ConsistencyError(
            f"[{serotype}] internal: stride table missing columns {missing}"
        )
    return out[list(STRIDE_TABLE_COLUMNS)]


def assemble_stride_table(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate per-serotype STRIDE frames and assert key uniqueness."""
    if not frames:
        return pd.DataFrame(columns=list(STRIDE_TABLE_COLUMNS))
    master = pd.concat(frames, ignore_index=True)
    key = ["serotype", "canon_label", "scale_level"]
    dup = master.duplicated(subset=key, keep=False)
    if dup.any():
        examples = master.loc[dup, key].drop_duplicates().head(5).to_dict("records")
        raise ConsistencyError(
            f"stride key {key} not unique; duplicate examples: {examples}"
        )
    return master.sort_values(
        ["serotype", "scale_index", "canon_label"]
    ).reset_index(drop=True)
