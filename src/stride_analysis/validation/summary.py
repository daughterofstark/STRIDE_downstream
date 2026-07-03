"""Level-2 (summary) schema validation and profile↔mechanism consistency.

``validate_profile_schema`` checks column presence, dtypes, value ranges,
scale agreement, hierarchy well-formedness, per-locus structure, and the
one-gated-row-per-locus invariant. ``check_profile_mechanism_consistency``
verifies the cross-file invariants (gated ρ match, exact loci partition, no
orphans on either side).
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pandas as pd

from ..models import MechanismFile
from ..models.errors import ConsistencyError, SchemaError
from ..models.schema import (
    FLOAT_TOL,
    N_SCALES,
    PROFILE_BOOL_COLUMNS,
    PROFILE_COLUMNS,
    PROFILE_FLOAT_COLUMNS,
    PROFILE_INT_COLUMNS,
    PROFILE_STR_COLUMNS,
    SCALE_INDEX_TO_LEVEL,
    SCALE_LEVEL_TO_INDEX,
)
from .hierarchy import validate_path_depth


def validate_profile_schema(df: pd.DataFrame, serotype: str) -> None:
    """Validate a Level-2 profile against the frozen schema."""
    where = f"[{serotype} profile]"

    actual = list(df.columns)
    missing = [c for c in PROFILE_COLUMNS if c not in actual]
    extra = [c for c in actual if c not in PROFILE_COLUMNS]
    if missing:
        raise SchemaError(f"{where} missing required columns: {missing}")
    if extra:
        raise SchemaError(f"{where} unexpected extra columns (drift): {extra}")

    for col in PROFILE_STR_COLUMNS:
        if df[col].isna().any():
            raise SchemaError(f"{where} column {col!r} has missing value(s)")
        if (df[col].astype(str).str.strip() == "").any():
            raise SchemaError(f"{where} column {col!r} has blank value(s)")
    if df["canon_label"].isna().any():
        raise SchemaError(f"{where} canon_label missing on some rows")

    for col in PROFILE_FLOAT_COLUMNS:
        coerced = pd.to_numeric(df[col], errors="coerce")
        bad = coerced.isna() & df[col].notna()
        if bad.any():
            raise SchemaError(
                f"{where} column {col!r} has non-numeric value(s) at rows "
                f"{list(df.index[bad][:5])}"
            )
    for col in PROFILE_INT_COLUMNS:
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.isna().any() or (coerced % 1 != 0).any():
            raise SchemaError(f"{where} column {col!r} must be integer-valued")
    for col in PROFILE_BOOL_COLUMNS:
        if df[col].dtype != bool and not set(df[col].dropna().unique()) <= {
            True, False,
        }:
            raise SchemaError(f"{where} column {col!r} must be boolean")

    for col in ("rho", "coherence"):
        vals = df[col].to_numpy(dtype=float)
        out = (vals < 0.0) | (vals > 1.0)
        if out.any():
            raise SchemaError(
                f"{where} column {col!r} has value(s) outside [0, 1] at rows "
                f"{list(df.index[out][:5])}"
            )
    for col in ("beta", "beta_se", "tau2", "sigma2_bar", "a_signed"):
        vals = df[col].to_numpy(dtype=float)
        if not np.isfinite(vals).all():
            raise SchemaError(
                f"{where} column {col!r} has non-finite value(s) at rows "
                f"{list(df.index[~np.isfinite(vals)][:5])}"
            )
    for col in ("beta_se", "tau2", "sigma2_bar"):
        vals = df[col].to_numpy(dtype=float)
        if (vals < 0).any():
            raise SchemaError(
                f"{where} column {col!r} has negative value(s) at rows "
                f"{list(df.index[vals < 0][:5])}"
            )

    for idx, row in df[["scale_index", "scale_level"]].iterrows():
        si = int(row["scale_index"])
        if si not in SCALE_INDEX_TO_LEVEL:
            raise SchemaError(f"{where} row {idx}: scale_index={si} unknown")
        if SCALE_INDEX_TO_LEVEL[si] != row["scale_level"]:
            raise SchemaError(
                f"{where} row {idx}: scale_index={si} implies "
                f"{SCALE_INDEX_TO_LEVEL[si]!r} but scale_level="
                f"{row['scale_level']!r}"
            )

    for idx, row in df[["region_id", "scale_level"]].iterrows():
        try:
            validate_path_depth(row["region_id"], row["scale_level"])
        except Exception as exc:
            raise SchemaError(f"{where} row {idx}: {exc}") from exc

    dup = df.duplicated(subset=["locus", "scale_level"], keep=False)
    if dup.any():
        dups = (
            df.loc[dup, ["locus", "scale_level"]]
            .drop_duplicates().head(5).to_dict("records")
        )
        raise SchemaError(f"{where} duplicate (locus, scale_level) rows: {dups}")

    per_locus = df.groupby("locus")["scale_level"].nunique()
    bad = per_locus[per_locus != N_SCALES]
    if len(bad):
        examples = {k: int(v) for k, v in bad.head(3).items()}
        raise SchemaError(
            f"{where} {len(bad)} locus/loci not having exactly {N_SCALES} "
            f"scales, e.g. {examples}"
        )

    res = df[df["scale_level"] == "residue"]
    mism = res[res["region_id"] != res["locus"]]
    if len(mism):
        ex = mism.iloc[0]
        raise SchemaError(
            f"{where} {len(mism)} residue row(s) with region_id != locus, e.g. "
            f"region_id={ex['region_id']!r} locus={ex['locus']!r}"
        )

    gated_counts = df[df["gated"]].groupby("locus").size()
    loci = set(df["locus"].unique())
    ungated = loci - set(gated_counts.index)
    if ungated:
        raise SchemaError(
            f"{where} {len(ungated)} locus/loci with no gated row, e.g. "
            f"{sorted(ungated)[:3]}"
        )
    multi = gated_counts[gated_counts > 1]
    if len(multi):
        examples = {k: int(v) for k, v in multi.head(3).items()}
        raise SchemaError(
            f"{where} {len(multi)} locus/loci with >1 gated row, e.g. "
            f"{examples}"
        )

    proteins = set(df["protein"].unique())
    if proteins != {serotype}:
        raise SchemaError(
            f"{where} protein column {sorted(proteins)}; expected all "
            f"== {serotype!r}"
        )


def check_profile_mechanism_consistency(
    df: pd.DataFrame, mech: MechanismFile, serotype: str
) -> None:
    """Verify the profile↔mechanism invariants for one serotype."""
    where = f"[{serotype}]"
    gated = df[df["gated"]].copy()
    gated_by_locus = gated.set_index("locus")
    gated_loci = set(gated_by_locus.index)

    mech_loci: list[str] = []
    for m in mech.mechanisms:
        mech_loci.extend(m.loci)
    mech_loci_set = set(mech_loci)

    if len(mech_loci) != len(mech_loci_set):
        counts = Counter(mech_loci)
        dups = sorted(locus for locus, n in counts.items() if n > 1)
        raise ConsistencyError(
            f"{where} a locus appears in >1 mechanism: {dups[:5]}"
        )

    orphan_mech = mech_loci_set - gated_loci
    orphan_prof = gated_loci - mech_loci_set
    if orphan_mech:
        raise ConsistencyError(
            f"{where} {len(orphan_mech)} mechanism loci have no gated profile "
            f"row (orphan mechanisms), e.g. {sorted(orphan_mech)[:3]}"
        )
    if orphan_prof:
        raise ConsistencyError(
            f"{where} {len(orphan_prof)} gated loci not covered by any "
            f"mechanism (orphan profile rows), e.g. {sorted(orphan_prof)[:3]}"
        )

    for m in mech.mechanisms:
        for loc in m.loci:
            grow = gated_by_locus.loc[loc]
            if grow["region_id"] != m.region_id:
                raise ConsistencyError(
                    f"{where} locus {loc!r} gates to {grow['region_id']!r} but "
                    f"mechanism claims {m.region_id!r}"
                )
            if grow["scale_level"] != m.scale_level:
                raise ConsistencyError(
                    f"{where} locus {loc!r} gated scale {grow['scale_level']!r} "
                    f"!= mechanism scale {m.scale_level!r}"
                )
            if not math.isclose(
                float(grow["rho"]), m.rho, abs_tol=FLOAT_TOL, rel_tol=0.0
            ):
                raise ConsistencyError(
                    f"{where} locus {loc!r}: gated rho={grow['rho']} != "
                    f"mechanism rho={m.rho}"
                )

    n_profile_loci = df["locus"].nunique()
    if mech.summary.n_loci != n_profile_loci:
        raise ConsistencyError(
            f"{where} summary.n_loci={mech.summary.n_loci} != distinct "
            f"profile loci={n_profile_loci}"
        )

    for m in mech.mechanisms:
        expected = SCALE_LEVEL_TO_INDEX.get(m.scale_level)
        if expected is None:
            raise ConsistencyError(
                f"{where} mechanism {m.region_id!r} unknown scale_level "
                f"{m.scale_level!r}"
            )
        if m.scale_index != expected:
            raise ConsistencyError(
                f"{where} mechanism {m.region_id!r}: scale_index={m.scale_index}"
                f" but scale_level {m.scale_level!r} implies {expected}"
            )
