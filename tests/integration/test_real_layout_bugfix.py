"""Regression tests for the real-STRIDE-layout bug fixes (S0 → S1A).

These reproduce, on a synthetic-but-native-shaped dataset, the two failures the
real STRIDE dengue dataset exposed:

* **Bug 1** — discovery must accept the native ``per_run/`` container without the
  user renaming it to ``replicates/``.
* **Bug 2** — the replicate canonical join key must be rebuilt as
  ``"{chain}:{canon_resid}"`` (the real ``label`` is a residue-name like
  ``"HIS-51"``), so that S1A can map every replicate residue to a canonical one.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from stride_analysis import run_s0
from stride_analysis._synthetic import write_real_dataset
from stride_s1a import run_s1a


def test_s0_then_s1a_on_native_per_run_layout(tmp_path: Path) -> None:
    root = write_real_dataset(
        tmp_path / "stride_outs",
        ["DENV1", "DENV2"],
        ["1st_run", "2nd_run", "3rd_run"],
        container="per_run",
    )

    # Bug 1: S0 succeeds on the native per_run/ layout, no renaming.
    s0_out = tmp_path / "s0"
    _stride, replicate, report = run_s0(root, s0_out)
    assert report.all_passed

    # Bug 2: canon_label rebuilt from chain:canon_resid, not the resname label.
    rt = pd.read_parquet(s0_out / "replicate_table.parquet")
    st = pd.read_parquet(s0_out / "stride_table.parquet")
    assert set(rt["canon_label"]) == {"NS3:51", "NS3:200"}
    assert "HIS-51" not in set(rt["canon_label"])
    prof_labels = set(st.loc[st["scale_level"] == "residue", "canon_label"])
    assert set(rt["canon_label"]) <= prof_labels

    # Bug 2 (consumer): S1A now succeeds — no "map to no canonical residue".
    s1a_out = tmp_path / "s1a"
    _tables, s1a_report = run_s1a(
        s0_out / "stride_table.parquet",
        s0_out / "replicate_table.parquet",
        s1a_out,
    )
    assert s1a_report.all_passed
    assert (s1a_out / "canonical_residues.parquet").is_file()


def test_s0_strict_cross_level_passes_on_native_layout(tmp_path: Path) -> None:
    """With the key fixed, even strict cross-level alignment passes on real data."""
    root = write_real_dataset(
        tmp_path / "stride_outs",
        ["DENV1", "DENV2"],
        ["1st_run", "2nd_run", "3rd_run"],
        container="per_run",
    )
    s0_out = tmp_path / "s0"
    _stride, _replicate, report = run_s0(root, s0_out, strict_cross_level=True)
    assert report.all_passed
