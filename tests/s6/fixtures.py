"""Synthetic fixtures for S6 tests.

Builds tiny, S1A-shaped ``replicate_inventory`` and S0-shaped ``replicate_table``
frames directly (no dependence on the real DENV data and no need to run S0/S1A).
The default scenario has the four dengue serotypes over **six** shared positions
(so Kendall's *W* is non-degenerate, unlike the two-position example dataset),
with per-run θ crafted to exercise every concordance branch:

- ``DENV1`` — three runs whose per-position θ preserve an identical ordering
  (each run adds a tiny monotone offset) → Kendall's *W* = 1.0 → ``strong``.
- ``DENV2`` — three runs with a different but again order-preserving θ → ``strong``.
- ``DENV3`` — one run's ordering is reversed relative to the others → a lower *W*
  (``weak`` / ``moderate``), used to check the class tracks the coefficient.
- ``DENV4`` — the same as DENV1 but with position ``NS3:250`` missing from the
  third run (and flagged ``in_all_replicates == False`` in the inventory), so the
  complete-case intersection drops to five positions and completeness is < 1.

A separate :func:`make_empty_inventory` builds the design's *blocked* state — the
inventory S1A emits when no per-run correlation CSVs were supplied (every
``n_replicates == 0``, ``in_all_replicates == False``) — for which S6 must produce
empty per-run tables and mark the analyses blocked.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]

#: canonical positions (canon_label → chain / domain_label)
_POSITIONS: list[tuple[str, str, str]] = [
    ("NS2B:-1", "NS2B", "Cofactor Interface"),
    ("NS3:51", "NS3", "Catalytic Triad"),
    ("NS3:75", "NS3", "Catalytic Triad"),
    ("NS3:135", "NS3", "Catalytic Triad"),
    ("NS3:200", "NS3", "C-Terminal Tail"),
    ("NS3:250", "NS3", "Oxyanion Loop"),
]

#: a distinct base θ per position, ordered so the six positions have six distinct
#: values (a clean strict ranking).
_BASE_THETA: dict[str, float] = {
    "NS2B:-1": 0.05,
    "NS3:51": 0.42,
    "NS3:75": -0.30,
    "NS3:135": 0.15,
    "NS3:200": -0.18,
    "NS3:250": 0.60,
}
_RUNS = [1, 2, 3]


def _theta_for(serotype: str, canon: str, run: int) -> float:
    base = _BASE_THETA[canon]
    if serotype == "DENV1" or serotype == "DENV4":
        # strictly order-preserving across runs → perfect concordance
        return round(base + 0.001 * (run - 1), 6)
    if serotype == "DENV2":
        return round(2.0 * base + 0.001 * (run - 1), 6)
    # DENV3: runs 1 & 2 order-preserving, run 3 negates (reverses the order)
    if run == 3:
        return round(-base, 6)
    return round(base + 0.001 * (run - 1), 6)


def make_replicate_table() -> pd.DataFrame:
    """An S0-shaped replicate table with per-run θ for the default scenario."""
    rows = []
    for serotype in SEROTYPES:
        for canon, chain, domain in _POSITIONS:
            for run in _RUNS:
                # DENV4 loses NS3:250 in the third run (incomplete case)
                if serotype == "DENV4" and canon == "NS3:250" and run == 3:
                    continue
                rows.append(
                    {
                        "serotype": serotype,
                        "replicate": f"{run}_run",
                        "replicate_index": run,
                        "canon_label": canon,
                        "r": _theta_for(serotype, canon, run),
                        "abs_r": abs(_theta_for(serotype, canon, run)),
                        "domain_label": domain,
                        "chain": chain,
                    }
                )
    return pd.DataFrame(rows)


def make_replicate_inventory() -> pd.DataFrame:
    """An S1A-shaped replicate inventory matching :func:`make_replicate_table`."""
    rows = []
    for serotype in SEROTYPES:
        for canon, _chain, _domain in _POSITIONS:
            incomplete = serotype == "DENV4" and canon == "NS3:250"
            n_rep = 2 if incomplete else 3
            rows.append(
                {
                    "serotype": serotype,
                    "canon_label": canon,
                    "n_replicates": n_rep,
                    "replicates": [f"{r}_run" for r in _RUNS[:n_rep]],
                    "replicate_indices": _RUNS[:n_rep],
                    "available": True,
                    "in_all_replicates": not incomplete,
                }
            )
    return pd.DataFrame(rows)


def make_empty_inventory() -> pd.DataFrame:
    """The design's blocked state: an inventory with zero replicates everywhere."""
    rows = []
    for serotype in SEROTYPES:
        for canon, _chain, _domain in _POSITIONS:
            rows.append(
                {
                    "serotype": serotype,
                    "canon_label": canon,
                    "n_replicates": 0,
                    "replicates": [],
                    "replicate_indices": [],
                    "available": False,
                    "in_all_replicates": False,
                }
            )
    return pd.DataFrame(rows)


def make_licensed_inventory(n_replicates: int = 5) -> pd.DataFrame:
    """An inventory whose replicate count licenses residue-scale claims (K >= 5)."""
    inv = make_replicate_inventory()
    inv["n_replicates"] = n_replicates
    inv["in_all_replicates"] = True
    return inv


def write_inputs(
    directory: Path,
    inventory: pd.DataFrame | None = None,
    table: pd.DataFrame | None = None,
) -> tuple[Path, Path | None]:
    """Write inventory (+ optional table) parquet files; return their paths.

    When ``table`` is ``None`` no replicate table is written, reproducing the
    design's blocked state on disk.
    """
    directory.mkdir(parents=True, exist_ok=True)
    inv = inventory if inventory is not None else make_replicate_inventory()
    inv_path = directory / "replicate_inventory.parquet"
    inv.to_parquet(inv_path, index=False)
    if table is None:
        return inv_path, None
    tbl_path = directory / "replicate_table.parquet"
    table.to_parquet(tbl_path, index=False)
    return inv_path, tbl_path
