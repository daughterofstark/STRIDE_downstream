# stride-dengue-analysis

[![CI](https://github.com/OWNER/stride-dengue-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/OWNER/stride-dengue-analysis/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A **downstream analysis framework** for STRIDE reproducibility outputs on the
four dengue serotypes (DENV1–4). STRIDE (the upstream method that produces the
replicate correlations, profiles, and mechanism summaries) is treated as a black
box — this repository never modifies it and never recomputes anything it
produced.

> **Before publishing:** replace `OWNER` in the badge and clone URLs (and in
> `pyproject.toml`) with your GitHub org/user, and confirm the `LICENSE` holder.

> **Status: Stages S0, S1A, S1B, S2, S3, S4, S5 and S6 implemented.**
> **S0** = ingestion, validation, and canonical data-layer construction (the two
> canonical tables every later stage consumes). **S1A** = the reusable *biological
> data layer* built on S0: canonical residue objects, domain summaries, replicate
> availability, and cross-serotype conservation mapping. **S1B** = the *biological
> annotation layer* built on S1A: deterministic per-residue, per-domain,
> per-hierarchy, and per-serotype structural annotations. **S2** = the
> *per-serotype reduction layer* built on S0 + S1B: the achieved-resolution
> census, residue-scale landscape, domain-scale reproducibility, the
> signed/significant screen, and a per-serotype scorecard — all reported over a
> **ρ\* band** (the gate is uncalibrated) with `licensed`/`exploratory` tier
> labels and **no cross-serotype tests**. **S3** = the *hierarchy reduction layer*
> built on the S0 STRIDE table: the ρ-vs-scale curve per locus, the
> domain−residue reproducibility gap, the monotonicity (upward-closure) audit,
> and the chain-level contrast (NS2B vs NS3) — a sibling reduction to S2,
> ρ\*-independent, with the same tier labels and **no cross-serotype tests**.
> **S4** = the *uncertainty layer* built on the S0 STRIDE table: the per-domain
> variance-component budget (τ² vs σ̄²) with the replicate-vs-sampling regime
> diagnostic, the per-residue replicate-disagreement map (positions ranked by τ²),
> the CI-based significance screen with Benjamini–Hochberg FDR control within each
> serotype, and the β_se-weighted effect summary per domain — a sibling reduction
> to S2/S3, ρ\*-independent, with the same tier labels and **no cross-serotype
> tests**.
> **S5** = the *cross-serotype layer* (n = 4) built on the S0 STRIDE table + the
> S1A conservation table: the conservation of reproducibility across shared
> positions (all/majority/some/none, with serotype-divergent and Catalytic-Triad
> flags), the direction concordance of shared signed positions
> (agree/majority/conflict), the tidy-long ρ(domain × serotype) matrix over the
> NS3 domains + NS2B (catalytic domains flagged), and a per-serotype scorecard.
> Serotype is the unit of replication: per-serotype values are aggregated first,
> then compared, and results are **descriptive** (counts / effect sizes), not
> p-values across residues; the domain × serotype matrix is `licensed`, the
> residue-scale products `exploratory`.
> **S6** = the *replicate layer* — the only stage that reads the replicate-specific
> inputs (the S1A `replicate_inventory` and, when it exists, the S0
> `replicate_table` of per-run θ): the per-serotype replicate-regime ledger
> (replicate count K, completeness, residue-claim licensing at K ≥ 5, per-run-effect
> availability), the descriptive across-run per-run-θ spread, the per-serotype rank
> concordance across runs (Kendall's *W* / mean pairwise Spearman, design §3.1), and
> an explicit **blocked-analysis ledger**. It does **not** re-implement the τ²-based
> replicate-disagreement mapping or τ²/σ̄² regime diagnostic (those read the
> aggregate and are S4's); it computes the per-run concordance / spread **only when
> the per-run correlation CSVs were supplied** and otherwise records them blocked,
> and it records **leave-one-replicate-out** stability as permanently blocked
> (it needs a STRIDE re-run) — the design's blocked subset is documented, never
> approximated.
> None of these stages contain **statistics, ranking, mechanism inference, or
> figures** beyond S2/S3/S4's structural reduction, S5's descriptive
> cross-serotype tallies, and S6's replicate-axis diagnostics — publication
> figures and calibrated inference belong to later stages (S7+), which are
> intentionally not implemented here. See
> [`docs/s1a.md`](docs/s1a.md), [`docs/s1b.md`](docs/s1b.md),
> [`docs/s2.md`](docs/s2.md), [`docs/s3.md`](docs/s3.md),
> [`docs/s4.md`](docs/s4.md), [`docs/s5.md`](docs/s5.md), and
> [`docs/s6.md`](docs/s6.md) for the tables and their consumers.

This is a **framework**, not a dataset: real STRIDE outputs are user-supplied
inputs (git-ignored), and everything under `outputs/` is regenerated by running
the pipeline. The repository works immediately after cloning against the
committed synthetic example.

---

## Quick start

```bash
git clone https://github.com/OWNER/stride-dengue-analysis
cd stride-dengue-analysis
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# run S0 against the committed synthetic example
stride-s0 --data-root examples/small_synthetic_dataset --output-dir outputs

# ... or against your own data (see data/README.md for the expected layout)
stride-s0 --data-root data --output-dir outputs

# then build the S1A biological data layer on top of the S0 tables
stride-s1a --input-dir outputs --output-dir outputs_s1a

# then build the S1B biological annotation layer on top of the S1A tables
stride-s1b --input-dir outputs_s1a --output-dir outputs_s1b

# then build the S2 per-serotype reduction layer on top of the S0 + S1B tables
stride-s2 --stride-input-dir outputs --annotation-input-dir outputs_s1b --output-dir outputs_s2

# then build the S3 hierarchy reduction layer on top of the S0 STRIDE table
stride-s3 --input-dir outputs --output-dir outputs_s3

# then build the S4 uncertainty layer on top of the S0 STRIDE table
stride-s4 --input-dir outputs --output-dir outputs_s4

# then build the S5 cross-serotype layer on top of the S0 + S1A tables
stride-s5 --input-dir outputs --conservation-input-dir outputs_s1a --output-dir outputs_s5

# then build the S6 replicate layer on top of the S1A inventory (+ optional S0 replicate table)
stride-s6 --input-dir outputs --inventory-input-dir outputs_s1a --output-dir outputs_s6
```

Programmatic use:

```python
from stride_analysis import run_s0, build_tables

# build both canonical tables and write artifacts
replicate_table, stride_table, report = run_s0("data", "outputs")

# or build in memory without writing
replicate_table, stride_table, report = build_tables("data")
assert report.all_passed
```

---

## Three data levels

| Level | What | Committed? | Read by S0? |
|---|---|---|---|
| **0 — Raw MD inputs** | `xtc`, `gro`, `tpr`, … | no | **never** |
| **1 — Replicate observations** | `*_correlations_v5.csv` (per replicate) | no | yes — *primary* |
| **2 — STRIDE summaries** | `*_profile.csv`, `*_mechanism.json` | no | yes — *derived* |

S0 ingests **both** Level 1 and Level 2 and keeps them as **two separate
canonical tables** — replicate observations are never collapsed into the
summaries. See [`data/README.md`](data/README.md) for the expected directory
layout.

---

## The two canonical tables

### Replicate table (Level 1)

One row per **(serotype, replicate, residue)**, preserving the original per-run
quantities from `*_correlations_v5.csv`.

- Unique key: `(serotype, replicate, canon_label)`
- Written to: `outputs/replicate_table.{parquet,csv}`

### STRIDE table (Level 2)

One row per **(serotype, canon_label, scale_level)**, merging the profile, the
explicit parsed hierarchy (`h_complex … h_residue`), and the mechanism payload
(attached to gated rows only), plus provenance.

- Unique key: `(serotype, canon_label, scale_level)`
- Written to: `outputs/stride_table.{parquet,csv}`

Both tables plus a machine-readable `schema_report.json` and a human-readable
`validation_report.md` are written to the output directory.

---

## Architecture

The code is organised as **reusable subpackages**, not a monolithic pipeline.
Later stages import these directly.

```
src/stride_analysis/
├── models/        # frozen schema, pydantic models, typed error hierarchy
├── io/            # dataset discovery + raw loaders
├── validation/    # Level-1 + Level-2 schema and consistency checks, hierarchy parsing
├── canonical/     # the two canonical-table builders (kept separate)
├── reporting.py   # schema_report.json + validation_report.md writers
├── s0.py          # S0 orchestration (thin composition of the above)
└── __main__.py    # CLI
```

Design principles:

- **Fail loudly.** Every defect raises a typed `StrideAnalysisError` subclass
  (`DiscoveryError`, `SchemaError`, `HierarchyError`, `ConsistencyError`) with an
  actionable message. Nothing is silently dropped or coerced.
- **Frozen schema.** Contracts live in `models/schema.py`; validators assert
  against them. The Level-1 replicate schema is a *required core* + tolerated
  known-optional columns (STRIDE appends columns across milestones).
- **Separation of levels.** Replicate and summary tables are distinct.
- **No absolute paths, no local assumptions.** Discovery resolves whatever
  layout it finds under the data root you pass.

---

## Repository layout

```
stride-dengue-analysis/
├── src/
│   ├── stride_analysis/     # S0: the canonical data layer (installable package)
│   ├── stride_s1a/          # S1A: the reusable biological data layer
│   ├── stride_s1b/          # S1B: the biological annotation layer
│   ├── stride_s2/           # S2: the per-serotype reduction layer
│   └── stride_s3/           # S3: the hierarchy reduction layer
│   └── stride_s4/           # S4: the uncertainty layer
│   └── stride_s5/           # S5: the cross-serotype layer (n=4)
│   └── stride_s6/           # S6: the replicate layer
├── tests/                   # unit + integration tests (synthetic fixtures only)
│   ├── unit/
│   ├── integration/
│   ├── s1a/                 # S1A tests
│   ├── s1b/                 # S1B tests
│   ├── s2/                  # S2 tests
│   └── s3/                  # S3 tests
│   └── s4/                  # S4 tests
│   └── s5/                  # S5 tests
│   └── s6/                  # S6 tests
├── examples/
│   └── small_synthetic_dataset/   # a tiny, valid dataset (committed)
├── docs/                    # architecture, data model, usage, s1a, s1b, s2, s3, s4, s5, s6
├── data/                    # your datasets (git-ignored; see data/README.md)
├── outputs/                 # generated artifacts (git-ignored)
├── notebooks/               # optional exploratory work (not implementation)
├── pyproject.toml
├── LICENSE                  # placeholder — confirm before publishing
├── CONTRIBUTING.md
└── README.md
```

---

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — subpackages and boundaries
- [`docs/data_model.md`](docs/data_model.md) — the two canonical tables, in detail
- [`docs/usage.md`](docs/usage.md) — CLI, Python API, partial datasets
- [`docs/s1a.md`](docs/s1a.md) — the S1A biological data layer and its consumers
- [`docs/s1b.md`](docs/s1b.md) — the S1B biological annotation layer and its consumers
- [`docs/s2.md`](docs/s2.md) — the S2 per-serotype reduction layer and its consumers
- [`docs/s3.md`](docs/s3.md) — the S3 hierarchy reduction layer and its consumers
- [`docs/s4.md`](docs/s4.md) — the S4 uncertainty layer and its consumers
- [`docs/s5.md`](docs/s5.md) — the S5 cross-serotype layer and its consumers
- [`docs/s6.md`](docs/s6.md) — the S6 replicate layer and its consumers
- [`data/README.md`](data/README.md) — expected input layout and data levels

---

## Developer workflow

```bash
pytest                          # run the suite (no real data required)
pytest --cov=stride_analysis    # with coverage
ruff check src tests            # lint
mypy                            # type-check
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

---

## Scope boundary

S0 deliberately excludes — these are later stages:

- biological interpretation of ρ, β, coherence, or direction;
- any statistics, ranking, or clustering;
- cross-serotype comparison or conservation analysis;
- figures, tables, or manuscript text.

S0's single job is to hand later stages two trustworthy, validated tables.

## License

MIT (placeholder — see [`LICENSE`](LICENSE)).
