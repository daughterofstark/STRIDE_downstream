# stride-dengue-analysis

A deterministic, fully-typed downstream-analysis framework that reduces
[STRIDE](https://en.wikipedia.org/wiki/Molecular_dynamics) reproducibility outputs
for the four dengue serotypes (DENV1–4) into publication-ready figures and
manuscript tables — **without** making any calibrated scientific claim the upstream
data does not support.

> **Status:** Stages **S0, S1A, S1B, S2, S3, S4, S5, S6, S7** are implemented,
> tested, and verified. The pipeline is **feature-complete** end to end. The
> project version is `0.2.0`; all stage packages are import-versioned at `0.1.0`
> except the S0 package (`stride_analysis`), which carries the project version
> `0.2.0`.

---

## Project overview

STRIDE emits, per dengue serotype, a **reproducibility profile** of a
molecular-dynamics effect measured across a spatial hierarchy (residue → secondary
structure → motif → domain → chain → protein → complex) for the NS2B–NS3 protease,
together with a gated **mechanism** report. The central quantity is **ρ** (`rho`) ∈
[0, 1], the fraction of an effect's variance at a region/scale that reproduces
across replicate runs.

`stride-dengue-analysis` is the **downstream** framework: it ingests those STRIDE
outputs once (Stage S0), builds a canonical, typed data layer, and then reduces it
through a chain of small, single-responsibility stages into per-serotype summaries,
hierarchy diagnostics, uncertainty budgets, cross-serotype conservation, a
replicate-axis layer, and finally the publication figures and manuscript tables.
Every stage is **deterministic** and **reporting-honest**: it re-derives nothing it
cannot support and documents — rather than approximates — anything that is blocked
by missing upstream information.

## Scientific motivation

The upstream STRIDE gate is **uncalibrated**. In every input file
`calibrated = false`, the gate threshold is a *provisional* constant
`rho_star = 0.5`, and 100 % of emitted mechanisms are flagged `gate_uncertain`.
Two consequences shape the entire design:

1. **ρ is a continuous reproducibility score, not a pass/fail verdict.** The
   framework never asserts a calibrated resolution claim. Results are reported over
   a ρ\* *band* (S2) or at the provisional threshold with an explicit
   `licensed` / `exploratory` **tier** label, never as a hard gate.
2. **Serotype is the unit of biological replication (n = 4).** Each serotype is one
   biological system summarised from **K = 3** replicate runs. Cross-serotype
   statistics (S5) aggregate each serotype to one value *first* and then compare —
   the ~200 residues per serotype are **not** independent samples of "dengue".
   Residue-scale products are labelled `exploratory` (STRIDE licenses residue-scale
   claims only at K ≥ 5); domain-scale products are `licensed` at K = 3.

The framework's job is to extract every defensible, reproducible signal from the
STRIDE outputs while making the uncalibrated status impossible to overlook.

---

## Repository architecture

Every stage is an installable Python package that follows the **same** internal
shape, so the whole repository reads like one idea repeated at increasing analytical
depth:

| Layer | Responsibility |
| --- | --- |
| `models/` | `schema.py` — the **single source of truth** (column names, constants, filenames); `errors.py` — a typed exception hierarchy; report dataclasses. |
| `io/` | Loaders that read prior-stage tables and **assert the columns each stage depends on**; deterministic writers (Parquet/CSV/Markdown/JSON/SVG). |
| `build/` | Pure, deterministic **builders**: one module per output table. They select, join, order, and round — they never mutate inputs and never reach outside their declared sources. |
| `validation/` | **Structural** checks only (completeness, columns, filenames, provenance). Never validates a scientific conclusion. |
| `plotting/` | *(S7 only)* dependency-free SVG rendering, kept out of the builders. |
| `sN.py` | A **thin orchestrator**: `build_sN` (load → build → validate, in memory) and `run_sN` (additionally write artifacts + on-disk validation). |
| `__main__.py` | The stage CLI (`stride-sN`). |

Design rules that hold everywhere: **fail loudly** (a missing/malformed input raises
a typed error, never a silent coercion), **deterministic builders** (sorted rows,
fixed rounding, no wall-clock timestamps), **full typing** (mypy `disallow_untyped_defs`),
and **NumPy-style docstrings**.

---

## Pipeline overview (S0–S7)

```
                 STRIDE outputs (profile.csv + mechanism.json + per-run CSVs)
                                        │
                                        ▼
   S0   ─ ingest & validate → canonical Level-1 / Level-2 tables
    ↓
   S1A  ─ reusable biological data layer (canonical residues, domains, conservation, replicate inventory)
    ↓
   S1B  ─ biological annotation layer (per-residue / -domain / -hierarchy / -serotype annotations)
    ↓
   S2   ─ per-serotype reduction (landscape, resolution census, domain ρ, signed screen, scorecard)
    ↓
   S3   ─ hierarchy reduction (ρ-vs-scale curves, Δρ gap, monotonicity audit, chain contrast)
    ↓
   S4   ─ uncertainty layer (τ²/σ̄² budgets, CI+FDR significance, β_se-weighted effects)
    ↓
   S5   ─ cross-serotype layer, n = 4 (conservation, direction concordance, domain×serotype matrix)
    ↓
   S6   ─ replicate layer (per-run regime, effect spread, rank concordance, blocked-analysis ledger)
    ↓
   S7   ─ reporting layer (publication figures F1–F8, manuscript tables T1–T5)
```

Compact form:

```
S0
 ↓
S1A
 ↓
S1B
 ↓
S2
 ↓
S3
 ↓
S4
 ↓
S5
 ↓
S6
 ↓
S7
```

S2, S3, S4, S5 and S6 are **sibling reductions** — each reads a canonical layer and
produces an independent family of tables; S7 is the only stage that fans in across
S2–S6.

---

## Stage-by-stage summary

- **S0 — Ingest & validate** (`stride_analysis`). Discovers a dataset, reads the raw
  STRIDE `profile.csv` / `mechanism.json` and the per-run correlation CSVs, and emits
  two canonical tables: the **Level-2** `stride_table` (the flattened profile ⋈
  mechanism, one row per locus × scale) and the **Level-1** `replicate_table` (per-run
  observations). Asserts schema, profile↔mechanism consistency, and the
  `gate_uncertain` flags; fails loudly on drift.
- **S1A — Biological data layer** (`stride_s1a`). Canonical residue objects, domain
  summaries, the replicate-availability inventory, and the shared-position
  conservation index — the reusable substrate later stages join against.
- **S1B — Biological annotation layer** (`stride_s1b`). Deterministic per-residue,
  per-domain, per-hierarchy, and per-serotype structural annotations built on S1A.
- **S2 — Per-serotype reduction** (`stride_s2`). Reproducibility landscape,
  achieved-resolution census, domain-scale reproducibility, the signed/significant
  screen, and a per-serotype scorecard — all reported over a **ρ\* band** with
  tier labels and **no** cross-serotype tests.
- **S3 — Hierarchy reduction** (`stride_s3`). The ρ-vs-scale curve per locus, the
  domain−residue reproducibility gap (Δρ), the monotonicity / upward-closure audit,
  and the NS2B-vs-NS3 chain contrast. ρ\*-independent.
- **S4 — Uncertainty layer** (`stride_s4`). The per-domain τ² vs σ̄² variance-component
  budget with the replicate-vs-sampling regime diagnostic, the per-residue
  replicate-disagreement map (ranked by τ²), the CI-based significance screen with
  Benjamini–Hochberg FDR within each serotype, and β_se-weighted effect summaries.
- **S5 — Cross-serotype layer, n = 4** (`stride_s5`). Conservation of reproducibility
  across shared positions, direction concordance of shared signed positions, the
  ρ(domain × serotype) matrix, and a cross-serotype scorecard — **descriptive**,
  serotype as the replication unit.
- **S6 — Replicate layer** (`stride_s6`). The only stage that reads the
  replicate-specific inputs: a per-serotype replicate-regime ledger, a descriptive
  per-run effect spread, per-serotype rank concordance (Kendall's *W* / mean pairwise
  Spearman), and an explicit **blocked-analysis ledger**. It does not duplicate S4's
  τ²-based aggregate products, and records leave-one-replicate-out stability as
  permanently blocked (it needs a STRIDE re-run).
- **S7 — Reporting layer** (`stride_s7`). Assembles the design's eight publication
  figures (**F1–F8**) and five manuscript tables (**T1–T5**) from the S2–S6 outputs.
  **Reporting only** — no new statistics, no inference, never reads raw STRIDE.
  Figures are emitted as deterministic, dependency-free **SVG** (plus prepared-data
  CSV/Parquet); tables as CSV/Parquet/Markdown; with a provenance-stamped summary
  JSON. The replicate layer (S6) feeds no design figure/table, so its blocked-analysis
  ledger is surfaced under the report's `limitations` rather than fabricated into a
  figure.

---

## Repository tree

```
stride-dengue-analysis/
├── src/
│   ├── stride_analysis/     # S0: ingestion + canonical data layer
│   │   ├── canonical/       #   stride_table + replicate_table builders
│   │   ├── io/              #   dataset discovery + raw loaders
│   │   ├── models/          #   schema, errors
│   │   ├── validation/      #   cross-level, hierarchy, replicate, summary
│   │   ├── reporting.py     #   schema_report + validation_report writers
│   │   └── s0.py
│   ├── stride_s1a/          # S1A: biological data layer
│   ├── stride_s1b/          # S1B: biological annotation layer
│   ├── stride_s2/           # S2: per-serotype reduction layer
│   ├── stride_s3/           # S3: hierarchy reduction layer
│   ├── stride_s4/           # S4: uncertainty layer
│   ├── stride_s5/           # S5: cross-serotype layer (n = 4)
│   ├── stride_s6/           # S6: replicate layer
│   └── stride_s7/           # S7: reporting layer (figures & tables)
│       ├── build/           #   figure-data prep + table assembly
│       └── plotting/        #   deterministic SVG rendering
│           (each stage package = models/ io/ build/ validation/ sN.py __main__.py)
├── tests/                   # synthetic fixtures only
│   ├── unit/ integration/   #   S0 unit + integration
│   └── s1a/ s1b/ s2/ s3/ s4/ s5/ s6/ s7/
├── examples/
│   └── small_synthetic_dataset/   # a tiny, valid dataset (committed)
├── docs/                    # architecture, data_model, usage, s1a … s7
├── data/                    # your datasets (git-ignored; see data/README.md)
├── .github/workflows/ci.yml
├── pyproject.toml
├── CONTRIBUTING.md
├── LICENSE                  # MIT (placeholder holder — confirm before publishing)
├── STRIDE_MASTER_HANDOFF.md
└── README.md
```

---

## Installation

Requires **Python ≥ 3.10**. Runtime dependencies are `pandas`, `pyarrow`, and
`pydantic`; no plotting or scientific-stack extras are needed.

```bash
git clone https://github.com/OWNER/stride-dengue-analysis
cd stride-dengue-analysis
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"      # editable install with the dev toolchain
```

The install exposes **nine** console scripts: `stride-s0`, `stride-s1a`,
`stride-s1b`, `stride-s2`, `stride-s3`, `stride-s4`, `stride-s5`, `stride-s6`,
`stride-s7`.

## Quick start

Run the whole pipeline against the committed synthetic example:

```bash
stride-s0  --data-root examples/small_synthetic_dataset --output-dir outputs
stride-s1a --input-dir outputs           --output-dir outputs_s1a
stride-s1b --input-dir outputs_s1a       --output-dir outputs_s1b
stride-s2  --stride-input-dir outputs --annotation-input-dir outputs_s1b --output-dir outputs_s2
stride-s3  --input-dir outputs           --output-dir outputs_s3
stride-s4  --input-dir outputs           --output-dir outputs_s4
stride-s5  --input-dir outputs --conservation-input-dir outputs_s1a --output-dir outputs_s5
stride-s6  --input-dir outputs --inventory-input-dir outputs_s1a    --output-dir outputs_s6
stride-s7  --s2-input-dir outputs_s2 --s3-input-dir outputs_s3 \
           --s4-input-dir outputs_s4 --s5-input-dir outputs_s5 \
           --s6-input-dir outputs_s6 --output-dir outputs_s7
```

Point `stride-s0 --data-root` at your own dataset (see `data/README.md` for the
expected layout) to run against real data.

## Running every stage

Each stage reads the output directory of the stage(s) it depends on and writes to
its own. The dependency edges are:

- **S0** ← the raw dataset (`--data-root`)
- **S1A** ← S0
- **S1B** ← S1A
- **S2** ← S0 (STRIDE table) + S1B (annotations)
- **S3** ← S0
- **S4** ← S0
- **S5** ← S0 (STRIDE table) + S1A (conservation table)
- **S6** ← S1A (replicate inventory) + S0 (optional replicate table)
- **S7** ← S2 + S3 + S4 + S5 + S6

Every stage's orchestrator also has a Python API — `build_sN(...)` (returns the
tables + report without writing) and `run_sN(...)` (writes artifacts). See the
CLI reference below and `docs/usage.md`.

---

## CLI reference

All flags have sensible defaults; only `stride-s0 --data-root` typically needs
setting. Every stage additionally accepts explicit `--<table>` path overrides for
its inputs (see `--help`).

| Command | Key flags (default) |
| --- | --- |
| `stride-s0` | `--data-root` (`data`), `--output-dir` (`outputs`), `--no-require-replicates`, `--no-require-summaries`, `--allow-unequal-replicates`, `--strict-cross-level` |
| `stride-s1a` | `--input-dir` (`outputs`), `--output-dir` (`outputs_s1a`) |
| `stride-s1b` | `--input-dir` (`outputs_s1a`), `--output-dir` (`outputs_s1b`) |
| `stride-s2` | `--stride-input-dir` (`outputs`), `--annotation-input-dir` (`outputs_s1b`), `--rho-star` (band `0.50–0.90` step `0.05`), `--output-dir` (`outputs_s2`) |
| `stride-s3` | `--input-dir` (`outputs`), `--output-dir` (`outputs_s3`) |
| `stride-s4` | `--input-dir` (`outputs`), `--output-dir` (`outputs_s4`) |
| `stride-s5` | `--input-dir` (`outputs`), `--conservation-input-dir` (`outputs_s1a`), `--rho-star` (`0.5`), `--output-dir` (`outputs_s5`) |
| `stride-s6` | `--input-dir` (`outputs`), `--inventory-input-dir` (`outputs_s1a`), `--no-replicate-table`, `--output-dir` (`outputs_s6`) |
| `stride-s7` | `--s2-input-dir` (`outputs_s2`), `--s3-input-dir` (`outputs_s3`), `--s4-input-dir` (`outputs_s4`), `--s5-input-dir` (`outputs_s5`), `--s6-input-dir` (`outputs_s6`), `--output-dir` (`outputs_s7`) |

Each command exits non-zero and prints an actionable `… FAILED: <ErrorType>: <msg>`
on a bad input, and prints a one-line success summary otherwise.

## Output directories

| Stage | Default dir | Artifacts |
| --- | --- | --- |
| S0 | `outputs/` | `stride_table.{parquet,csv}`, `replicate_table.{parquet,csv}`, `schema_report.json`, `validation_report.md` |
| S1A | `outputs_s1a/` | `canonical_residues`, `domain_table`, `replicate_inventory`, `conservation_table` (`.parquet`), `dataset_summary.json` |
| S1B | `outputs_s1b/` | `residue_annotation`, `domain_annotation`, `hierarchy_annotation`, `serotype_annotation` (`.parquet`), `annotation_summary.json` |
| S2 | `outputs_s2/` | `residue_landscape`, `resolution_census`, `domain_reproducibility`, `signed_screen`, `serotype_summary` (`.parquet`), `reduction_summary.json` |
| S3 | `outputs_s3/` | `scale_curve`, `resolution_gap`, `monotonicity_audit`, `chain_contrast` (`.parquet`), `hierarchy_summary.json` |
| S4 | `outputs_s4/` | `variance_budget`, `residue_variance`, `significance_screen`, `domain_effect_summary` (`.parquet`), `uncertainty_summary.json` |
| S5 | `outputs_s5/` | `position_conservation`, `direction_concordance`, `domain_serotype_matrix`, `cross_serotype_scorecard` (`.parquet`), `conservation_summary.json` |
| S6 | `outputs_s6/` | `replicate_regime`, `replicate_effect_spread`, `replicate_concordance`, `replicate_blocked_analyses` (`.parquet`), `replicate_summary.json` |
| S7 | `outputs_s7/` | **F1–F8** (`.svg` + `.csv` + `.parquet` each), **T1–T5** (`.csv` + `.parquet` + `.md` each), `artifact_manifest.parquet`, `report_summary.json` |

All output directories are git-ignored; they are always regenerable from the inputs.

---

## Repository invariants

- **Frozen stages.** S0–S7 are complete. A stage is only ever modified for a genuine
  correctness bug (smallest possible fix + full re-verification); new analysis goes
  in a new stage.
- **Reporting honesty.** No stage fabricates or approximates a quantity it cannot
  support from its declared inputs. Blocked analyses (e.g. S6's per-run concordance
  when the per-run CSVs are absent, or leave-one-replicate-out stability) are
  **documented**, never guessed.
- **No raw-data or generated-output commits.** Only the committed synthetic example
  under `examples/` is tracked; real data under `data/` and every `outputs*/`
  directory are ignored.
- **Two data levels stay separate.** Level-1 replicate observations and Level-2
  STRIDE summaries are distinct canonical tables and are never merged.
- **Uncalibrated gate is never hidden.** Every reduction carries `calibrated = false`,
  the provisional ρ\*, and `licensed` / `exploratory` tier labels; no calibrated
  pass/fail claim is made anywhere.

## Determinism and reproducibility

Running any stage twice on the same input produces **byte-identical** output.
This is achieved by:

- sorting every emitted table by stable keys and rounding floats to a fixed
  precision;
- writing CSVs with a fixed line terminator and formatting, and Parquet with a
  fixed column order;
- embedding **no** wall-clock timestamp in any artifact (S7's provenance
  intentionally omits the design's "date" field for this reason);
- rendering S7 figures as **hand-built SVG** rather than via a plotting library —
  raster back-ends (matplotlib &c.) are version-sensitive and not byte-reproducible,
  and would add a heavyweight dependency, so PNG/PDF are intentionally not emitted.
  SVG is a vector, publication-quality, dependency-free format that converts to
  raster downstream if needed.

Each stage stamps a provenance header (input SHA-256 digests, ρ\*, `calibrated = false`,
the K = 3 note) into its summary JSON.

## Testing

The suite is **534 tests**, all using **synthetic fixtures only** (no real data, no
network, no reliance on another stage's tests):

```bash
pytest                                   # 534 passed
pytest tests/s7 -q                       # a single stage's suite
pytest --cov=stride_analysis --cov=stride_s7 --cov-report=term-missing
```

Per suite: `unit` 58, `integration` 10, `s1a` 44, `s1b` 72, `s2` 65, `s3` 50,
`s4` 61, `s5` 70, `s6` 55, `s7` 49. Each stage's tests cover every builder, the
orchestrator, structural validation, deterministic output, the CLI, and both
empty-input and missing-input handling.

## Documentation

The `docs/` directory holds the authoritative per-topic references:

- [`docs/architecture.md`](docs/architecture.md) — the shared package architecture and the stage graph
- [`docs/data_model.md`](docs/data_model.md) — the canonical tables and the hierarchy grammar
- [`docs/usage.md`](docs/usage.md) — end-to-end usage
- [`docs/s1a.md`](docs/s1a.md), [`docs/s1b.md`](docs/s1b.md), [`docs/s2.md`](docs/s2.md), [`docs/s3.md`](docs/s3.md), [`docs/s4.md`](docs/s4.md), [`docs/s5.md`](docs/s5.md), [`docs/s6.md`](docs/s6.md), [`docs/s7.md`](docs/s7.md) — each stage's tables and consumers

`STRIDE_MASTER_HANDOFF.md` is the exhaustive, self-contained maintenance handoff.

## Development workflow

```bash
pip install -e ".[dev]"
ruff check .        # lint (E, F, I, UP, B, W; line length by formatter)
mypy                # strict typing across src + tests
pytest              # full suite
```

CI (GitHub Actions) runs the same four gates plus the full S0→S7 CLI smoke chain on
Python 3.10, 3.11, and 3.12. Keep changes additive, mirror the existing package
shape, and re-run all four gates before opening a PR.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md). In short: synthetic data only in git, no
generated outputs committed, fail loudly with typed errors, keep the two data levels
separate, put reusable logic in the subpackages and keep orchestrators thin, and
respect the frozen stages.

## Citation

A citation entry will be added on publication. Until then, please cite the
repository URL and the commit hash you used. *(Placeholder — update before public
release.)*

## License

Released under the **MIT License** — see [`LICENSE`](LICENSE). Confirm the copyright
holder before making the repository public.
