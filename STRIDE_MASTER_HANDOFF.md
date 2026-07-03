# STRIDE Downstream Analysis — Master Handoff Document

**Purpose.** This document is a complete, self-contained technical handoff for the
`stride-dengue-analysis` repository (a.k.a. `STRIDE_downstream-main`). It is
written so that a **brand-new Claude conversation with zero prior context** can
pick up development and implement the next stage without any additional
knowledge. Attach this file (and the repository ZIP) to the new conversation.

It documents the exact current state, the architecture every stage must follow,
the invariants that must never be broken, the development and verification
workflow, every pitfall encountered and its resolution, precisely what is done
and what remains, and the exact next stage to build.

> **Scope note.** The repository implements a *reproducibility / data-analysis
> framework* for STRIDE outputs on the four dengue serotypes (DENV1–4). Every
> stage is a deterministic table-reduction pipeline. There is **no wet-lab, no
> networked service, and no ML**. All work is pure Python (pandas + pyarrow +
> pydantic), fully typed, fully tested.

---

## 0. TL;DR orientation

- **Implemented and shipped:** S0, S1A, S1B, S2, S3, S4. All green.
- **Not yet implemented:** S5 (cross-serotype, n=4), S6 (replicate layer — blocked
  on data that is not available), S7 (figures & tables).
- **The next stage to build is S5.**
- **Baseline health:** `ruff check .` ✅, `mypy` ✅ (155 source files), `pytest`
  ✅ **360 passed**, editable install ✅, all six CLIs ✅, CI smoke S0→S4 ✅.
- **Git-tracked files:** 193.
- **Golden rule:** never modify a completed stage except for a genuine,
  explained, re-verified bug fix. Each new stage is purely additive.

---

# 1. Current repository state

## 1.1 What the project is

STRIDE produces, for each dengue serotype, a per-residue *reproducibility profile*
across a 7-level structural hierarchy, plus a *mechanism* file describing gated
(significant) effects. This repository ingests those raw outputs (S0), builds
annotation layers (S1A/S1B), and then performs a series of independent,
deterministic *reduction* analyses (S2/S3/S4, and future S5+) that read the S0
tidy master table and emit reusable parquet tables consumed by a future
figures/tables stage (S7).

Key scientific framing (baked into the code and docstrings):

- The **reproducibility statistic ρ** and the **variance components τ² / σ̄²** are
  produced upstream by STRIDE and are only ever **read**, never recomputed.
- The **gate is uncalibrated**: a provisional threshold ρ\* = 0.5 ships in every
  mechanism file. No stage makes a calibrated pass/fail resolution claim.
- **Two output tiers:** *licensed* (domain-scale and coarser — the claim level at
  K = 3 replicates) and *exploratory* (residue-scale, caveated).
- **Serotype is the unit of biological replication (n = 4).** No cross-serotype
  test may treat residues or frames as independent samples. Cross-serotype work
  is deferred to S5.

## 1.2 Implemented stages (summary table)

| Stage | Package | Version | Input → Output | One-line purpose |
|---|---|---|---|---|
| S0 | `stride_analysis` | 0.2.0 | 8 raw files → typed canonical tables | Ingest, validate schema, build the tidy STRIDE master table |
| S1A | `stride_s1a` | 0.1.0 | S0 tables → canonical residue/domain/conservation tables | Canonical join layer + 199-shared-position index |
| S1B | `stride_s1b` | 0.1.0 | S1A tables → biological annotation tables | Residue/domain/hierarchy/serotype annotation |
| S2 | `stride_s2` | 0.1.0 | S0 tidy + S1B annotation → per-serotype summaries | Reproducibility landscape / census / domain ρ / signed screen over a ρ\* band |
| S3 | `stride_s3` | 0.1.0 | S0 tidy → hierarchy (scale-axis) reduction | ρ-vs-scale curve, Δρ gap, monotonicity audit, chain contrast |
| S4 | `stride_s4` | 0.1.0 | S0 tidy → uncertainty layer | τ²/σ̄² budgets, τ² ranking, CI+FDR significance screen, β_se-weighted summaries |

**Dependency graph (what each stage reads):**

```
raw STRIDE files
   │
   ▼
  S0 ── stride_table.parquet (the tidy master; 45 cols) ─┬──────────────┬───────────────┐
   │                                                     │              │               │
   ├── replicate_table.parquet                           │              │               │
   ▼                                                     │              │               │
  S1A ── canonical_residues / domain_table /             │              │               │
   │      conservation_table / replicate_inventory       │              │               │
   ▼                                                     │              │               │
  S1B ── residue/domain/hierarchy/serotype_annotation ───┘ (S2 only)    │               │
                                                                        │               │
  S2 reads  S0 stride_table + S1B annotation ───────────────────────────┘               │
  S3 reads  S0 stride_table ONLY ─────────────────────────────────────────────┐         │
  S4 reads  S0 stride_table ONLY ─────────────────────────────────────────────┴─────────┘

  S2, S3, S4 are SIBLING reductions: each reads S0 (S2 also reads S1B) independently.
  No reduction stage consumes another reduction stage's outputs.
```

## 1.3 Package layout (every source file)

All packages live under `src/` (src-layout, importable only after editable
install). Every stage-package uses the identical five-part internal layout:
`models/`, `io/`, `build/`, `validation/`, plus a thin orchestrator module and a
`__main__.py` CLI.

```
src/
  stride_analysis/                 # S0 — the exception to the reduction layout (it is the ingester)
    __init__.py  __main__.py  s0.py  reporting.py  _synthetic.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  discovery.py  loaders.py
    canonical/ __init__.py  replicate_table.py  stride_table.py
    validation/ __init__.py  cross_level.py  hierarchy.py  replicate.py  summary.py

  stride_s1a/                      # S1A — canonical join layer
    __init__.py  __main__.py  s1a.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  loaders.py  writers.py
    build/     __init__.py  canonical_residues.py  domain_table.py
               replicate_inventory.py  conservation_table.py
    validation/ __init__.py  checks.py

  stride_s1b/                      # S1B — biological annotation layer
    __init__.py  __main__.py  s1b.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  loaders.py  writers.py
    build/     __init__.py  _classify.py  residue_annotation.py  domain_annotation.py
               hierarchy_annotation.py  serotype_annotation.py
    validation/ __init__.py  checks.py

  stride_s2/                       # S2 — per-serotype reduction (ρ* band)
    __init__.py  __main__.py  s2.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  loaders.py  writers.py
    build/     __init__.py  _frames.py  _screens.py  resolution_census.py
               residue_landscape.py  domain_reproducibility.py  signed_screen.py
               serotype_summary.py
    validation/ __init__.py  checks.py

  stride_s3/                       # S3 — hierarchy (scale-axis) reduction
    __init__.py  __main__.py  s3.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  loaders.py  writers.py
    build/     __init__.py  _curves.py  _frames.py  scale_curve.py
               resolution_gap.py  monotonicity_audit.py  chain_contrast.py
    validation/ __init__.py  checks.py

  stride_s4/                       # S4 — uncertainty layer
    __init__.py  __main__.py  s4.py
    models/    __init__.py  errors.py  schema.py
    io/        __init__.py  loaders.py  writers.py
    build/     __init__.py  _stats.py  _frames.py  variance_budget.py
               residue_variance.py  significance_screen.py  domain_effect_summary.py
    validation/ __init__.py  checks.py
```

**Non-source top-level:** `.github/workflows/ci.yml`, `.gitignore`,
`CONTRIBUTING.md`, `LICENSE`, `README.md`, `pyproject.toml`, `docs/`, `data/`
(git-ignored real data; only `data/README.md` committed), `examples/`
(`small_synthetic_dataset/` — a committed valid dataset), `notebooks/`,
`outputs/` (git-ignored except `.gitkeep`).

## 1.4 Architecture (the layered contract)

Each stage package is composed of five responsibilities. This is the pattern to
replicate for every new stage.

- **`models/schema.py`** — the *single source of truth* for that stage: every
  output column tuple (`*_COLUMNS`), every primary-key tuple (`*_KEY`), every
  constant (thresholds, tier labels, decimals, provisional ρ\*), every input
  required-column tuple, and every output filename (`OUT_*`, `IN_*`). All values
  are module-level `Final`. No logic.
- **`models/errors.py`** — the stage's exception hierarchy. A single base
  (`S{N}Error`), and typically `InputError` (bad/missing input) and
  `ConsistencyError` (an internal invariant failed). S2 also has `ConfigError`.
- **`models/__init__.py`** — the report carriers: a `ValidationCheck` dataclass
  `(name: str, scope: str, passed: bool, detail: str = "")` and an `S{N}Report`
  dataclass with per-table `n_*` counts, `serotypes: list[str]`,
  `provenance: dict`, `facts: dict`, an `.all_passed` property, and an
  `.add(name, scope, passed, detail)` method.
- **`io/loaders.py`** — pure input loading: read each input parquet, assert the
  required columns are present, and provide `file_digest(path)` (SHA-256 hex, or
  `""` if missing) for the provenance header. Never mutates.
- **`io/writers.py`** — write each output parquet (`index=False`) and the summary
  JSON. The JSON payload is a fixed shape (see §1.7).
- **`build/`** — the **pure builders**. One module per output table, each exposing
  a `build_<table>(...) -> pd.DataFrame`. Pure: no IO, no global state, no
  mutation of inputs, deterministic. Shared internal helpers live in
  underscore-prefixed modules (`_frames.py`, `_stats.py`, `_screens.py`,
  `_curves.py`, `_classify.py`) and are also pure. `build/__init__.py` re-exports
  the public builders.
- **`validation/checks.py`** — structural/arithmetic validators. Each takes the
  built tables plus the report, appends a `ValidationCheck`, and raises
  `ConsistencyError` on failure. **No statistical or biological assertions.**
- **`s{n}.py`** (orchestrator) — thin composition only: load → build → validate →
  (optionally) write. Exposes `build_s{n}(...) -> (S{N}Tables, S{N}Report)` (no
  writes) and `run_s{n}(..., output_dir) -> (S{N}Tables, S{N}Report)` (writes
  artifacts). `S{N}Tables` is a `@dataclass` holding the built DataFrames.
- **`__main__.py`** (CLI) — argparse wrapper calling `run_s{n}`. Prints
  `S{N} OK: …` on success (returns 0) or `S{N} FAILED: <ErrorType>: <msg>` to
  stderr (returns 1). Registered as `stride-s{n}` in `pyproject.toml`.
- **`__init__.py`** — the public API: re-exports `run_s{n}`, `build_s{n}`,
  `S{N}Tables`, the individual builders, `S{N}Report`, and the error classes,
  with `__all__` and `__version__`.

## 1.5 Public APIs (exact, per package)

Each `__init__.py` exports (confirmed from the installed packages):

```python
# stride_analysis (S0)  __version__ = "0.2.0"
run_s0, build_tables, discover_dataset, Dataset, SerotypeDataset,
MechanismFile, Report, StrideAnalysisError, DiscoveryError, SchemaError,
HierarchyError, ConsistencyError

# stride_s1a  __version__ = "0.1.0"
run_s1a, build_s1a, S1ATables,
build_canonical_residues, build_domain_table, build_replicate_inventory,
build_conservation_table,
S1AReport, S1AError, InputError, ConsistencyError

# stride_s1b  __version__ = "0.1.0"
run_s1b, build_s1b, S1BTables,
build_residue_annotation, build_domain_annotation, build_hierarchy_annotation,
build_serotype_annotation,
S1BReport, S1BError, InputError, ConsistencyError

# stride_s2  __version__ = "0.1.0"
run_s2, build_s2, S2Tables,
build_resolution_census, build_residue_landscape, build_domain_reproducibility,
build_signed_screen, build_serotype_summary,
S2Report, S2Error, InputError, ConfigError, ConsistencyError

# stride_s3  __version__ = "0.1.0"
run_s3, build_s3, S3Tables,
build_scale_curve, build_resolution_gap, build_monotonicity_audit,
build_chain_contrast,
S3Report, S3Error, InputError, ConsistencyError

# stride_s4  __version__ = "0.1.0"
run_s4, build_s4, S4Tables,
build_variance_budget, build_residue_variance, build_significance_screen,
build_domain_effect_summary,
S4Report, S4Error, InputError, ConsistencyError
```

Orchestrator signatures follow one of two shapes:

```python
build_s{n}(<input paths…>) -> tuple[S{N}Tables, S{N}Report]     # no writes
run_s{n}(<input paths…>, output_dir) -> tuple[S{N}Tables, S{N}Report]  # writes
```

S3/S4 take a single `stride_table_path`; S2 takes stride + annotation inputs;
S1A/S1B take their upstream dirs/paths; S0 takes a data-root.

## 1.6 CLI entry points (exact flags & defaults)

All are registered under `[project.scripts]` and print `S{N} OK`/`S{N} FAILED`.

| CLI | Flags | Output default |
|---|---|---|
| `stride-s0` | `--data-root` (def `data`), `--output-dir` (def `outputs`), `--no-require-replicates`, `--no-require-summaries`, `--allow-unequal-replicates`, `--strict-cross-level` | `outputs/` |
| `stride-s1a` | `--input-dir` (def `outputs`), `--stride-table`, `--replicate-table`, `--output-dir` (def `outputs_s1a`) | `outputs_s1a/` |
| `stride-s1b` | `--input-dir` (def `outputs_s1a`), `--canonical-residues`, `--domain-table`, `--conservation-table`, `--replicate-inventory`, `--output-dir` (def `outputs_s1b`) | `outputs_s1b/` |
| `stride-s2` | `--stride-input-dir` (def `outputs`), `--annotation-input-dir` (def `outputs_s1b`), `--stride-table`, `--residue-annotation`, `--domain-annotation`, `--rho-star`, `--output-dir` (def `outputs_s2`) | `outputs_s2/` |
| `stride-s3` | `--input-dir` (def `outputs`), `--stride-table`, `--output-dir` (def `outputs_s3`) | `outputs_s3/` |
| `stride-s4` | `--input-dir` (def `outputs`), `--stride-table`, `--output-dir` (def `outputs_s4`) | `outputs_s4/` |

Canonical local run order:

```bash
stride-s0  --data-root examples/small_synthetic_dataset --output-dir outputs
stride-s1a --input-dir outputs      --output-dir outputs_s1a
stride-s1b --input-dir outputs_s1a  --output-dir outputs_s1b
stride-s2  --stride-input-dir outputs --annotation-input-dir outputs_s1b --output-dir outputs_s2
stride-s3  --input-dir outputs      --output-dir outputs_s3
stride-s4  --input-dir outputs      --output-dir outputs_s4
```

## 1.7 Outputs (artifacts per stage)

Every reduction stage writes four/five parquet tables + one summary JSON. The
JSON payload shape is uniform:

```json
{
  "stage": "S4",
  "all_checks_passed": true,
  "serotypes": ["DENV1", "..."],
  "provenance": { "calibrated": false, "provisional_rho_star": 0.5, "...": "...",
                  "inputs": { "<name>": { "path": "...", "sha256": "..." } } },
  "n_<table1>": 8, "n_<table2>": 8, "...": 0,
  "facts": { "...": "..." },
  "checks": [ { "name": "...", "scope": "...", "passed": true, "detail": "..." } ]
}
```

| Stage | Parquet outputs | Summary JSON |
|---|---|---|
| S0 | `stride_table.parquet` (+`.csv`), `replicate_table.parquet` (+`.csv`) | `schema_report.json`, `validation_report.md` |
| S1A | `canonical_residues.parquet`, `domain_table.parquet`, `replicate_inventory.parquet`, `conservation_table.parquet` | `dataset_summary.json` |
| S1B | `residue_annotation.parquet`, `domain_annotation.parquet`, `hierarchy_annotation.parquet`, `serotype_annotation.parquet` | `annotation_summary.json` |
| S2 | `resolution_census.parquet`, `residue_landscape.parquet`, `domain_reproducibility.parquet`, `signed_screen.parquet`, `serotype_summary.parquet` | `reduction_summary.json` |
| S3 | `scale_curve.parquet`, `resolution_gap.parquet`, `monotonicity_audit.parquet`, `chain_contrast.parquet` | `hierarchy_summary.json` |
| S4 | `variance_budget.parquet`, `residue_variance.parquet`, `significance_screen.parquet`, `domain_effect_summary.parquet` | `uncertainty_summary.json` |

All output directories are git-ignored (`outputs/`, `outputs_s1a/` … `outputs_s4/`;
future stages must add `outputs_s5/` etc.).

## 1.8 The S0 STRIDE master table (the tidy frame every reduction reads)

`stride_table.parquet` — **45 columns**, one row per (serotype, canon_label,
scale_level); i.e. every locus is present at all 7 hierarchy scales. This is the
"tidy" / "profile" input referenced throughout the design.

```
serotype, canon_label, scale_level, scale_index, locus, region_id,
rho, gated, beta, beta_se, tau2, sigma2_bar, a_signed, coherence, method, status,
h_complex, h_protein, h_chain, h_domain, h_motif, h_secondary_structure, h_residue,
is_gated_scale,
mech_label, mech_direction, mech_beta_signed, mech_beta_ci_lower, mech_beta_ci_upper,
mech_beta_se, mech_coherence, mech_reproducible_magnitude_energy, mech_rho_star,
mech_calibrated, mech_gate_uncertain, mech_status, mech_region_id, mech_n_loci,
profile_source, mechanism_source,
gate_rho_star, gate_alpha, gate_coherence_threshold, mechanism_calibrated,
mechanism_schema_version
```

**Scale grammar** (7 levels, finest→coarsest, `scale_index` 0→6):
`residue(0), secondary_structure(1), motif(2), domain(3), chain(4), protein(5),
complex(6)`.

**Empirically verified facts about the synthetic example** (used to design
fixtures — a new stage can rely on these being the shape of real data too):

- Variance components (`tau2`, `sigma2_bar`, `beta`, `beta_se`) at the **domain
  scale** are **region-constant** per `(serotype, h_chain, h_domain)` (nunique=1),
  because a domain region's aggregate values are shared by its member loci.
- `gate_alpha` = 0.05 everywhere; `gate_rho_star` = 0.5 (provisional);
  `mech_calibrated`/`mechanism_calibrated` = false.
- Signed gated mechanisms carry real `mech_beta_signed`/`mech_beta_se`/CI; **mixed**
  mechanisms carry NaN for those.
- `approx identity`: `rho ≈ beta² / (beta² + tau2 + sigma2_bar)` holds closely but
  not exactly (the estimator applies shrinkage) — so **never reconstruct ρ from
  the variance components**; read `rho` directly.
- Cross-serotype join key is `canon_label` (canonical numbering). In real data:
  199 residue labels are shared by all four serotypes; union is 248. The S1A
  `conservation_table` already encodes this (`n_serotypes`, `serotypes_present`),
  and S1B assigns `conservation_class ∈ {pan_serotype, partial, serotype_unique}`.

## 1.9 Documentation

`docs/` contains one file per stage plus cross-cutting docs. Each stage doc
follows an identical structure (mirror it for new stages):

- `architecture.md` — the layered contract + a paragraph per implemented stage +
  a forward-reference to the next unimplemented stage (currently "S5+").
- `data_model.md` — the S0 canonical data model and the 7-level hierarchy grammar.
- `usage.md` — end-to-end run instructions.
- `s1a.md`, `s1b.md`, `s2.md`, `s3.md`, `s4.md` — per-stage: intro, "what it is
  not", input table, core concepts, each output (schema + key + producer +
  downstream consumers + invariants), validation, usage (bash + Python), and a
  consumption-map quick-reference table.
- `README.md` (top-level) — status line, per-stage description block, quick-start,
  repo tree, docs list.
- `CONTRIBUTING.md` — the additive-stage rules.

When adding a stage you MUST: create `docs/s{n}.md`; add a paragraph +
forward-reference bump in `docs/architecture.md`; and update `README.md`
(status line, description, quick-start, repo tree, docs list).

## 1.10 Tests

`tests/` layout mirrors the packages. Per-stage suites live in `tests/s{n}/` with
`__init__.py`, `conftest.py` (pytest fixtures, usually a `stride_table` fixture),
`fixtures.py` (synthetic-data builders — no dependence on real data or on running
upstream stages), and one `test_*.py` per concern (builders, validation,
pipeline/CLI, and stage-specific helpers).

**Current counts (all passing, total = 360):**

| Suite | Tests |
|---|---|
| `tests/unit` (S0 unit) | 58 |
| `tests/integration` (S0 pipeline + CLI) | 10 |
| `tests/s1a` | 44 |
| `tests/s1b` | 72 |
| `tests/s2` | 65 |
| `tests/s3` | 50 |
| `tests/s4` | 61 |
| **Total** | **360** |

Each per-stage suite covers: pure helpers against hand-computed values; each
builder's schema/keys/values; every validation failure path (assert it raises
`ConsistencyError`); loader errors (missing file, missing columns, unreadable);
orchestration end-to-end; **determinism** (run twice, `assert_frame_equal`);
empty-input tolerance; and the CLI (success rc=0 + "S{N} OK"; failure rc=1 +
"S{N} FAILED").

## 1.11 CI

`.github/workflows/ci.yml` — matrix on Python 3.10 / 3.11 / 3.12. Steps:
`pip install -e ".[dev]"` → `ruff check .` → `mypy` → `pytest --cov=… (all six
packages) --cov-report=term-missing` → smoke-run S0→S4 in sequence on the
committed example (`stride-s0 … /tmp/s0_ci`, then each downstream stage reads the
prior `/tmp/*_ci`). Adding a stage requires appending its `--cov=stride_s{n}` and
a smoke-run step.

## 1.12 Tooling configuration (from `pyproject.toml`)

- **Build:** setuptools, src-layout, `packages.find where=["src"]`.
- **Runtime deps:** `pandas>=2.0`, `pyarrow>=12.0`, `pydantic>=2.5`.
- **Dev deps:** `pytest>=7.4`, `pytest-cov>=4.1`, `mypy>=1.8`, **`ruff>=0.15,<0.16`
  (pinned — see pitfalls §5)**.
- **pytest:** `testpaths=["tests"]`, `addopts="-q"`, `pythonpath=["src"]`.
- **mypy:** `python_version="3.10"`, `files=[all six src pkgs + "tests"]`,
  `ignore_missing_imports=true`, `disallow_untyped_defs=true`,
  `disallow_incomplete_defs=true`, `warn_redundant_casts=true`,
  `warn_unused_ignores=true`, `no_implicit_optional=true`. Override for
  `numpy.*/pandas.*/pyarrow.*`: `follow_imports="skip"`.
- **ruff:** line-length 88, target py310, `select=["E","F","I","UP","B","W"]`,
  `ignore=["E501"]`, **`known-first-party` lists all six packages + `tests`**
  (critical — see pitfalls §5), `per-file-ignores` allows `B011` in tests.
- **coverage:** branch, source = all six packages, omit `_synthetic.py`.

---

# 2. Design philosophy — conventions every future stage MUST follow

## 2.1 Pure builders

Every table is produced by a function `build_<table>(inputs…) -> pd.DataFrame`
that is **pure**: it performs no file IO, reads no globals, mutates none of its
inputs (copy before modifying), uses no randomness or wall-clock, and returns the
same output for the same input every time. Builders live one-per-table in
`build/`, with shared logic factored into pure `_underscore.py` helpers. Builders
raise `ConsistencyError` on an internal invariant breach but otherwise do not do
validation (that is the `validation/` layer's job). Builders are re-exported from
`build/__init__.py` so later stages can import and reuse them.

## 2.2 Deterministic outputs

Determinism is a hard requirement and is tested. Rules:

- Sort every output by its primary key before returning
  (`.sort_values(KEY).reset_index(drop=True)`).
- Use `sort=True` (or explicit sorting) in every `groupby`.
- Round floats to a fixed, schema-declared number of decimals (`RHO_DECIMALS`,
  `P_DECIMALS`, etc.) to eliminate binary-float noise across platforms.
- Never rely on dict/set iteration order for output ordering.
- Ranking uses a deterministic tie policy (e.g. pandas `rank(method="min")` over a
  stably-sorted frame).
- No timestamps or hostnames in outputs. The only "external" value permitted is
  the input file **SHA-256 digest** in the provenance header (a pure function of
  bytes).

## 2.3 `models / io / build / validation` layout

Non-negotiable. `models/` holds the frozen schema + errors + report carriers;
`io/` holds loaders + writers; `build/` holds pure builders; `validation/` holds
structural checks. No cross-contamination: builders never do IO; loaders never
compute derived columns; validators never mutate.

## 2.4 Orchestrator pattern

The `s{n}.py` orchestrator is thin: it composes load → build → validate → write
and assembles the `S{N}Report` (serotypes, provenance header, facts, checks). It
holds no business logic beyond wiring. Two entry points: `build_s{n}` (in-memory,
no writes) and `run_s{n}` (writes artifacts). Both return
`(S{N}Tables, S{N}Report)`.

## 2.5 CLI conventions

`argparse` in `__main__.py`, `prog="stride_s{n}"`. Input flags default to the
prior stage's conventional output directory; an explicit per-file override flag is
provided for each input; `--output-dir` defaults to `outputs_s{n}`. On success:
print `S{N} OK: <counts>, <k> checks passed.` and a second line naming the output
dir + summary JSON; return 0. On any `S{N}Error`: print
`S{N} FAILED: <ErrorType>: <msg>` to **stderr**; return 1. `main(argv=None)`
returns the int; `if __name__ == "__main__": raise SystemExit(main())`.

## 2.6 Typing conventions

- `from __future__ import annotations` at the top of every module.
- Full annotations on every function (mypy `disallow_untyped_defs`); return types
  always explicit.
- Schema constants are `Final` (`typing.Final`), tuples for column/key lists,
  `dict[str,int]` for maps.
- Prefer explicit `str(...)`/`float(...)`/`int(...)` coercions when pulling values
  out of DataFrame cells (pandas cell types are `Any`), so downstream types are
  concrete.
- When a pandas expression's static type is ambiguous to mypy (notably
  boolean-mask `DataFrame.__getitem__`), use `typing.cast(pd.DataFrame, …)` — it
  is a runtime no-op and keeps behaviour identical. (This exact issue bit S4; see
  §5.)
- No `Any` in public signatures. `# type: ignore` only with a specific error code
  and only when unavoidable (mypy `warn_unused_ignores` will flag stale ones).

## 2.7 Documentation style

NumPy-style docstrings on every public function and module. Module docstrings
open with a one-paragraph statement of purpose and cite the relevant design
sections (e.g. "§3.5, §5.2"). Each stage gets a `docs/s{n}.md` mirroring the
existing per-stage docs exactly (intro; "what it is not"; input; core concepts;
per-output schema/key/producer/consumers/invariants; validation; usage;
consumption map). Keep the "what it is not" box — it encodes the guardrails
(uncalibrated gate, no cross-serotype tests unless this IS S5, no figures, ρ read
not recomputed).

## 2.8 Testing strategy

For each new stage, add `tests/s{n}/` with a self-contained `fixtures.py` that
builds a synthetic S0-shaped (or upstream-shaped) frame directly — **do not run
upstream stages and do not depend on real data**. Cover: pure helpers vs
hand-computed numbers; each builder's columns == the schema tuple, key
uniqueness, and specific values; every validator failure branch; loader error
paths; orchestrator end-to-end with provenance/facts assertions; determinism
(double-run `assert_frame_equal`); empty-input tolerance; and CLI success/failure.
Aim for parity with existing suites (≈45–70 tests). All tests must be fully typed
(they are type-checked too).

---

# 3. Repository invariants (must always hold)

1. **Do not modify completed stages.** S0, S1A, S1B, S2, S3, S4 are frozen. A new
   stage is purely additive: new `src/stride_s{n}/`, new `tests/s{n}/`, new
   `docs/s{n}.md`, and additive edits only to `pyproject.toml`, `.gitignore`,
   `.github/workflows/ci.yml`, `README.md`, `docs/architecture.md`.
2. **No behaviour changes except genuine bug fixes.** If (and only if) a real
   correctness bug is found in a completed stage, fix the minimal code, explain
   exactly what was wrong, and re-run the full suite to prove no regression.
   Otherwise every prior-stage file stays **byte-identical**.
3. **Deterministic execution.** Same inputs ⇒ byte-identical outputs, on every
   platform and Python 3.10–3.12. Enforced by sorting, fixed rounding, and
   double-run tests.
4. **Schema stability.** Never rename/reorder/drop a column or key in a shipped
   table; never change an output filename. New stages add new tables; they do not
   alter existing ones. `*_COLUMNS`/`*_KEY`/`OUT_*` are frozen contracts.
5. **API stability.** Public `__all__` symbols, orchestrator signatures, builder
   signatures, and CLI flags of shipped stages do not change.
6. **Reproducibility / provenance.** Every stage stamps a provenance header
   (input SHA-256 digest(s), `calibrated: false`, provisional ρ\*, and any
   stage-specific params like the FDR α or ρ\* band) into its summary JSON.
7. **Validation philosophy = structural only.** Validators assert *arithmetic and
   structural* invariants (key uniqueness, counts partition, fractions in range
   and summing to 1, flags consistent with recomputed values, monotonicity where
   claimed). They make **no statistical or biological claims** and never mutate
   data. On failure they raise `ConsistencyError` and the run produces **no
   partial outputs**.
8. **ρ and variance components are read, never recomputed.** Do not reconstruct ρ
   from τ²/σ̄²/β; read the `rho` column.
9. **The gate is uncalibrated.** No stage emits a calibrated pass/fail resolution
   verdict. Resolution statements are reported over a ρ\* band (S2) or as
   ρ\*-independent descriptions (S3/S4). Provisional ρ\* = 0.5 is always labelled
   provisional.
10. **Serotype is the unit of replication (n = 4).** No stage before S5 performs a
    cross-serotype test. Within-serotype residue counts (~200) are never treated
    as independent samples. When multiple positions are screened within a
    serotype, control FDR (Benjamini–Hochberg), never raw counts.
11. **Two-tier labelling.** Domain-scale and coarser = `licensed`; residue-scale =
    `exploratory`. Every relevant row carries a `tier`.
12. **Sibling independence.** Reduction stages read S0 (and, for S2, S1B). They do
    **not** read each other's outputs. A new reduction stage should take its input
    from S0/S1 tables, not from S2/S3/S4.
13. **No figures before S7, no ML ever, no networked calls in the pipeline.**
14. **Generated outputs are never committed.** `outputs*/` and `*.parquet` are
    git-ignored (except `outputs/.gitkeep`). Real datasets under `data/` are not
    committed (except `data/README.md`). The `examples/small_synthetic_dataset/`
    IS committed.
15. **`src/stride_s{n}/build/` must be a tracked package.** The `.gitignore` has a
    **root-anchored** `/build/` rule (for the setuptools build dir) that must NOT
    shadow `src/stride_s{n}/build/`. Always verify `git check-ignore
    src/stride_s{n}/build/__init__.py` returns exit code 1 (not ignored). (See §5.)
16. **No new runtime dependencies without cause.** S4 deliberately implemented the
    normal CDF (via `math.erf`) and Benjamini–Hochberg by hand to avoid adding
    SciPy. Prefer stdlib/pandas; justify any new dependency.

---

# 4. Development workflow (exact, per stage)

Working assumptions: the repo is extracted to a working dir; a virtualenv lives at
`.venv`; activate it (`. .venv/bin/activate`) before any `python/pip/ruff/mypy/
pytest/stride-*` command.

**Step 1 — Read the design + baseline.** Read the relevant design section(s) in
`STRIDE_Downstream_Analysis_Design.docx` (convert with `pip install python-docx
--break-system-packages`, then iterate `Document(...).paragraphs`/`.tables`).
Confirm the attached baseline already contains all prior stages (list
`src/`, check `pyproject` scripts, run the suite). Record the baseline as clean:
`ruff check .`, `mypy`, `pytest` all green.

**Step 2 — Implement (additive only).** Create `src/stride_s{n}/` with the full
`models/io/build/validation` + orchestrator + CLI + `__init__`/`__main__`, mirroring
an existing stage (S3/S4 are the closest templates for an S0-reading reduction).
Write `schema.py` first (it is the contract).

**Step 3 — Integrate config (additive).** In `pyproject.toml` add: the
`stride-s{n}` script; `stride_s{n}` to coverage source; `src/stride_s{n}` to mypy
`files`; `stride_s{n}` to ruff `known-first-party`. In `.gitignore` add
`outputs_s{n}/`. In `.github/workflows/ci.yml` add `--cov=stride_s{n}` and a
smoke-run step. In `README.md` + `docs/architecture.md` add the stage description
and bump the forward-reference. Add `docs/s{n}.md`.

**Step 4 — Editable install.** `pip install -e ".[dev]"` and confirm the new
`stride-s{n}` entry point resolves (`which stride-s{n}`).

**Step 5 — Smoke the CLI.** Generate S0 output
(`stride-s0 --data-root examples/small_synthetic_dataset --output-dir /tmp/s0`),
then run the new stage against it; hand-verify a couple of numeric outputs.

**Step 6 — ruff.** `ruff check .` — fix all findings (auto-fix import order /
trailing whitespace with `ruff check . --fix`; resolve `B007` unused-loop-vars by
renaming to `_name`, etc.). Must end "All checks passed!".

**Step 7 — mypy.** `mypy` (whole repo). Must end "Success: no issues found".
Common fixes: type DataFrame-cell reads with explicit coercions; type dict-heavy
fixtures with a small frozen dataclass rather than a heterogeneous dict; use
`cast(pd.DataFrame, …)` for ambiguous boolean-mask indexing. (See §5.)

**Step 8 — pytest.** Write `tests/s{n}/`, then `pytest`. Expect
`<prior_total> + <new> passed` (prior baseline before S5 = 360). Zero
regressions.

**Step 9 — Git-tracking verification.** If not already a git repo: `git init`,
`git config user.email/name`, `git add -A`. Then verify **every** new file is
staged (`git status --short src/stride_s{n} tests/s{n} docs/s{n}.md` all show
`A`), and specifically confirm the build subpackage is not ignored:
```
git check-ignore src/stride_s{n}/build/__init__.py   # MUST exit 1 (not ignored)
```
Also confirm no new file is ignored: `git check-ignore <all new files>` returns
nothing.

**Step 10 — GitHub Actions compatibility.** Reproduce CI locally: fresh venv,
`pip install -e ".[dev]"`, `ruff check .`, `mypy`,
`pytest --cov=… --cov-report=term-missing`, then the full smoke sequence
S0→S{n} on the committed example. A stronger check: materialise a clean tree from
the git index (`git checkout-index -a -f --prefix=/tmp/fresh/`), install and run
everything there to prove nothing depends on untracked files.

**Step 11 — Byte-identical prior stages.** Extract the original baseline ZIP to a
temp dir and `diff -rq` the prior-stage `src/` packages and `tests/` suites
against the working copy (ignoring `__pycache__`). The only differences allowed
are the additive config/doc edits from Step 3 (and any explained bug fix).

**Step 12 — Package.** Build the deliverable ZIP from the git index (so it is
exactly the tracked tree, no caches/venv): `git checkout-index -a -f
--prefix=/tmp/pkg/STRIDE_downstream-main/`, verify the file count equals
`git ls-files | wc -l`, install+test the packaged tree once, then `zip -r`.

---

# 5. Known pitfalls (encountered + resolved)

**P1 — The uploaded baseline ZIP will not open with `unzip`.** The archive has a
valid End-Of-Central-Directory record followed by a trailing 40-char git-hash
comment (e.g. `090fca0765bdfa634182df14ce2ecc5a039e0769`). The `unzip` utility
reports "End-of-central-directory signature not found", but the archive is intact.
**Resolution:** extract with Python instead —
`python3 -c "import zipfile; zipfile.ZipFile('X.zip').extractall('.')"`.
`ZipFile.testzip()` returns `None` (no corruption). (ZIPs you *create* with `zip
-r` are normal and open with `unzip` fine.)

**P2 — `/build/` in `.gitignore` shadowing `src/stride_s{n}/build/`.** A naive
`build/` ignore rule would exclude every stage's `build/` package, silently
dropping pure builders from the repo. The `.gitignore` uses a **root-anchored**
`/build/` (only the top-level setuptools build dir). **Always verify** with
`git check-ignore src/stride_s{n}/build/__init__.py` → must exit **1** (not
ignored). This is an explicit deliverable requirement each stage.

**P3 — Git-tracking verification is mandatory, and the repo may not be a git repo.**
The delivered tree sometimes arrives without a `.git`. Initialise it, `git add
-A`, and confirm each new source/test/doc file shows status `A`. Count: baseline
before S5 = 193 tracked files.

**P4 — Ruff version mismatch makes CI non-reproducible.** Ruff's isort
import-ordering can differ between releases, so an unpinned Ruff (CI installs
latest) can disagree with a local install and emit `I001`. **Resolution:** Ruff is
pinned `>=0.15,<0.16` in `[project.optional-dependencies].dev`, and
`known-first-party` explicitly lists all packages + `tests` so ordering does not
depend on Ruff's src auto-detection or on whether the package is installed. Keep
both when adding a stage (append `stride_s{n}` to `known-first-party`).

**P5 — mypy on test fixtures (heterogeneous dict unpacking).** A fixtures dict
mapping names to heterogeneous tuples/dicts made values `object`, so unpacking
(`a, b, c = spec["x"]`) produced `has-type` errors and untypable variables.
**Resolution:** model the fixture spec with a small `@dataclass(frozen=True)`
holding concretely-typed fields, and type the row accumulator as
`dict[str, object]`. Do this from the start for any non-trivial fixture.

**P6 — pandas typing: boolean-mask `__getitem__` inferred as `bool`.** Under some
pandas-stubs versions,
`signed = df[df["is_signed"].astype(bool)]` is typed as `bool`, so `signed.groupby`
/ `signed.shape` raise `"bool" has no attribute …`. My local stubs accepted it;
stricter stubs did not (this shipped as a two-error mypy failure and was patched).
**Resolution:** wrap in `typing.cast(pd.DataFrame, df[mask])` — a runtime no-op
that fixes the inference regardless of stub version. Do not "fix" it by changing
the filter logic (behaviour must not change). More generally: prefer explicit
scalar coercions (`str()/float()/int()`) when reading cells, and `cast` for
ambiguous frame expressions.

**P7 — Editable install / import visibility.** The project is src-layout;
`stride_*` packages are importable only after `pip install -e ".[dev]"` **and**
with the venv active. Forgetting to activate the venv yields
`ModuleNotFoundError`. (pytest still finds them via `pythonpath=["src"]`, but the
CLIs and ad-hoc `python -c` imports need the install + active venv.)

**P8 — Domain-scale values are region-constant; assert it, don't assume silently.**
When a builder reads a domain-scale row per `(serotype, chain, domain)`, multiple
member loci carry identical domain rows. Collapse them to one row and **assert**
the region-constant invariant (raise `ConsistencyError` if τ²/σ̄²/ρ/β differ
within a region), as S2's `domain_reproducibility` and S4's `domain_regions` do.

**P9 — Avoid new heavy dependencies.** S4 needed a normal-CDF p-value and BH-FDR;
adding SciPy would bloat the dependency set. **Resolution:** implement
`normal_cdf` via `math.erf` and Benjamini–Hochberg by hand (both in
`stride_s4/build/_stats.py`), fully tested against known values. Follow this
preference.

**P10 — Deterministic ranking and grouping.** Unsorted `groupby` or default
ranking tie-handling can reorder outputs across runs/platforms. Always
`sort=True`, sort by key before returning, and choose an explicit `rank(method=…)`.
Tested via double-run `assert_frame_equal`.

**P11 — Summary JSON must be stable.** Only include deterministic content (counts,
sorted serotypes, digests). No timestamps/hostnames. Keep the payload key order
identical to the existing writers.

---

# 6. Current progress

## 6.1 Complete (shipped, all green)

- **S0 — Ingest & validate (`stride_analysis`, v0.2.0).** Discovers the 8 raw
  input files, loads and schema-validates them (profile + mechanism + hierarchy +
  correlations), checks profile↔mechanism consistency and gate-uncertain flags,
  and builds the canonical **`stride_table.parquet`** (45-col tidy master, one row
  per serotype×locus×scale) and **`replicate_table.parquet`**, plus
  `schema_report.json` and a human-readable `validation_report.md`. This is the
  single source of truth every downstream stage reads.
- **S1A — Canonical join layer (`stride_s1a`, v0.1.0).** From S0 tables, builds
  `canonical_residues`, `domain_table`, `replicate_inventory`, and the
  `conservation_table` (the 199-shared-position index: `n_serotypes`,
  `serotypes_present` per `canon_label`), + `dataset_summary.json`.
- **S1B — Biological annotation layer (`stride_s1b`, v0.1.0).** From S1A, builds
  `residue_annotation` (adds `conservation_class ∈ {pan_serotype, partial,
  serotype_unique}`, domain/SS status, availability), `domain_annotation`,
  `hierarchy_annotation`, `serotype_annotation`, + `annotation_summary.json`.
- **S2 — Per-serotype reduction (`stride_s2`, v0.1.0).** From S0 tidy + S1B
  annotation, builds `resolution_census`, `residue_landscape`,
  `domain_reproducibility`, `signed_screen`, `serotype_summary`, +
  `reduction_summary.json`. Every resolution quantity is reported **over a ρ\*
  band** (default 0.50→0.90 step 0.05), never a single threshold; rows are
  tier-labelled; no cross-serotype tests.
- **S3 — Hierarchy reduction (`stride_s3`, v0.1.0).** From S0 tidy only, reduces
  *along the scale axis*: `scale_curve` (ρ-vs-scale per locus), `resolution_gap`
  (Δρ domain−residue, distributed-effect flag), `monotonicity_audit` (I2
  upward-closure), `chain_contrast` (NS2B vs NS3), + `hierarchy_summary.json`.
  ρ\*-independent; sibling to S2.
- **S4 — Uncertainty layer (`stride_s4`, v0.1.0).** From S0 tidy only:
  `variance_budget` (τ²/σ̄² per domain + replicate/sampling regime, Tier A),
  `residue_variance` (per-residue decomposition + τ² rank, Tier B),
  `significance_screen` (CI-exclusion + two-sided Wald p via `math.erf` +
  Benjamini–Hochberg FDR within each serotype), `domain_effect_summary`
  (1/β_se²-weighted mean effect + per-domain CI-exclusion fraction, Tier A), +
  `uncertainty_summary.json`. ρ\*-independent; sibling to S2/S3; no SciPy.

## 6.2 Remaining (not implemented)

- **S5 — Cross-serotype (n = 4).** Conservation of reproducibility & direction
  across serotypes over the 199 shared positions; domain×serotype matrices.
  **This is the next stage.** (Details in §7.)
- **S6 — Replicate layer.** Rank concordance / per-run θ correlation /
  leave-one-replicate-out stability. **Blocked:** requires per-run correlation
  CSVs that are not in the current dataset (only the K=3 aggregate τ²/σ̄² is
  available). Do not attempt until that data exists.
- **S7 — Figures & tables.** Publication outputs F1–F8 / T1–T5, deterministic and
  provenance-stamped. The consumer of all prior tables.

---

# 7. Next stage — S5 (Cross-serotype, n = 4)  — SPECIFICATION ONLY, DO NOT IMPLEMENT

**Design row (verbatim):** *S5 — Cross-serotype (n=4) — shared-position frame →
conservation — "Conservation of reproducibility & direction across serotypes,
with serotype as the replication unit. Domain×serotype matrices."*

**Goal.** Produce the cross-serotype conservation layer, treating **serotype as
the unit of replication (n = 4)** and never treating residues/frames as
independent. This is the first and only stage permitted to compare across
serotypes.

**Inputs (read-only; sibling to S2/S3/S4).** The S0 `stride_table.parquet`
(for ρ, direction, gated scale, variance components per serotype×position) and the
S1A `conservation_table.parquet` and/or S1B `residue_annotation.parquet` (for the
199-shared-position index and `conservation_class`). S5 should **aggregate to one
value per serotype per region first**, then compare across the four serotypes. It
must not consume S2/S3/S4 outputs.

**Products the design calls for (§3.3, §5.2):**
- **Conservation of reproducibility** across the 199 shared `canon_label`
  positions: for each shared position, is it reproducible (gated at the
  provisional ρ\*) in all four / some / none of the serotypes — with an explicit
  `n_serotypes_reproducible` count and a conservation label.
- **Serotype-divergent positions:** reproducible-and-signed in one serotype but
  absent/mixed in another (candidate serotype-specific mechanism differences).
- **Direction concordance:** for shared, signed positions, do serotypes agree on
  increase vs decrease? (e.g. all-agree / majority / conflicting.)
- **Domain × serotype matrix:** ρ(domain × serotype) over the NS3 domains + NS2B
  — a tidy long table keyed by (domain, serotype), suitable for a future heatmap.
- **Conserved catalytic-machinery check:** behaviour of the Catalytic Triad
  (NS3:51/75/135) and Oxyanion Loop across serotypes (expected most conserved).
- **Cross-serotype scorecard** (T4/T1 flavour): per-serotype
  (n_loci, %residue-gated, %signed, %mixed, ρ median) aggregated to one row per
  serotype.

**Statistical guardrails (mandatory, §5.2):** the replication unit is the serotype
(n = 4). Prefer **descriptive statistics and effect sizes over p-values at n = 4**;
if any test is reported, it is across the four serotype-level aggregates, never
across residues. No calibrated verdicts (gate still uncalibrated; use provisional
ρ\* = 0.5 and stamp `calibrated: false`). Two-tier labelling still applies
(domain-scale conservation is licensed; residue-scale is exploratory). No figures
(that is S7).

**Suggested outputs** (to be finalised against the design during implementation;
4 tables + summary JSON, mirroring the established pattern):
`position_conservation.parquet` (per shared canon_label), `direction_concordance.
parquet` (per shared signed canon_label), `domain_serotype_matrix.parquet` (per
domain×serotype), `cross_serotype_scorecard.parquet` (per serotype), +
`conservation_summary.json`. Package `stride_s5`, CLI `stride-s5` (inputs default
to `outputs` for the STRIDE table and `outputs_s1a`/`outputs_s1b` for the shared
index; `--output-dir` default `outputs_s5`).

**Do NOT implement now.** This section is the brief for the next conversation.

---

# 8. Verification checklist (complete before committing ANY stage)

Run every item; all must pass. (Activate the venv first.)

**Build / install**
- [ ] `pip install -e ".[dev]"` succeeds.
- [ ] `which stride-s{n}` resolves (new entry point installed).

**Static analysis**
- [ ] `ruff check .` → "All checks passed!" (no `I001`/`B007`/`W291`/…).
- [ ] `mypy` → "Success: no issues found in N source files" (N grows by the new
      package's file count; before S5 it is 155).

**Tests**
- [ ] `pytest` → `<prior_total + new> passed`, **0 failed** (prior baseline = 360).
- [ ] New `tests/s{n}/` covers: pure helpers vs hand-computed values; each
      builder's columns == schema tuple + key uniqueness + specific values; every
      validator failure path raises `ConsistencyError`; loader error paths;
      orchestrator end-to-end (provenance `calibrated:false`, digests, facts);
      **determinism** (double-run `assert_frame_equal`); empty-input tolerance;
      CLI success (rc 0, "S{N} OK") and failure (rc 1, "S{N} FAILED").

**CLI / functional**
- [ ] Full smoke sequence runs: S0 → … → S{n} on
      `examples/small_synthetic_dataset`, producing all expected artifacts.
- [ ] Hand-verify at least two numeric outputs against a manual calculation.

**Git tracking**
- [ ] `git add -A`; every new `src/stride_s{n}/*.py`, `tests/s{n}/*.py`, and
      `docs/s{n}.md` shows status `A`.
- [ ] `git check-ignore src/stride_s{n}/build/__init__.py` → **exit code 1** (NOT
      ignored).
- [ ] `git check-ignore <all new files>` → prints nothing (none ignored).
- [ ] Tracked-file count increased only by the new files (record the new total).

**GitHub Actions compatibility**
- [ ] `.github/workflows/ci.yml` updated: `--cov=stride_s{n}` added to the pytest
      line; a "Smoke-run S{n}" step added after the previous stage.
- [ ] Reproduce CI locally in a **fresh** venv (or from a clean
      `git checkout-index` tree): install → ruff → mypy →
      `pytest --cov=… --cov-report=term-missing` → full S0→S{n} smoke — all green.

**Additive-only integrity**
- [ ] `pyproject.toml` changes are additive: new script; `stride_s{n}` appended to
      coverage source, mypy `files`, ruff `known-first-party`.
- [ ] `.gitignore` adds only `outputs_s{n}/`; the root-anchored `/build/` rule
      untouched.
- [ ] `README.md` + `docs/architecture.md` updated (status/description/quick-start/
      tree/docs-list + forward-reference bump); `docs/s{n}.md` added.
- [ ] `diff -rq` of every **prior** stage's `src/` package and `tests/` suite
      against the original baseline shows **no differences** (byte-identical),
      except any explicitly explained + re-verified bug fix.

**Packaging (for delivery)**
- [ ] Build the ZIP from the git index (`git checkout-index -a -f --prefix=…`),
      confirm file count == `git ls-files | wc -l`, install + `pytest` the
      packaged tree once, strip caches/venv, `zip -r`.

---

*End of master handoff. Attach this document plus the repository ZIP to a new
conversation to continue with S5.*
