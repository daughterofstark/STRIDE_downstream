# Architecture

The framework is organised as **reusable subpackages** with a thin orchestration
layer on top. There is deliberately **no monolithic pipeline**: the S0
orchestrator (`s0.py`) only *composes* the subpackages, and later stages import
the same subpackages directly.

```
                 ┌──────────────┐
                 │  __main__.py │  CLI
                 └──────┬───────┘
                        │
                 ┌──────▼───────┐
                 │    s0.py     │  orchestration (compose, collect report)
                 └──────┬───────┘
        ┌───────────────┼───────────────┬──────────────┐
        ▼               ▼               ▼              ▼
   ┌─────────┐    ┌───────────┐   ┌───────────┐  ┌──────────┐
   │   io/   │    │validation/│   │ canonical/│  │reporting │
   └────┬────┘    └─────┬─────┘   └─────┬─────┘  └──────────┘
        │               │               │
        └───────────────┴──────┬────────┘
                               ▼
                          ┌─────────┐
                          │ models/ │  schema, pydantic models, errors
                          └─────────┘
```

## Subpackages

### `models/`
The foundation. Contains:
- `schema.py` — every frozen schema constant (Level-1 required/optional columns,
  Level-2 profile columns, mechanism fields, the hierarchy grammar, and the two
  canonical-table column lists). Single source of truth.
- `__init__.py` — pydantic models for the mechanism JSON, dataclass descriptors
  for discovered datasets (`Dataset`, `SerotypeDataset`, `ReplicateInput`,
  `SummaryInput`), and the `Report` carrier.
- `errors.py` — the typed exception hierarchy rooted at `StrideAnalysisError`.

### `io/`
- `discovery.py` — resolves a data root into a `Dataset`, auto-detecting the
  nested run-dir layout. Enforces structural rules (matched pairs, equal
  replicate counts, no duplicate serotypes) without hardcoding paths.
- `loaders.py` — thin readers (CSV → DataFrame, JSON → validated `MechanismFile`).

### `validation/`
- `hierarchy.py` — `region_id` path parsing against the frozen grammar.
- `replicate.py` — Level-1 schema validation (required core + tolerated optional).
- `summary.py` — Level-2 profile schema + profile↔mechanism consistency.
- `cross_level.py` — advisory replicate↔summary residue-alignment check.

### `canonical/`
- `replicate_table.py` — builds/assembles the replicate table.
- `stride_table.py` — builds/assembles the STRIDE table.
The two builders are independent; the tables are never merged.

### `reporting.py`
Serialises the `Report` into `schema_report.json` and `validation_report.md`.

### `s0.py`
`build_tables()` (in memory) and `run_s0()` (writes artifacts). Pure composition.

## Dependency direction

`models/` depends on nothing else in the package. `io/`, `validation/`, and
`canonical/` depend only on `models/`. `s0.py` depends on all of them.
Nothing depends on `s0.py` except the CLI. This keeps every layer importable and
testable in isolation, and lets S1+ reuse `io`/`validation`/`canonical` without
importing any S0 orchestration.

## Error model

Every failure is a typed exception:

| Exception | Raised when |
|---|---|
| `DiscoveryError` | dataset layout problems (missing/duplicate/unmatched/inconsistent) |
| `SchemaError` | a file violates its frozen schema |
| `HierarchyError` | a `region_id` cannot be parsed |
| `ConsistencyError` | a cross-file/cross-row invariant fails |

All subclass `StrideAnalysisError`, so callers can catch everything with one
`except`.

## Downstream stages

**S1A** (`src/stride_s1a/`) is implemented and follows the same architecture as
S0 — reusable `models/`, `io/`, `build/`, and `validation/` subpackages with a
thin `s1a.py` orchestrator and a CLI. It consumes **only** the S0 canonical
parquet tables (never the raw STRIDE files) and builds the biological data layer
(canonical residues, domain summaries, replicate availability, conservation
mapping). See [`s1a.md`](s1a.md).

**S1B** (`src/stride_s1b/`) is implemented and follows the same architecture. It
consumes **only** the four S1A parquet tables (never S0 outputs, raw STRIDE
files, or trajectories) and builds the biological *annotation* layer: per-residue,
per-domain, per-hierarchy, and per-serotype deterministic structural annotations.
Every derived column is a categorical or boolean label from a closed vocabulary —
no statistics, ranking, clustering, or figures. See [`s1b.md`](s1b.md).

**S2** (`src/stride_s2/`) is implemented and follows the same architecture. It
consumes **only** the S0 STRIDE table and the two S1B annotation tables, and
builds the per-serotype *reduction* layer: the achieved-resolution census, the
residue-scale reproducibility landscape, domain-scale reproducibility, the
signed/significant screen, and a per-serotype scorecard. Every resolution
quantity is reported over a **ρ\* band** (the gate is uncalibrated), rows are
labelled `licensed` (domain-scale) or `exploratory` (residue-scale), and there
are **no cross-serotype tests** (those are S5). See [`s2.md`](s2.md).

**S3** (`src/stride_s3/`) is implemented and follows the same architecture. It
consumes **only** the S0 STRIDE table (the profile) and reduces it *along the
scale axis* into the hierarchy layer: the ρ-vs-scale curve per locus, the
domain−residue reproducibility gap Δρ, the monotonicity (I2 upward-closure)
audit, and the chain-level contrast (e.g. NS2B vs NS3). It is a **sibling
reduction to S2** — both read the S0 profile independently; S3 does not consume
S2's outputs. Its products are ρ\*-independent descriptions of the profile, rows
are labelled `licensed` (domain-scale and coarser) or `exploratory`
(residue-scale), and there are **no cross-serotype tests** (those are S5). See
[`s3.md`](s3.md).

Further stages (S4+) import the canonical tables and/or the S1A/S1B/S2/S3 tables
(or their builders) and add their own module. They must not modify S0, S1A, S1B,
S2, or S3, must keep the data levels separate, and must not commit generated
outputs. See `CONTRIBUTING.md`.
