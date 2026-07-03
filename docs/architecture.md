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

## Extending to S1 (not implemented here)

A later stage imports the canonical tables (or the builders) and adds its own
module. It must not modify S0, must keep the two data levels separate, and must
not commit generated outputs. See `CONTRIBUTING.md`.
