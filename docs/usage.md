# Usage

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Command line

```bash
stride-s0 --data-root <path> --output-dir <path>
# equivalently:
python -m stride_analysis --data-root <path> --output-dir <path>
```

Against the committed example:

```bash
stride-s0 --data-root examples/small_synthetic_dataset --output-dir outputs
```

### Flags

| Flag | Effect |
|---|---|
| `--data-root PATH` | dataset root (default `data`). Any path works; nothing is hardcoded. |
| `--output-dir PATH` | where artifacts are written (default `outputs`). |
| `--no-require-replicates` | allow serotypes with no Level-1 replicate data (summaries-only). |
| `--no-require-summaries` | allow serotypes with no Level-2 summaries (replicates-only). |
| `--allow-unequal-replicates` | do not require equal replicate counts across serotypes. |
| `--strict-cross-level` | require replicate residue labels to be a subset of profile labels. |

The command exits `0` on success and `1` on any `StrideAnalysisError`, printing
an actionable message to stderr — suitable for CI.

## Python API

```python
from stride_analysis import run_s0, build_tables, discover_dataset

# discover only (inspect structure before building)
dataset = discover_dataset("data")
for s in dataset.serotypes:
    print(s.serotype, s.n_replicates, s.summary is not None)

# build both tables in memory (no writes)
replicate_table, stride_table, report = build_tables("data")
assert report.all_passed

# build + write artifacts
replicate_table, stride_table, report = run_s0("data", "outputs")
```

`build_tables` / `run_s0` accept the same toggles as the CLI:

```python
build_tables(
    "data",
    require_replicates=False,        # summaries-only OK
    require_summaries=True,
    enforce_equal_replicate_counts=True,
    strict_cross_level=False,
)
```

## Reusing subpackages in a later stage

The framework is not a monolith; import what you need:

```python
from stride_analysis.io import discover_dataset, load_correlations
from stride_analysis.validation import validate_correlations_schema
from stride_analysis.canonical import assemble_replicate_table
```

## Partial datasets

Real projects arrive incrementally. Both levels are independently optional:

```bash
# only summaries so far
stride-s0 --data-root data --output-dir outputs --no-require-replicates

# only replicate correlations so far
stride-s0 --data-root data --output-dir outputs --no-require-summaries
```

When a level is absent, its canonical table is written only if non-empty; the
validation report records what was ingested.

## Outputs

Written to `--output-dir`:

| File | Contents |
|---|---|
| `replicate_table.parquet` / `.csv` | Level-1 canonical table |
| `stride_table.parquet` / `.csv` | Level-2 canonical table |
| `schema_report.json` | machine-readable schema facts + check outcomes |
| `validation_report.md` | human-readable validation summary |

All outputs are regenerable and git-ignored.
