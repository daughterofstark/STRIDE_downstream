# Stage S1A — the reusable biological data layer

S1A is the first biological data layer. It consumes **only** the S0 canonical
tables (`stride_table.parquet` and `replicate_table.parquet`) and builds reusable
structural annotations and derived datasets that every later biological stage
will import.

> **What S1A is not.** S1A performs no statistical testing, no hypothesis
> generation, no clustering, no machine learning, no ranking, no mechanism
> inference, and produces no figures. It builds the data layer only. Those
> analyses belong to S1B and later.

S1A never re-reads the raw STRIDE CSVs/JSONs — once S0 has built the canonical
tables, they are the sole upstream input.

---

## Inputs

| Input | Produced by | Columns S1A reads |
|---|---|---|
| `stride_table.parquet` | S0 | `serotype, canon_label, scale_level, locus, region_id, h_complex … h_residue` (residue-scale rows) |
| `replicate_table.parquet` | S0 | `serotype, replicate, replicate_index, canon_label` |

A **missing** `replicate_table.parquet` is treated as a valid "summaries-only"
state (S0 only writes it when non-empty): the replicate inventory then records
every residue as unavailable rather than failing.

---

## Outputs

Five artifacts are written to the S1A output directory.

### 1. `canonical_residues.parquet` — Task 1

The canonical residue object: one row per **(serotype, canon_label)**.

| Column | Meaning |
|---|---|
| `serotype` | serotype the residue was measured in |
| `canon_label` | canonical residue label (cross-serotype identity key) |
| `chain`, `domain`, `motif`, `secondary_structure` | structural annotations (from the S0 hierarchy) |
| `hierarchy_path` | full `complex/…/residue` path (the residue-scale `region_id`) |
| `complex`, `protein`, `residue` | the remaining hierarchy levels, materialised |

- **Unique key:** `(serotype, canon_label)`.
- Hierarchy values are preserved verbatim, including STRIDE's `unassigned` /
  `none` / `unknown` sentinels — S1A does not clean or reinterpret them.
- **Consumed by:** every later stage that needs residue identity or structural
  context; it is the spine the other three tables and all of S1B join against.

### 2. `domain_table.parquet` — Task 3

Structural domain summaries: one row per **(serotype, chain, domain)**.

| Column | Meaning |
|---|---|
| `serotype`, `complex`, `protein`, `chain`, `domain` | domain identity + hierarchy membership |
| `n_residues` | number of canonical residues in the domain (a count, not a score) |
| `canon_labels` | sorted list of member residue identifiers |
| `hierarchy_path` | `complex/protein/chain/domain` |

- **Unique key:** `(serotype, chain, domain)`.
- Structural counts only — no biological scoring, no ranking.
- **Consumed by:** later stages that aggregate or summarise at domain
  granularity (e.g. domain-level rollups in S1B), and any reporting that needs a
  domain roster per serotype.

### 3. `replicate_inventory.parquet` — Task 4

Per-residue replicate availability: one row per **(serotype, canon_label)**.

| Column | Meaning |
|---|---|
| `n_replicates` | how many replicates observed this residue |
| `replicates` | sorted list of replicate names (e.g. `1st_run`) |
| `replicate_indices` | sorted list of 1-based replicate indices |
| `available` | `n_replicates > 0` |
| `in_all_replicates` | observed in *every* replicate of the serotype |

- **Unique key:** `(serotype, canon_label)`.
- **Availability only.** S1A never averages replicates, never aggregates effect
  values, never computes significance — it records presence and nothing else.
- Residues with no replicate observations are recorded with `n_replicates == 0`
  and `available == False`, never dropped.
- **Consumed by:** any later stage that must gate an analysis on replicate
  availability (e.g. "only analyse residues seen in all replicates"). S1B reads
  this to decide which residues are eligible before it does any statistics.

### 4. `conservation_table.parquet` — Task 2

Cross-serotype presence map: one row per **`canon_label`** (the union).

| Column | Meaning |
|---|---|
| `n_serotypes` | in how many serotypes the residue is present |
| `serotypes_present` | sorted list of serotypes containing it |
| `serotypes_absent` | sorted list of serotypes lacking it |
| `in_all_serotypes` | present in every serotype (the intersection) |
| `in_any_serotype` | always `True` by construction (union membership) |
| `chain`, `domain` | structural annotation (consistent across serotypes) |

- **Unique key:** `(canon_label,)`.
- Pure set membership — no sequence alignment, no conservation *scoring*, no
  interpretation of why a residue is present or absent.
- **Consumed by:** later cross-serotype stages that compare conserved vs.
  serotype-specific residues. S1B and comparison stages read this to restrict
  attention to the intersection (or to study the set difference) before any
  testing.

### 5. `dataset_summary.json`

Machine-readable run summary: serotypes seen, table sizes, conservation counts,
per-serotype residue/domain facts, and the outcome of every Task-5 validation
check. Consumed by CI smoke checks and by humans auditing a run.

---

## Validation (Task 5)

Before the tables are trusted, S1A verifies:

- **every STRIDE locus maps to exactly one canonical residue** (and vice-versa —
  no orphan residues);
- **every canonical residue has exactly one hierarchy path**;
- **every replicate row maps to a canonical residue** (no orphan replicate rows);
- **a shared `canon_label` carries consistent structural annotation** across
  serotypes (so the conservation map never conflates distinct positions).

Any failure raises a `ConsistencyError` (a subclass of `S1AError`) with an
actionable message; the run does not produce partial outputs on failure.

---

## Usage

```bash
# after S0 has written its tables to outputs/
stride-s1a --input-dir outputs --output-dir outputs_s1a
# or point at the tables explicitly:
stride-s1a --stride-table outputs/stride_table.parquet \
           --replicate-table outputs/replicate_table.parquet \
           --output-dir outputs_s1a
```

Python API:

```python
from stride_s1a import run_s1a, build_s1a

# build + write artifacts
tables, report = run_s1a("outputs/stride_table.parquet",
                         "outputs/replicate_table.parquet",
                         "outputs_s1a")

# build in memory only (no writes)
tables, report = build_s1a("outputs/stride_table.parquet",
                           "outputs/replicate_table.parquet")
assert report.all_passed
tables.canonical_residues   # the four tables are attributes of S1ATables
tables.domain_table
tables.replicate_inventory
tables.conservation_table
```

Reusing individual builders in a later stage:

```python
from stride_s1a.build import build_canonical_residues, build_conservation_table
```

---

## Consumption map (quick reference)

| Table | Primary later consumers |
|---|---|
| `canonical_residues` | all stages — residue identity & structural context |
| `domain_table` | domain-level rollups / reporting |
| `replicate_inventory` | availability gating before any statistics |
| `conservation_table` | cross-serotype comparison (conserved vs specific) |

S1A hands later stages a trustworthy, purely-structural biological layer. No
statistics, ranking, or interpretation happen until S1B.
