# Stage S1B — the biological annotation layer

S1B is the biological interpretation layer. It consumes **only** the four S1A
parquet tables and derives reusable, deterministic biological annotation tables
that later stages consume.

> **What S1B is not.** S1B performs no statistical testing, no ranking, no
> clustering, no machine learning, no hypothesis generation, and produces no
> figures. Every derived column is a categorical or boolean label produced by a
> fixed rule from S1A facts. "Interpretation" here means deterministic
> structural classification, not inference.

S1B never reads raw STRIDE outputs, `profile.csv`/`mechanism.json`, or MD
trajectories. Its only inputs are the S1A tables.

---

## Inputs

| Input | Produced by | Columns S1B reads |
|---|---|---|
| `canonical_residues.parquet` | S1A | `serotype, canon_label, chain, domain, motif, secondary_structure, hierarchy_path, complex, protein, residue` |
| `domain_table.parquet` | S1A | `serotype, complex, protein, chain, domain, n_residues, hierarchy_path` |
| `replicate_inventory.parquet` | S1A | `serotype, canon_label, n_replicates, available, in_all_replicates` |
| `conservation_table.parquet` | S1A | `canon_label, n_serotypes, …` |

---

## Deterministic annotation vocabularies

Every derived label comes from one of these closed vocabularies, defined in
`stride_s1b.models.schema`:

| Vocabulary | Values | Rule |
|---|---|---|
| domain status | `assigned`, `unassigned` | `assigned` unless the domain value is a STRIDE sentinel (`unassigned`/`none`/`unknown`/empty) |
| secondary-structure status | `resolved`, `unresolved` | `resolved` unless the SS value is a sentinel |
| conservation class | `pan_serotype`, `partial`, `serotype_unique` | present in all / some-but-not-all / exactly one serotype |
| availability class | `all_replicates`, `some_replicates`, `no_replicates` | observed in every / some / no replicates of its serotype |

The sentinel comparison is case-insensitive and whitespace-tolerant; the
underlying S1A values are never rewritten.

---

## Outputs

Five artifacts are written to the S1B output directory.

### 1. `residue_annotation.parquet`

Per-residue biological annotation.

- **Schema:** `serotype, canon_label, chain, domain, hierarchy_path, domain_status, secondary_structure_status, conservation_class, n_serotypes_present, availability_class, n_replicates`
- **Primary key:** `(serotype, canon_label)`
- **Producer:** S1B (`build_residue_annotation`), from canonical residues +
  conservation table + replicate inventory.
- **Downstream consumers:** any stage that filters or groups residues by
  structural category — e.g. "annotated residues in an assigned domain that are
  pan-serotype and seen in all replicates". This is the central S1B table.
- **Invariants:** exactly one row per canonical residue (bijection with the S1A
  canonical residues — no missing, duplicate, or orphan annotations); every
  category value is drawn from its closed vocabulary.

### 2. `domain_annotation.parquet`

Per-domain structural annotation.

- **Schema:** `serotype, complex, protein, chain, domain, hierarchy_path, domain_status, n_residues, n_pan_serotype_residues, fully_conserved`
- **Primary key:** `(serotype, chain, domain)`
- **Producer:** S1B (`build_domain_annotation`), from the S1A domain table +
  the residue annotation.
- **Downstream consumers:** stages that operate at domain granularity (domain
  rollups, domain-level reporting) and want a domain's conservation composition
  without recomputing it.
- **Invariants:** unique key; `n_residues` equals the number of residue
  annotations mapping to the domain; every residue's `(serotype, chain, domain)`
  exists here; `fully_conserved` is true iff every member residue is
  `pan_serotype`.

### 3. `hierarchy_annotation.parquet`

Per-residue hierarchy resolution annotation.

- **Schema:** `serotype, hierarchy_path, canon_label, complex, protein, chain, domain, motif, secondary_structure, n_levels_total, n_levels_resolved, fully_resolved`
- **Primary key:** `(serotype, hierarchy_path)`
- **Producer:** S1B (`build_hierarchy_annotation`), from the canonical residues.
- **Downstream consumers:** stages that need structural-path completeness (how
  fully a residue's position is resolved through the 7-level grammar) without
  re-parsing paths.
- **Invariants:** unique key (one residue per `(serotype, hierarchy_path)`);
  `n_levels_total == 7`; `n_levels_resolved` counts non-sentinel levels;
  `fully_resolved` iff all seven levels are resolved; 1:1 with the residue
  annotation on `(serotype, hierarchy_path)`.

### 4. `serotype_annotation.parquet`

Per-serotype structural composition.

- **Schema:** `serotype, n_residues, n_domains, n_assigned_domain_residues, n_unassigned_domain_residues, n_pan_serotype_residues, n_partial_residues, n_serotype_unique_residues, n_residues_all_replicates, n_residues_some_replicates, n_residues_no_replicates`
- **Primary key:** `(serotype,)`
- **Producer:** S1B (`build_serotype_annotation`), by tallying the residue and
  domain annotations.
- **Downstream consumers:** dataset-level summaries and any stage needing a
  serotype's category composition at a glance.
- **Invariants:** unique key; every serotype exists in the S1A canonical
  residues; `n_residues` equals the serotype's residue-annotation row count; the
  three conservation counts sum to `n_residues`; likewise the three availability
  counts.

### 5. `annotation_summary.json`

Machine-readable run summary: serotypes, table sizes, category-count facts
(conservation, availability, domain-status distributions; number of
fully-resolved hierarchy paths), and every structural validation outcome.
Consumed by CI smoke checks and humans auditing a run.

---

## Validation (structural only)

Before the tables are trusted, S1B verifies — with no statistical assertions:

- **every residue has exactly one biological annotation** (bijection with the
  S1A canonical residues; no missing, duplicate, or orphan annotations);
- **hierarchy paths remain unique** (one row per `(serotype, hierarchy_path)`);
- **every domain membership is internally consistent** (unique domain key; every
  residue's domain exists; declared `n_residues` matches the mapped count);
- **every serotype references existing canonical residues** (unique key; known
  serotypes; declared counts match);
- **referential integrity between every generated table** (residue ↔ hierarchy
  1:1 on path; residue domains ⊆ domain annotation; the serotype annotation
  covers exactly the serotypes present in the other tables).

Any failure raises a `ConsistencyError` (a subclass of `S1BError`); the run does
not produce partial outputs on failure.

---

## Usage

```bash
# after S1A has written its tables to outputs_s1a/
stride-s1b --input-dir outputs_s1a --output-dir outputs_s1b
# or point at the four tables explicitly:
stride-s1b --canonical-residues outputs_s1a/canonical_residues.parquet \
           --domain-table outputs_s1a/domain_table.parquet \
           --replicate-inventory outputs_s1a/replicate_inventory.parquet \
           --conservation-table outputs_s1a/conservation_table.parquet \
           --output-dir outputs_s1b
```

Python API:

```python
from stride_s1b import run_s1b, build_s1b

# build + write artifacts
tables, report = run_s1b(
    "outputs_s1a/canonical_residues.parquet",
    "outputs_s1a/domain_table.parquet",
    "outputs_s1a/replicate_inventory.parquet",
    "outputs_s1a/conservation_table.parquet",
    "outputs_s1b",
)

# build in memory only (no writes)
tables, report = build_s1b(cr_path, dt_path, ri_path, ct_path)
assert report.all_passed
tables.residue_annotation     # the four tables are attributes of S1BTables
tables.domain_annotation
tables.hierarchy_annotation
tables.serotype_annotation
```

Reusing individual builders in a later stage:

```python
from stride_s1b.build import build_residue_annotation, build_domain_annotation
```

---

## Consumption map (quick reference)

| Table | Primary later consumers |
|---|---|
| `residue_annotation` | residue-level filtering/grouping by structural category |
| `domain_annotation` | domain-level rollups and conservation composition |
| `hierarchy_annotation` | structural-path resolution completeness |
| `serotype_annotation` | dataset-level category composition per serotype |

S1B hands later stages a trustworthy, purely-structural annotation layer. No
statistics, ranking, or inference happen here.
