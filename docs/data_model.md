# Data model

S0 produces **two separate canonical tables**. They are never merged — replicate
observations (Level 1) and STRIDE summaries (Level 2) describe the same residues
at different levels of derivation, and later stages consume whichever they need.

---

## Replicate table (Level 1)

The primary per-run observations, ingested from every
`*_correlations_v5.csv`.

- **Grain:** one row per **(serotype, replicate, residue)**.
- **Unique key:** `(serotype, replicate, canon_label)`.
- **Output:** `outputs/replicate_table.{parquet,csv}`.

### Columns

Leading identity + provenance columns are prepended; the original STRIDE
correlation columns are preserved verbatim after them.

| Column | Meaning |
|---|---|
| `serotype` | serotype tag (from the directory name) |
| `replicate` | run-directory name (e.g. `1st_run`) |
| `replicate_index` | 1-based replicate order (sorted run dirs) |
| `canon_label` | canonical residue label — the join key to the STRIDE table |
| `file_resid`, `canon_resid`, `name` | residue identity as reported by STRIDE |
| `label` | STRIDE's own residue label (== `canon_label`) |
| `r`, `abs_r` | the per-replicate effect field and its magnitude |
| *(known-optional)* | `p_raw`, `p_bonf`, `rmsf`, `theta_se`, `n_eff`, `chain`, `domain`, … when present |
| `source_path`, `run_dir` | provenance |

Because STRIDE appends columns across its milestones, the replicate schema is
validated as a **required core** (`file_resid, canon_resid, name, label, r,
abs_r`) plus **known-optional** columns (type-checked when present). Genuinely
unknown extra columns are preserved and reported, never dropped. When replicates
carry different column sets, the assembled table is the column union with NA in
the gaps.

---

## STRIDE table (Level 2)

The triplicate summaries, ingested from each `*_profile.csv` +
`*_mechanism.json` pair.

- **Grain:** one row per **(serotype, canon_label, scale_level)**.
- **Unique key:** `(serotype, canon_label, scale_level)`.
- **Output:** `outputs/stride_table.{parquet,csv}`.

### Column groups

| Group | Columns |
|---|---|
| identity | `serotype, canon_label, scale_level, scale_index, locus, region_id` |
| profile payload | `rho, gated, beta, beta_se, tau2, sigma2_bar, a_signed, coherence, method, status` |
| explicit hierarchy | `h_complex, h_protein, h_chain, h_domain, h_motif, h_secondary_structure, h_residue` |
| gate flag | `is_gated_scale` |
| mechanism payload | `mech_label, mech_direction, mech_beta_signed, mech_beta_ci_lower, mech_beta_ci_upper, mech_beta_se, mech_coherence, mech_reproducible_magnitude_energy, mech_rho_star, mech_calibrated, mech_gate_uncertain, mech_status, mech_region_id, mech_n_loci` |
| provenance | `profile_source, mechanism_source, gate_rho_star, gate_alpha, gate_coherence_threshold, mechanism_calibrated, mechanism_schema_version` |

### The hierarchy

Each locus's residue-scale `region_id` is a 7-segment path that is parsed once
into the explicit `h_*` columns:

```
complex / protein / chain / domain / motif / secondary_structure / residue
```

A region at a coarser scale is a prefix of this path. Every profile row for a
locus carries the same parsed hierarchy (broadcast from the locus's residue
path).

### The mechanism payload

A locus is **gated** at exactly one scale (its settled resolution ℓ̂*). The
mechanism payload is attached **only to that gated row**; all other scale rows
for the locus carry `null` in the `mech_*` columns. This is why those columns are
object-typed — they legitimately hold NA off the gated row.

A mechanism with `direction == "mixed"` carries `null` signed-beta fields
(`mech_beta_signed`, `mech_beta_ci_lower`, `mech_beta_ci_upper`, `mech_beta_se`)
— there is no coherent signed effect. S0 validates this "null iff mixed" rule.

### Invariants S0 enforces

- exactly one gated row per locus;
- gated profile `rho` equals the mechanism `rho`;
- each mechanism's `loci` all gate to that mechanism's `region_id` and scale;
- the mechanism↔gated-locus mapping is an exact partition (no orphan mechanisms,
  no orphan gated loci);
- `(serotype, canon_label, scale_level)` is unique.

---

## Relationship between the tables

The two tables share the `canon_label` join key at residue scale, but are kept
separate on purpose. The replicate table holds *raw per-run* effects (`r`); the
STRIDE table holds the *reproducibility-summarised* quantities (`rho`, `beta`,
`direction`). A later stage that wants to relate them can join on
`(serotype, canon_label)` — S0 does not do this join, and does not collapse the
levels.
