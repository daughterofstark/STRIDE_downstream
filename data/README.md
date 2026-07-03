# `data/` — your datasets go here (not committed)

This directory is where you place **your own STRIDE outputs**. Its contents are
**git-ignored** (except this README): real datasets are user-supplied inputs, not
part of the source. The repository works after cloning without any data here —
run against `examples/small_synthetic_dataset/` to see it work.

## Three data levels

The framework recognises three data levels. **Level 0 is never read** by this
repository.

| Level | What | Committed? | Read by S0? |
|---|---|---|---|
| **0 — Raw MD inputs** | `xtc`, `gro`, `tpr`, … | no | **no** (never a dependency) |
| **1 — Replicate observations** | `*_correlations_v5.csv` per replicate | no | yes (primary) |
| **2 — STRIDE summaries** | `*_profile.csv`, `*_mechanism.json` | no | yes (derived) |

## Expected layout

The discovery logic accepts a nested run-directory layout and locates the
Level-2 summaries either in a `summaries/` subfolder or at the data root. A
complete dataset looks like:

```
data/                                  (or any path you pass to --data-root)
├── 1st_run/
│   ├── DENV1/analysis_output/DENV1_correlations_v5.csv
│   ├── DENV2/analysis_output/DENV2_correlations_v5.csv
│   └── ...
├── 2nd_run/
│   └── ... (same serotypes)
├── 3rd_run/
│   └── ...
└── summaries/
    ├── DENV1_profile.csv
    ├── DENV1_mechanism.json
    ├── DENV2_profile.csv
    └── ...
```

- Run directories define replicate order (`1st_run` → replicate index 1). The
  names are discovered, not hardcoded — any sorted set of run dirs works.
- Each serotype must appear under every run directory (equal replicate counts),
  unless you pass `--allow-unequal-replicates`.
- Level-2 summaries may live in `summaries/` **or** directly at the data root.

## Running S0 on your data

```bash
stride-s0 --data-root data --output-dir outputs
# or
python -m stride_analysis --data-root data --output-dir outputs
```

Partial datasets are supported:

```bash
# summaries only (no replicate correlations yet):
stride-s0 --data-root data --output-dir outputs --no-require-replicates

# replicates only (no STRIDE summaries yet):
stride-s0 --data-root data --output-dir outputs --no-require-summaries
```

See `../examples/small_synthetic_dataset/` for a tiny, valid dataset in exactly
this layout.
