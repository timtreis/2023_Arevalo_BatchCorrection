# F4 Reference Mapping — Wave 2 Query Experiment

**Reference**: scenario_5, all 5 sources (2, 3, 6, 8, 10), COMPOUND + TARGET2  
**Query**: one wave2 source at a time (source_5, source_9, source_11)  
**Methods**: Symphony (linear), scPoli (non-linear)  
**Ablation**: T+ / T- / T-_matched per query source

---

## The experiment in one sentence

Train a full 5-source reference atlas on scenario_5. Independently map each wave2 source
into it under three query-arm conditions (T+/T-/T-_matched). Evaluate cross-source compound
retrieval precision in the joint embedding.

---

## What changes from the smoke test

| Dimension | Smoke test | F4 experiment |
|---|---|---|
| Reference sources | 4 (excl. source_8) | All 5 |
| Query source | source_8 (in ref training data) | wave2 source (never seen) |
| Compound overlap | ~100% | ~5,392 shared compounds |
| Atlas name convention | `..._excl_{source}` | `..._scenario5_full` |
| labeled_indices | all cells (default) | shared-compound cells only |
| Evaluation | within-query compound ASW (wrong) | cross-source retrieval precision |

---

## Phase 0 — Prerequisite: run wave2 through the Snakemake pipeline

**Blocker — nothing else can proceed without this.**

The pipeline needs to run scenario_wave2 to produce:

```
outputs/scenario_wave2/mad_int_featselect.parquet   ← query input for nb01
```

The scenario config already exists at `inputs/conf/scenario_wave2.json`:
- sources: source_5, source_9, source_11
- plate_types: COMPOUND + TARGET2
- same IQR/clip params as scenario_5

Steps:
1. Confirm pipeline rules support scenario_wave2 (no scenario-specific hard-coding).
2. Launch: `pixi run pipeline --config scenario=scenario_wave2` (or equivalent Snakemake invocation).
3. The optimization step (Optuna HPO) is NOT needed for the query scenario — only preproc through `mad_int_featselect` is needed. If possible, run only up to that rule.

**Note**: wave2 needs only the preprocessing rules, not training. Check whether there's a Snakemake target for just preprocessing.

---

## Phase 1 — Full reference atlas (no source holdout)

### Design decision

The atlas is trained once on all 5 scenario_5 sources. It is then reused for all three wave2
query sources. We do NOT do a separate leave-one-out model per query source — the reference is
fixed, and only the query-side condition embeddings are learned during mapping.

### 1a — Update nb20 (scPoli atlas)

**Current**: splits out `QUERY_SOURCE` from scenario_5.  
**New**: when `QUERY_SOURCE = ""` (empty string), use all sources.

Parameter block changes:
```python
QUERY_SOURCE = ""          # empty = use all sources; non-empty = exclude that source
REFERENCE_SCENARIO = "scenario_5"
ATLAS_NAME = (
    f"scpoli_atlas_{REFERENCE_SCENARIO}_full"
    if not QUERY_SOURCE
    else f"scpoli_atlas_excl_{QUERY_SOURCE}"
)
```

Path 1 (pipeline model symlink) is unchanged — will still be missing (scenario_5 was
run before model-saving was added), so it falls through to inline training.

Path 2 (inline training) changes:
```python
# Before: adata = ...; ref, _ = split_by_source(adata, SOURCE_COL, QUERY_SOURCE)
# After:
adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
ref = adata if not QUERY_SOURCE else split_by_source(adata, SOURCE_COL, QUERY_SOURCE)[0]
```

Everything else (HPO params load, model init, train, save) is unchanged.

Expected training time: somewhat longer than the 4-source smoke test (~20% more data).
Use the same `n_epochs=600`, `use_early_stopping=True` as the smoke test.

### 1b — Update nb30 (Symphony atlas)

**Current**: reads a pre-split h5ad from `data/query_arms/{QUERY_SOURCE}_reference.h5ad`.  
**New**: when `QUERY_SOURCE = ""`, reads from the scenario_5 parquet directly.

Parameter block changes:
```python
QUERY_SOURCE = ""
REFERENCE_SCENARIO = "scenario_5"
BATCH_KEY = "Metadata_Source"
N_PCS = 50

ATLAS_DIR = (
    MODEL_OUT / f"symphony_atlas_{REFERENCE_SCENARIO}_full"
    if not QUERY_SOURCE
    else MODEL_OUT / f"symphony_atlas_excl_{QUERY_SOURCE}"
)
ATLAS_DIR.mkdir(parents=True, exist_ok=True)
```

Data loading changes:
```python
# Before: ref = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad")
# After:
if QUERY_SOURCE:
    ref = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad")
else:
    from src.data_io import load_parquet_as_anndata
    from src.paths import scenario_input_parquet
    ref = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
```

Scale + PCA + Harmony cell is unchanged.

---

## Phase 2 — Query arm preparation (wave2 sources)

### Update nb01

Current nb01 takes a source from scenario_5. New version takes a source from any scenario.

New parameters:
```python
QUERY_SOURCE = "source_5"         # which wave2 source to prepare
QUERY_SCENARIO = "scenario_wave2" # which scenario parquet to slice it from
```

The query arm construction logic stays the same (make_tplus / make_tminus / make_tminus_matched),
but the input parquet changes:
```python
# Before: adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
# After:
from src.paths import scenario_input_parquet
adata = load_parquet_as_anndata(scenario_input_parquet(QUERY_SCENARIO))
query = adata[adata.obs["Metadata_Source"] == QUERY_SOURCE].copy()

# Reference arm: scenario_5 full (all 5 sources, for context in Symphony atlas)
# Note: only needed if Symphony nb30 reads from query_arms; for the full atlas it reads
# the parquet directly, so the reference h5ad is only needed for the excl. variant.
ref = load_parquet_as_anndata(scenario_input_parquet("scenario_5"))
```

Outputs written to `data/query_arms/`:
- `{QUERY_SOURCE}_tplus.h5ad`
- `{QUERY_SOURCE}_tminus.h5ad`
- `{QUERY_SOURCE}_tminus_matched.h5ad`
- `{QUERY_SOURCE}_swap_ids.json`
- `{QUERY_SOURCE}_reference.h5ad`  ← full scenario_5 (for any variant that still needs it)

Run nb01 three times (once per wave2 source: source_5, source_9, source_11), OR parameterize
with a sweep.

---

## Phase 3 — Mapping

### 3a — scPoli mapping (nb21)

Key change: pass `labeled_indices` = indices of wells whose compound exists in the reference.

```python
QUERY_SOURCE = "source_5"         # wave2 source
T_ARM = "tplus"
REFERENCE_ATLAS = "scenario_5_full"
ATLAS_NAME = f"scpoli_atlas_{REFERENCE_ATLAS}"
N_EPOCHS = 50
PRETRAIN_EPOCHS = 0   # atlas encoder already pretrained; skip pretraining
```

**Why `PRETRAIN_EPOCHS = 0`**: The encoder is frozen. Pretraining epochs do nothing useful
(no prototype loss fires, nothing to learn). Setting to 0 saves wall-clock time. Only the
condition embedding epochs matter.

**Computing `labeled_indices`**:
```python
# shared_compounds = compounds present in both reference training data and query
ref_compounds = set(
    pd.read_parquet(
        scenario_input_parquet("scenario_5"),
        columns=["Metadata_JCP2022"]
    )["Metadata_JCP2022"].unique()
)
labeled_mask = query.obs["Metadata_JCP2022"].isin(ref_compounds)
labeled_indices = list(np.where(labeled_mask)[0])
```

Passing to `load_query_data`:
```python
model_query = scPoli.load_query_data(
    adata=query,
    reference_model=str(atlas_dir),
    labeled_indices=labeled_indices,
)
```

This focuses the prototype pull loss on the ~5,392 shared compound wells. Wave2-only compound
wells still flow through the frozen encoder and get the batch embedding correction, but they
don't generate noisy prototype gradients.

### 3b — Symphony mapping (nb31)

No structural change needed — Symphony is a linear projection, it handles any query without
knowing which compounds are shared. The atlas name changes to `symphony_atlas_scenario5_full`.

---

## Phase 4 — Evaluation (nb40 redesign)

### The singleton problem

98% of wave2 query wells are singleton per compound (1 well/compound/source). Within-query
compound ASW is undefined for singletons and measures nothing useful.

### New primary metric: cross-source compound precision@k

For each query well, find its k nearest neighbors in the **reference** embedding. Compute
the fraction where the neighbor has the same compound ID. Average over all query wells.

```python
# In nb40, after loading the mapped query and reference:
from sklearn.neighbors import NearestNeighbors

ref_emb = ref.obsm[ref_key]       # (n_ref_wells, latent_dim)
query_emb = eval_query.obsm[EMBEDDING_KEY]  # (n_eval_wells, latent_dim)

nbrs = NearestNeighbors(n_neighbors=K, metric="cosine").fit(ref_emb)
dists, indices = nbrs.kneighbors(query_emb)

ref_compounds = ref.obs[COMPOUND_COL].values
query_compounds = eval_query.obs[COMPOUND_COL].values

# For each query well: fraction of top-K ref neighbors with same compound
match = (ref_compounds[indices] == query_compounds[:, None])  # (n_eval, K)
precision_at_k = match.mean()  # mean fraction across all wells and all k positions
```

Choose `K = 10` initially. Report also `precision@1` (nearest neighbor accuracy) as a
stricter check.

This metric works for singleton query wells (only 1 query well per compound is fine — we're
looking up its neighbors in the reference which has 3-64 wells per compound).

### Secondary metrics (retained from original nb40)
- `batch_asw` on joint embedding: cross-source mixing within compound class
- `kbet` and `ilisi` on joint embedding via pynndescent

### Drop
- `compound_asw` on query-alone (meaningless for singletons)

---

## Phase 5 — Run order and papermill sweep

For each wave2 query source (source_5, source_9, source_11):
  For each T-arm (tplus, tminus, tminus_matched):
    For each method (symphony, scpoli):
      Run nb40 → appends one row to `results/f4/{QUERY_SOURCE}_metrics.csv`

That's 3 × 3 × 2 = 18 nb40 runs. Before any of those:
- 1 nb20 run (full scPoli atlas) — slow, GPU required
- 1 nb30 run (full Symphony atlas) — fast
- 3 nb01 runs (one per wave2 source)
- 3 × 3 = 9 nb21 runs (scPoli mapping, one per source × arm)
- 3 × 3 = 9 nb31 runs (Symphony mapping, one per source × arm)

Total: 1 + 1 + 3 + 9 + 9 + 18 = 41 notebook executions after phase 0 (wave2 preprocessing).

---

## Files to modify

| File | Change |
|---|---|
| `notebooks/20_train_atlas_scpoli.py` | Add `QUERY_SOURCE = ""` mode; conditional split |
| `notebooks/30_train_atlas_harmony_symphony.py` | Add `QUERY_SOURCE = ""` mode; parquet reading path |
| `notebooks/01_prepare_query_arms.py` | Add `QUERY_SCENARIO` param; slice by `QUERY_SOURCE` |
| `notebooks/21_map_scpoli.py` | Add `labeled_indices` for shared compounds; `PRETRAIN_EPOCHS=0` |
| `notebooks/40_metrics.py` | Replace `compound_asw` with cross-source precision@k |

No new files needed. `src/paths.py`, `src/metrics.py`, `src/data_io.py` need no changes.

---

## Open questions before implementation

1. **Does the pipeline Snakemake support scenario_wave2 without modification?** Check whether
   any rule hard-codes scenario names or makes assumptions about source lists.

2. **scPoli memory with full 5-source reference**: scenario_5 has 370K wells. Memory budget
   for scPoli training needs checking — the 4-source smoke test used 296K wells. 25% more.

3. **Evaluation: which compound set to use for labeled_indices in nb21?**
   - Option A: all compounds in scenario_5 (any compound that appears in the reference)
   - Option B: only compounds appearing in ≥2 reference sources (higher-quality prototypes)
   - Start with Option A; Option B is a sensitivity check.

4. **Symphony N_PCS for full 5-source atlas**: 50 PCs is fine; may want to verify the
   variance explained doesn't change dramatically with the 5th source added.

5. **Should we run the T+ ablation per-source or jointly?**
   Start with one source (source_5) for end-to-end validation, then sweep all three.
