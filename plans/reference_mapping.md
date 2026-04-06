# Reference Mapping Experiments Plan

**Date**: 2026-04-05
**Goal**: Evaluate reference mapping as a practical deployment mode for Cell Painting batch correction. Train a reference atlas on Wave 1 data, map Wave 2 query data into it, measure whether compound identity and biological signal are preserved.

**Paper angle**: Reference mapping answers the practical question: "I have a new experiment from a different lab/microscope. Can I map it into an existing atlas and get meaningful results without retraining?"

---

## 0. Theoretical Foundations

### 0.1 Three Paradigms, Increasing Complexity

| | Symphony | scVI | scPoli |
|--|---------|------|--------|
| **Model class** | Linear mixture of experts | Deep generative (VAE) | Deep generative (CVAE) + prototypes |
| **Batch effect assumption** | Piecewise linear (per-cluster shift) | Nonlinear (neural net decoder) | Nonlinear + supervised anchors |
| **Likelihood** | Implicit (PCA + linear) | ZINB/NB/Gaussian (configurable) | MSE (Gaussian) |
| **Correct for Cell Painting?** | Partially (linear approx) | Yes if `gene_likelihood="normal"` | Yes (MSE = Gaussian) |
| **Surgery mechanism** | Project + assign to ref clusters | Add new batch embedding column | Add new condition embedding vector |
| **What is frozen** | PCA loadings + cluster params | All weights except new batch cols | All weights except new embedding |
| **Query parameters** | ~K×d correction factors | ~2×n_hidden new weights | 1 embedding vector (d_emb dims) |
| **Mapping speed** | ~1-2 seconds (CPU) | ~minutes (GPU) | ~minutes (GPU) |
| **Supervised signal** | None | None (unsupervised) | Prototype loss on labels |
| **TARGET2 anchor benefit** | Moderate (better clusters) | Modest (implicit) | Large (explicit prototype loss) |

### 0.2 Mathematical Models

**Symphony/Harmony**: Operates in PCA space. Alternates between (1) soft k-means clustering with a maximum-diversity penalty (information-theoretic, penalizes clusters dominated by a single batch via KL divergence) and (2) per-cluster linear correction: `Z_corrected = Z - sum_k R_ik * (Y_k_batch - Y_k)`. The correction is piecewise-linear — within each soft cluster, batch effects are simple shifts; across clusters, different shifts are allowed. Symphony compresses the reference into PCA loadings + cluster centroids + correction factors (~100KB). Query mapping is just project → assign → correct.

**scVI**: Generative model `p(x|z,s)` with latent `z ~ N(0,I)` and batch-conditional neural network decoder. The encoder `q(z|x,s)` maps cells to latent space conditioned on batch. "Surgery" = freeze all weights, add a new column to the first encoder/decoder layers for the query batch embedding, fine-tune only those new parameters (50-200 epochs). **Critical for Cell Painting**: default likelihood is ZINB (count data) — should test `gene_likelihood="normal"` which is the correct Gaussian likelihood for continuous CellProfiler features. The current `adata.X -= min_value` shift to non-negative is a hack around the ZINB assumption.

**scPoli**: CVAE with two key innovations: (1) **learnable condition embeddings** of fixed dimension d_emb (instead of one-hot batch vectors), enabling addition of new batches by adding a single embedding vector, and (2) **cell-type prototypes** — learned anchor points in latent space, one per label, with loss `L_proto = sum_k sum_{i in class_k} d(z_i, p_k)`. Total loss: `L = L_recon + beta*L_KL + eta*L_proto`. Uses MSE reconstruction (correct Gaussian likelihood for Cell Painting). Surgery = freeze all weights, initialize new embedding as mean of existing embeddings, optionally fine-tune.

### 0.3 When Each Method is Optimal

**Linear/additive batch effects** (intensity shifts, background differences):
→ **Symphony wins.** Piecewise-linear model is exactly right. No wasted capacity. Fastest, most stable.

**Non-linear, feature-dependent batch effects** (confocal vs widefield PSF differences distort texture/shape features non-linearly):
→ **scVI or scPoli wins.** Neural network decoders can learn arbitrary non-linear batch-conditional transformations.

**Batch effects that interact with biology** (compound-specific distortions across microscopes):
→ **scPoli wins.** Prototype loss explicitly anchors biological identity across batches. Without it, the model has no incentive to preserve compound-specific structure.

**Cell Painting reality**: Likely a mixture of all three. Intensity features ≈ linear. Texture/shape features (Haralick, Zernike, granularity) ≈ non-linear. Compound-specific distortions across different microscope PSFs ≈ biology-interacting.

**Prediction for F4**: scPoli fine-tuned > scVI online > Symphony > scPoli zero-shot > PCA-only. Gap between Symphony and deep learning largest for texture/shape features. CV8000 query (seen microscope) maps better than Operetta query (unseen).

### 0.4 scVI Likelihood Issue

The current pipeline uses ZINB (default) for scVI, but Cell Painting features are continuous, can be negative, and are not counts. scVI supports `gene_likelihood="normal"` which is theoretically correct. The current `adata.X -= min_value` shift is a workaround. **Action**: Test `gene_likelihood="normal"` vs ZINB as part of the reference mapping experiments. scPoli already uses MSE (Gaussian) correctly.

---

## 1. Experiment Design: F4 (Wave 1 → Wave 2)

### 1.1 Why F4

F4 is the best-powered reference mapping test in JUMP-CP:
- **5,392 evaluable compounds** (5,086 COMPOUND + 306 TARGET2) shared between train and test
- **Cross-wave generalization**: completely different production campaigns
- **Microscope heterogeneity**: train has CV8000+Opera Phenix+ImageXpress; test has CV8000+Opera Phenix+Operetta (partial overlap)
- **Practical utility proof**: exactly the workflow a real lab would use

### 1.2 Reference (Train): C4

| Property | Value |
|----------|-------|
| Sources | source_2, source_3, source_6, source_8, source_10, source_15 |
| Plate types | COMPOUND + TARGET2 (384-well only, excludes source_1's 1536-well) |
| Batch key | Metadata_Source |
| Microscopes | CV8000 (S2,S6,S10), Opera Phenix (S3,S8), ImageXpress (S15) |
| Compounds | ~82K COMPOUND + 302 TARGET2 |
| Estimated cells | ~500K+ |
| Label key | Metadata_JCP2022 |

Config to create: `inputs/conf/scenario_c4.json`

### 1.3 Query (Test): Wave 2

| Property | Value |
|----------|-------|
| Sources | source_5, source_9, source_11 |
| Plate types | COMPOUND + TARGET2 |
| Batch key | Metadata_Source |
| Microscopes | CV8000 (S5), Operetta (S9), Opera Phenix (S11) |
| Compounds | ~36K COMPOUND + 302 TARGET2 |
| Estimated cells | ~150-200K |
| Label key | Metadata_JCP2022 |

Config exists: `inputs/conf/scenario_wave2.json`

### 1.4 Evaluation Set

- **Shared compounds**: 5,086 COMPOUND + 306 TARGET2 = 5,392 total
- **14.1% of Wave 2 test compounds** appear in Wave 1 train set
- Evaluation is performed ONLY on shared compounds (fair comparison)
- Both negcon mAP (compound-level) and MOA mAP (mechanism-level) where DRH annotations exist

---

## 2. Methods

### 2.1 scPoli / scArches — Transfer Learning (No Retraining)

**How it works**: scPoli learns disentangled batch and biology representations. Query data is embedded into the reference's learned latent space without retraining the model. The model's learned batch effects are applied to the new batch.

**Implementation**:
```python
# Reference training (already in pipeline)
model = scPoli(adata_ref, condition_keys=["Metadata_Source"], cell_type_keys=["Metadata_JCP2022"])
model.train(n_epochs=400, pretraining_epochs=pretrain_epochs, early_stopping=True)
model.save("reference_model/")

# Query mapping (NEW)
model = scPoli.load("reference_model/", adata_query)
query_latent = model.get_latent(adata_query, mean=True)
```

**Key decisions**:
- Use best HPO params from C4 scenario (or from S4/S5 if C4 not ready)
- `mean=True` for deterministic embeddings (no sampling noise)
- No fine-tuning on query — pure zero-shot transfer

**Script**: `scripts/map_with_scpoli.py` (new)

**Variants to test**:
1. **Zero-shot**: Load reference model, embed query directly
2. **Fine-tuned**: Load reference model, run 10-20 additional epochs on query data with frozen encoder (scArches surgery)

### 2.2 scVI — Online Update via scArches

**How it works**: scVI learns a VAE with batch-conditional decoder. Reference mapping uses `load_query_data()` which adds a new batch embedding for the query batch and fine-tunes for a few epochs while keeping most weights frozen.

**Implementation**:
```python
# Reference training (already in pipeline)
scvi.model.SCVI.setup_anndata(adata_ref, batch_key="Metadata_Source")
vae = scvi.model.SCVI(adata_ref, n_layers=2, n_latent=30)
vae.train(max_epochs=999999, early_stopping=True)
vae.save("reference_model/")

# Query mapping (NEW)
scvi.model.SCVI.prepare_query_anndata(adata_query, "reference_model/")
query_model = scvi.model.SCVI.load_query_data(adata_query, "reference_model/")
query_model.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})
query_latent = query_model.get_latent_representation(adata_query)
```

**Key decisions**:
- Use `scvi_single` architecture (single batch key = Metadata_Source)
- Online update epochs: test 50, 100, 200 (hyperparameter)
- Weight decay 0.0 during fine-tuning (prevent catastrophic forgetting)
- Use best HPO architecture from C4 (n_layers, n_latent, n_hidden, dropout, gene_likelihood)

**Script**: `scripts/map_with_scvi.py` (new)

**Variants to test**:
1. **Online update (default)**: Fine-tune on query with frozen reference weights
2. **Zero-shot**: Skip fine-tuning, just encode query through reference encoder (baseline)

### 2.3 Symphony / Harmony — Linear Reference Compression

**How it works**: Symphony compresses a Harmony-integrated reference into summary statistics (PCA loadings + Harmony mixture model parameters). Query cells are projected into reference PCA space, then batch-corrected using the stored Harmony model. Purely linear — no deep learning.

**Implementation**:
```python
# Reference construction
import harmonypy
import symphonypy  # or implement manually

# Step 1: PCA on reference
pca = PCA(n_components=50).fit(X_ref)
Z_ref = pca.transform(X_ref)

# Step 2: Harmony on reference PCA
ho = harmonypy.run_harmony(Z_ref, meta_ref, ["Metadata_Source"])
# Save: pca.components_, ho.R (cluster assignments), ho.Z_corr, ho.Phi_moe

# Query mapping
# Step 1: Project query into reference PCA space
Z_query = pca.transform(X_query)

# Step 2: Apply Harmony correction using stored reference model
# symphonypy.map_embedding(Z_query, reference_object)
Z_query_corrected = symphony_map(Z_query, reference_harmony_model)
```

**Key decisions**:
- PCA components: 50 (match harmony_v1 HPO range)
- Use harmony_v1 (standard) not harmony_v2 (harmonypy fork)
- Symphony is an R package (`symphony`) — check if Python wrapper exists, otherwise use `symphonypy` or implement the linear algebra directly
- Alternative: manually re-implement the projection (it's just PCA + ridge regression, ~50 lines)

**Script**: `scripts/map_with_symphony.py` (new)

**Dependency check needed**: Is `symphonypy` available? If not, the mapping is:
```python
# Manual Symphony mapping (no package needed):
Z_query_pca = pca.transform(X_query)
# Assign query cells to Harmony clusters using reference centroids
cluster_assignments = assign_clusters(Z_query_pca, ref_centroids)
# Apply per-cluster corrections learned from reference
Z_query_corrected = Z_query_pca - cluster_corrections[cluster_assignments]
```

### 2.4 Baselines (for comparison)

1. **No correction**: Project query into reference PCA space, no batch correction
2. **Naive Harmony**: Run Harmony on combined ref+query (NOT a reference mapping — requires full retraining). This is the upper bound for what Harmony could achieve if you had access to all data.
3. **Combat**: Apply reference Combat parameters to query (linear correction)

---

## 3. Evaluation Metrics

### 3.1 Primary: Compound Reproducibility (mAP)

Already in pipeline (`metrics/map.py`). Compute on shared compounds only.

| Metric | Description | Good value |
|--------|-------------|------------|
| `negcon_mean_map` | mAP: do same-compound replicates cluster? | >0.3 |
| `negcon_frac_p` | Fraction of compounds with significant mAP | >0.5 |

Compute separately for:
- Query-only mAP: How well do query replicates cluster in mapped space?
- Cross-set mAP: Do query cells land near reference cells of the same compound?

### 3.2 Secondary: MOA Preservation

| Metric | Description | Good value |
|--------|-------------|------------|
| `moa_mean_map` | mAP on Metadata_DRH_MOA grouping | >0.2 |
| `moa_frac_p` | Fraction of MOAs with significant mAP | >0.3 |

### 3.3 Batch Integration Quality

From scib-metrics (already in pipeline):

| Metric | What it measures |
|--------|-----------------|
| `iLISI` | Are ref+query cells interleaved in neighborhoods? |
| `graph_connectivity` | Can you traverse from query to reference? |
| `pcr_comparison` | Is batch variance removed? |
| `kBET` | Local batch mixing |

### 3.4 Reference Mapping-Specific

| Metric | Description | Implementation |
|--------|-------------|----------------|
| Cross-set kNN accuracy | For shared compounds: what % of query cells have a same-compound reference cell in their k nearest neighbors? | Custom, ~20 lines |
| Procrustes disparity | Geometric alignment between ref and query embeddings of shared compounds | `scipy.spatial.procrustes` |
| Label transfer accuracy | Train kNN classifier on reference labels, predict query labels, measure accuracy on shared compounds | `sklearn.neighbors.KNeighborsClassifier` |

---

## 4. Pipeline Integration

### 4.1 New Scripts

| Script | Purpose | Container |
|--------|---------|-----------|
| `scripts/map_with_scpoli.py` | scPoli zero-shot + fine-tuned mapping | scpoli.sif |
| `scripts/map_with_scvi.py` | scVI online update mapping | scvi.sif |
| `scripts/map_with_symphony.py` | Symphony linear projection | base.sif or harmony_v1.sif |
| `scripts/eval_reference_mapping.py` | Cross-set mAP, kNN accuracy, Procrustes | scibmetrics-gpu.sif or base.sif |

### 4.2 New Snakemake Rules

```python
# In rules/projection.smk or new rules/reference_mapping.smk

rule train_reference:
    """Train reference model on C4 data."""
    input: adata=f"outputs/scenario_c4/{preproc}_all_methods.h5ad"
    output: model_dir=directory(f"outputs/f4_reference_mapping/models/{{method}}/")
    # For each method, save trained model checkpoint

rule map_query:
    """Map Wave 2 query data into reference space."""
    input:
        query=f"outputs/scenario_wave2/{preproc}.parquet",
        model=f"outputs/f4_reference_mapping/models/{{method}}/"
    output: f"outputs/f4_reference_mapping/embeddings/{{method}}_query.parquet"

rule eval_reference_mapping:
    """Evaluate mapping quality on shared compounds."""
    input:
        ref_emb=f"outputs/f4_reference_mapping/embeddings/{{method}}_reference.parquet",
        query_emb=f"outputs/f4_reference_mapping/embeddings/{{method}}_query.parquet"
    output: f"outputs/f4_reference_mapping/metrics/{{method}}_metrics.parquet"
```

### 4.3 Output Structure

```
outputs/f4_reference_mapping/
    models/
        scpoli/          # Saved scPoli reference model
        scvi/            # Saved scVI reference model
        symphony/        # Saved PCA + Harmony parameters
    embeddings/
        scpoli_reference.parquet
        scpoli_query_zeroshot.parquet
        scpoli_query_finetuned.parquet
        scvi_reference.parquet
        scvi_query_online.parquet
        scvi_query_zeroshot.parquet
        symphony_reference.parquet
        symphony_query.parquet
        baseline_pca_query.parquet
    metrics/
        all_methods_map.parquet       # mAP per method
        all_methods_moa_map.parquet   # MOA mAP per method
        all_methods_scib.parquet      # iLISI, graph_conn, etc.
        all_methods_crossset.parquet  # kNN accuracy, Procrustes
        all_methods_summary.parquet   # Pivoted summary table
    plots/
        reference_mapping_comparison.pdf
        umap_per_method.pdf
        map_by_microscope.pdf         # Breakdown by query microscope
```

---

## 5. Execution Plan

### Phase 0: Prerequisites
- [ ] C4 scenario config created and run (Wave 1 reference atlas)
- [ ] Wave 2 scenario run (needed for preprocessing, already have config)
- [ ] Verify compound overlap: load both preprocessed parquets, intersect Metadata_JCP2022

### Phase 1: Reference Model Training (~2 days)
- [ ] Train scPoli on C4 data, save model checkpoint
- [ ] Train scVI on C4 data, save model checkpoint
- [ ] Run Harmony on C4, save PCA loadings + Harmony parameters
- [ ] Verify reference models by computing within-reference mAP

### Phase 2: Query Mapping (~1-2 days)
- [ ] scPoli zero-shot: load model, embed Wave 2 query
- [ ] scPoli fine-tuned: load model, 10-20 epochs on query, embed
- [ ] scVI online update: load_query_data, fine-tune 50/100/200 epochs, embed
- [ ] scVI zero-shot: encode query through reference encoder (no update)
- [ ] Symphony: project query into reference Harmony space
- [ ] Baselines: PCA-only, naive Harmony (combined ref+query)

### Phase 3: Evaluation (~1-2 days)
- [ ] Compute negcon mAP on shared compounds (query-only + cross-set)
- [ ] Compute MOA mAP where DRH annotations exist
- [ ] Compute scib-metrics on combined ref+query embeddings
- [ ] Compute cross-set kNN accuracy and Procrustes
- [ ] Breakdown by query microscope type (CV8000 vs Operetta vs Opera Phenix)

### Phase 4: Figures (~1 day)
- [ ] Method comparison bar chart (mAP per method)
- [ ] UMAP per method (ref colored by source, query overlaid)
- [ ] mAP by microscope heatmap (method × microscope)
- [ ] Scatter: reference mAP vs query mAP per compound (shows which compounds transfer)

---

## 6. TARGET2 Anchor Ablation

### 6.1 Theoretical Basis

TARGET2 compounds (302 shared across all sources) act as "spike-in controls" — known biological identities that should produce the same morphological signature regardless of lab/microscope. The principle: **if two cells from different batches have the same treatment, any remaining feature difference is attributable to batch effects.**

Each method benefits differently:

**Harmony/Symphony**: TARGET2 compounds naturally satisfy the maximum-diversity clustering constraint (they exist in every batch). They provide ground-truth examples of "same biology, different batch" that constrain the per-cluster correction factors. In effect, they are **calibration points** for the linear correction.

**scVI**: The decoder learns `p(x|z,s)`. When the model sees the same compound across multiple batches, it must explain differences as batch effects (s) not biology (z). TARGET2 provides ~40K labeled cells of explicit paired examples — same z, different s — making disentanglement well-posed. This is implicit (scVI doesn't use labels), but the data structure helps.

**scPoli**: Benefits most directly. Each TARGET2 compound gets a prototype in latent space. The prototype loss pulls all replicates across all batches toward the same point. With 302 compounds × ~21 reps × 6-7 sources = ~40K supervised cross-batch anchor cells. This is **explicit supervised disentanglement signal**.

### 6.2 Ablation Design

| Condition | Reference data | Query data | Evaluation |
|-----------|---------------|------------|------------|
| **Full (with TARGET2)** | C4: COMPOUND + TARGET2 | Wave 2: COMPOUND + TARGET2 | mAP on ~5,086 shared COMPOUND (TARGET2 excluded from eval) |
| **No TARGET2** | C4: COMPOUND only | Wave 2: COMPOUND only | mAP on same ~5,086 shared COMPOUND |
| **TARGET2 in ref only** | C4: COMPOUND + TARGET2 | Wave 2: COMPOUND only | mAP on same ~5,086 shared COMPOUND |

Key: evaluation is always on shared COMPOUND compounds only — TARGET2 is the treatment variable, not the evaluation variable.

### 6.3 Expected Outcomes

- **scPoli**: Largest improvement with TARGET2 (direct prototype supervision)
- **scVI**: Modest improvement (implicit disentanglement help)
- **Symphony**: Moderate improvement (better-constrained clusters/correction)
- **"Ref only" condition**: If reference-only TARGET2 still helps, it proves the anchors improve the learned model itself, not just the mapping step — important practical finding (query labs don't need to run TARGET2 plates)

### 6.4 Implementation

Two additional reference models per method (3 methods × 3 conditions = 9 total model trainings). Preprocessing already supports filtering by plate type via the config's `plate_types` field.

---

## 7. Risks and Mitigations

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| C4 too large to train (500K+ cells) | Medium | Use scVI subsampling; scPoli handles large datasets well; Harmony is linear |
| Symphony Python package unavailable | High | Implement manually (~50 lines of linear algebra) |
| All methods perform similarly | Medium | Include baselines (PCA-only) to show methods add value |
| Operetta microscope (S9) too different | Medium | Report per-microscope breakdown; if Operetta fails, it's a finding |
| Wave 2 source_9 1536-well format issues | Medium | Can exclude S9 as sensitivity analysis |
| scVI online update overfits on small query batches | Low | Test multiple epoch counts; compare to zero-shot |

---

## 8. Paper Presentation Strategy

### 8.1 Narrative Arc

1. **Problem**: Consortium-scale Cell Painting generates data across many labs and microscopes. As the atlas grows, retraining from scratch becomes impractical. New labs need to map data into the existing atlas.
2. **Why reference mapping**: Unlike de novo integration, reference mapping (a) scales to unlimited query data, (b) preserves the reference embedding, (c) doesn't require sharing raw reference data (privacy/IP), (d) enables a "submit and map" workflow.
3. **Three paradigms**: Linear compression (Symphony), deep generative surgery (scVI), prototype-guided transfer (scPoli) — representing increasing complexity and decreasing assumptions.
4. **Experiment F4**: Train on Wave 1 (6 sources, 3 microscope types), map Wave 2 (3 sources, 3 microscope types including one unseen). Evaluate on 5,392 shared compounds — **strongest evaluation in the Cell Painting reference mapping literature** (most benchmarks use cell-type labels, we have compound identity ground truth).
5. **TARGET2 ablation**: Do shared control compounds ("anchors") improve reference mapping quality? Which methods benefit most?
6. **Results + practical implications**: Which method should a new JUMP-CP lab use? When is linear sufficient? When do you need deep learning?

### 8.2 Recommended Figures

1. **Method comparison bar chart** (main result): mAP (negcon + MOA) per method, with PCA-only baseline. Error bars from compound-level variance.
2. **UMAP panels** (2×3 grid): Reference in gray, query colored by source. One panel per method. Shows qualitative integration.
3. **Microscope-stratified heatmap**: Methods (rows) × query microscope types (columns), colored by mAP. Reveals which method handles unseen microscopes best.
4. **TARGET2 ablation**: Grouped bar chart showing mAP with/without TARGET2 anchors, per method. The differential reveals which paradigm benefits most from supervision.
5. **Compound-level scatter**: Reference mAP (x) vs mapped query mAP (y) per compound. Points near diagonal = successful transfer. Deviations reveal which compounds/morphologies fail.
6. **Computational comparison**: Log-scale bar chart of mapping wall-clock time. Symphony at seconds vs deep learning at minutes — practically important.

### 8.3 Key Differentiator vs Prior Work

- Symphony paper (Kang 2021) evaluates on cell-type label transfer accuracy
- scArches paper (Lotfollahi 2022) evaluates on kNN classification + LISI
- **Our paper**: 5,392 compounds with ground-truth identity across train/test — orders of magnitude more evaluation power than cell-type labels
- Also: first systematic comparison of reference mapping paradigms on morphological profiling data (not scRNA-seq)
- Note CellPainTR (arxiv 2509.06986) as concurrent work — Transformer approach, different paradigm, doesn't do reference mapping per se

---

## 9. Expected Results (Hypotheses)

1. **scPoli fine-tuned > scPoli zero-shot**: Fine-tuning should help adapt to new microscopes
2. **scVI online > Symphony**: Non-linear methods should capture complex microscope effects better
3. **Symphony > PCA-only**: Even linear batch correction improves mapping
4. **Per-microscope variation**: CV8000 query cells (shared microscope type with reference) should map better than Operetta (unseen in reference)
5. **MOA mAP lower than compound mAP**: MOA grouping is coarser, but preservation indicates functional signal maintained
6. **Cross-set mAP < within-set mAP**: Mapping is harder than within-batch correction, so absolute scores will be lower
