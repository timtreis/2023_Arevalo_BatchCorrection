# Arevalo et al. (2024) — Published Baseline Results

**Paper**: "Evaluating batch correction methods for image-based cell profiling"
**Journal**: Nature Communications 15, 6516 (2024)
**DOI**: 10.1038/s41467-024-50613-5
**GitHub**: github.com/carpenter-singh-lab/2023_Arevalo_BatchCorrection

---

## Methods Evaluated (10 + baseline)

| Method | Type | Hyperparameter notes |
|--------|------|---------------------|
| Harmony | Linear embedding | clusters=300, iterations=20 (non-default) |
| Seurat CCA | Linear (R) | Default |
| Seurat RPCA | Linear (R) | Default |
| Scanorama | Linear embedding | KNN=20, alpha=0.1, sigma=15 (default) |
| fastMNN | Linear (R) | Default |
| MNN | Linear (Python/mnnpy) | neighbor_size=20 (default) |
| Combat | Linear | scanpy implementation, no hyperparams |
| DESC | Deep learning | Convergence hyperparams adjusted (defaults collapsed to -1/1 tanh) |
| scVI | Deep learning (VAE) | n_latent=30 (non-default, was 10), n_hidden=128, dropout=0.1; data shifted: x_i = x_i - min(x) + 1 |
| Sphering | Statistical normalization | ZCA-cor |
| Baseline | None (uncorrected) | — |

**Not in our pipeline**: MNN (replaced by fastMNN), BBKNN (excluded — modifies kNN graph, not profiles)
**New in our pipeline**: harmony_v1, harmony_v2, scvi_single, scvi_multi, scanvi_single, scanvi_multi, sysvi, scpoli

---

## Scoring Formula

```
Overall = 0.4 × mean(batch_metrics) + 0.6 × mean(bio_metrics)
```

### Batch correction metrics (4):
- `graph_conn` — Graph connectivity
- `kbet` — k-BET rejection rate
- `lisi_batch` — Batch LISI
- `silhouette_batch` — Batch silhouette (cosine)

### Bio conservation metrics (6):
- `lisi_label` — Label LISI
- `ari` — Leiden ARI
- `nmi` — Leiden NMI
- `asw` — Label silhouette
- `negcon_mean_map` — mAP for negative controls (DMSO)
- `nonrep_mean_map` — mAP non-replicate retrieval

### Excluded from scoring (computed but redlisted):
- `pcr`, `pcr_batch`, `il_f1`, `il_asw`, fraction-below-p variants

---

## Published Aggregate Scores (Supplementary Table 2, mean across all scenarios)

| Method | Batch | Bio | Overall | Rank |
|--------|-------|-----|---------|------|
| Seurat CCA | 0.54 | 0.48 | 0.50 | 1 |
| Seurat RPCA | 0.50 | 0.48 | 0.49 | 2 |
| Harmony | 0.49 | 0.47 | 0.48 | 3 |
| Scanorama | 0.55 | 0.46 | 0.49 | 4 |
| fastMNN | 0.51 | 0.46 | 0.48 | 5 |
| scVI | 0.45 | 0.46 | 0.46 | 6 |
| Combat | 0.40 | 0.46 | 0.44 | 7 |
| MNN | 0.41 | 0.45 | 0.44 | 8 |
| Sphering | 0.39 | 0.45 | 0.43 | 9 |
| Baseline | 0.40 | 0.45 | 0.43 | 10 |
| DESC | 0.37 | 0.43 | 0.40 | 11 |

---

## Per-Scenario Detailed Scores

### Scenario 1 — Single source (source_6), TARGET2, batch_key=Metadata_Batch

| Method | GC | kBET | LISI_b | Sil_b | LISI_l | ARI | NMI | Sil_l | mAP_c | mAP_n |
|--------|-----|------|--------|-------|--------|-----|-----|-------|-------|-------|
| fastMNN | 0.81 | 0.46 | 0.44 | 0.76 | 0.98 | 0.05 | 0.58 | 0.55 | 0.95 | 0.63 |
| Seurat RPCA | 0.81 | 0.44 | 0.43 | 0.78 | 0.98 | 0.06 | 0.60 | 0.53 | 0.95 | 0.60 |
| Seurat CCA | 0.80 | 0.42 | 0.46 | 0.77 | 0.98 | 0.07 | 0.59 | 0.52 | 0.95 | 0.61 |
| Harmony | 0.79 | 0.42 | 0.40 | 0.74 | 0.98 | 0.07 | 0.60 | 0.55 | 0.95 | 0.60 |
| Combat | 0.78 | 0.44 | 0.41 | 0.84 | 0.98 | 0.05 | 0.56 | 0.53 | 0.94 | 0.55 |
| MNN | 0.78 | 0.44 | 0.40 | 0.84 | 0.98 | 0.07 | 0.56 | 0.52 | 0.94 | 0.54 |
| DESC | 0.78 | 0.39 | 0.42 | 0.74 | 0.98 | 0.06 | 0.58 | 0.55 | 0.95 | 0.56 |
| Baseline | 0.77 | 0.41 | 0.39 | 0.85 | 0.98 | 0.06 | 0.55 | 0.52 | 0.94 | 0.52 |
| Sphering | 0.74 | 0.39 | 0.37 | 0.87 | 0.98 | 0.06 | 0.54 | 0.52 | 0.94 | 0.51 |
| Scanorama | 0.54 | 0.72 | 0.53 | 0.59 | 0.98 | 0.09 | 0.56 | 0.45 | 0.94 | 0.49 |
| scVI | 0.72 | 0.39 | 0.43 | 0.70 | 0.98 | 0.06 | 0.58 | 0.53 | 0.94 | 0.56 |

### Scenario 2 — 3 sources, TARGET2, batch_key=Metadata_Source

| Method | GC | kBET | LISI_b | Sil_b | LISI_l | ARI | NMI | Sil_l | mAP_c | mAP_n |
|--------|-----|------|--------|-------|--------|-----|-----|-------|-------|-------|
| Scanorama | 0.53 | 0.66 | 0.67 | 0.79 | 0.98 | 0.04 | 0.42 | 0.47 | 0.92 | 0.46 |
| Seurat CCA | 0.64 | 0.28 | 0.55 | 0.86 | 0.98 | 0.04 | 0.45 | 0.49 | 0.93 | 0.47 |
| Harmony | 0.64 | 0.26 | 0.55 | 0.86 | 0.98 | 0.04 | 0.42 | 0.48 | 0.92 | 0.46 |
| Seurat RPCA | 0.65 | 0.24 | 0.46 | 0.84 | 0.98 | 0.05 | 0.46 | 0.49 | 0.93 | 0.47 |
| fastMNN | 0.56 | 0.27 | 0.52 | 0.80 | 0.97 | 0.04 | 0.39 | 0.45 | 0.92 | 0.45 |
| scVI | 0.59 | 0.16 | 0.48 | 0.74 | 0.98 | 0.04 | 0.45 | 0.47 | 0.90 | 0.47 |
| MNN | 0.62 | 0.10 | 0.19 | 0.76 | 0.98 | 0.04 | 0.42 | 0.48 | 0.90 | 0.42 |
| Combat | 0.62 | 0.05 | 0.16 | 0.76 | 0.98 | 0.03 | 0.42 | 0.48 | 0.91 | 0.42 |
| DESC | 0.59 | 0.11 | 0.32 | 0.67 | 0.97 | 0.03 | 0.38 | 0.42 | 0.89 | 0.42 |
| Baseline | 0.60 | 0.04 | 0.17 | 0.75 | 0.98 | 0.04 | 0.42 | 0.48 | 0.90 | 0.41 |
| Sphering | 0.58 | 0.03 | 0.13 | 0.76 | 0.98 | 0.03 | 0.40 | 0.48 | 0.91 | 0.39 |

### Scenario 3 — 3 sources, TARGET2+COMPOUND, batch_key=Metadata_Source

| Method | GC | kBET | LISI_b | Sil_b | LISI_l | ARI | NMI | Sil_l | mAP_c | mAP_n |
|--------|-----|------|--------|-------|--------|-----|-----|-------|-------|-------|
| Seurat CCA | 0.55 | 0.17 | 0.51 | 0.79 | 1.00 | 0.02 | 0.30 | 0.29 | 0.28 | 0.06 |
| Scanorama | 0.54 | 0.33 | 0.36 | 0.76 | 1.00 | 0.01 | 0.29 | 0.27 | 0.20 | 0.05 |
| fastMNN | 0.54 | 0.14 | 0.48 | 0.72 | 1.00 | 0.02 | 0.30 | 0.21 | 0.24 | 0.06 |
| Seurat RPCA | 0.55 | 0.15 | 0.28 | 0.79 | 1.00 | 0.02 | 0.32 | 0.30 | 0.25 | 0.05 |
| Harmony | 0.54 | 0.06 | 0.25 | 0.79 | 1.00 | 0.02 | 0.32 | 0.31 | 0.24 | 0.05 |
| scVI | 0.54 | 0.07 | 0.27 | 0.74 | 1.00 | 0.02 | 0.32 | 0.23 | 0.24 | 0.05 |
| Sphering | 0.54 | 0.02 | 0.00 | 0.80 | 1.00 | 0.00 | 0.26 | 0.35 | 0.21 | 0.03 |
| Combat | 0.54 | 0.03 | 0.02 | 0.76 | 1.00 | 0.01 | 0.29 | 0.30 | 0.22 | 0.04 |
| MNN | 0.54 | 0.03 | 0.01 | 0.75 | 1.00 | 0.01 | 0.29 | 0.29 | 0.20 | 0.04 |
| Baseline | 0.54 | 0.03 | 0.00 | 0.74 | 1.00 | 0.00 | 0.28 | 0.29 | 0.19 | 0.03 |
| DESC | 0.53 | 0.00 | 0.15 | 0.43 | 1.00 | 0.02 | 0.38 | 0.04 | 0.14 | 0.02 |

### Scenario 4 — 5 sources, TARGET2, batch_key=Metadata_Source

Not individually extractable from the HTML supplementary (shown as Figure 5 in the main paper). Intermediate difficulty between Scenarios 2 and 5.

### Scenario 5 — 5 sources, TARGET2+COMPOUND, batch_key=Metadata_Source

| Method | GC | kBET | LISI_b | Sil_b | LISI_l | ARI | NMI | Sil_l | mAP_c | mAP_n |
|--------|-----|------|--------|-------|--------|-----|-----|-------|-------|-------|
| Seurat CCA | 0.35 | 0.13 | 0.49 | 0.85 | 1.00 | 0.01 | 0.25 | 0.31 | 0.27 | 0.07 |
| Scanorama | 0.34 | 0.14 | 0.42 | 0.83 | 1.00 | 0.01 | 0.29 | 0.27 | 0.23 | 0.05 |
| Seurat RPCA | 0.35 | 0.10 | 0.31 | 0.84 | 1.00 | 0.02 | 0.28 | 0.31 | 0.27 | 0.06 |
| fastMNN | 0.34 | 0.08 | 0.45 | 0.79 | 1.00 | 0.01 | 0.25 | 0.22 | 0.26 | 0.07 |
| Harmony | 0.34 | 0.03 | 0.25 | 0.84 | 1.00 | 0.02 | 0.28 | 0.32 | 0.26 | 0.05 |
| scVI | 0.34 | 0.04 | 0.23 | 0.77 | 1.00 | 0.01 | 0.28 | 0.25 | 0.24 | 0.05 |
| Sphering | 0.34 | 0.01 | 0.00 | 0.84 | 1.00 | 0.00 | 0.21 | 0.38 | 0.22 | 0.03 |
| Combat | 0.34 | 0.01 | 0.02 | 0.81 | 1.00 | 0.00 | 0.24 | 0.31 | 0.24 | 0.04 |
| MNN | 0.34 | 0.03 | 0.01 | 0.81 | 1.00 | 0.00 | 0.23 | 0.31 | 0.22 | 0.04 |
| Baseline | 0.34 | 0.01 | 0.01 | 0.80 | 1.00 | 0.00 | 0.23 | 0.31 | 0.21 | 0.03 |
| DESC | 0.33 | 0.01 | 0.27 | 0.58 | 1.00 | 0.02 | 0.34 | 0.01 | 0.22 | 0.04 |

---

## Key Observations for Reproduction

1. **Rankings are tight**: Overall scores range 0.40–0.50. Small differences matter.
2. **Scenario 1 is easiest**: Most methods score well (batch metrics 0.37–0.72, bio ~0.50–0.63).
3. **Scenarios 3 & 5 are hardest**: Adding COMPOUND plates dramatically drops bio scores (mAP drops to 0.14–0.28).
4. **DESC is consistently worst**: Often collapses, especially on multi-source scenarios.
5. **Seurat CCA ranks #1 overall**: Driven by strong batch correction across scenarios.
6. **scVI ranks #6**: Despite being deep learning — defaults were suboptimal (paper notes n_latent=30 was non-default).

## Reproduction Criteria

From the research plan: correlation with published results must be **>0.95** for the original 11 methods. Since we don't have MNN in our pipeline, compare the 9 overlapping methods:
- Harmony, Seurat CCA, Seurat RPCA, Scanorama, fastMNN, Combat, DESC, scVI, Sphering (+ Baseline)

## Potential Discrepancies to Watch

1. **scVI hyperparams**: Arevalo used n_latent=30. Our pipeline tunes this via Optuna — default run should match their settings.
2. **Harmony**: Arevalo used clusters=300, iterations=20. Check our defaults.
3. **DESC**: Known to collapse. Need to verify our convergence hyperparams match.
4. **Metric versions**: scib-metrics API has changed (bras vs silhouette_batch). Our compatibility shim should handle this.
5. **Data preprocessing**: Must use identical MAD → INT → featselect pipeline with same parameters.
