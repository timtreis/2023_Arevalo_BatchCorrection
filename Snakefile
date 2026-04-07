# Define processing workflows and correction methods to run
WORKFLOWS = [
    "mad_int_featselect",
]
METHODS = [
    "scanorama",
    # "scanorama_pca",
    # "mnn", # skipped: ancient deps, replaced by fastMNN
    "fastMNN",
    "harmony_v1",
    "harmony_v2",
    # "harmony_pca", # performs the same as normal harmony
    "combat",
    "desc",
    "scvi_single",
    "scvi_multi",
    # "scvi_normal",  # dropped: scVI's gene_likelihood="normal" has unstable exp() variance
    #                  # parameterization (NaN crashes, GitHub #2103/#516). sysVI already
    #                  # provides a proper Gaussian VAE with softplus + clamping + nan_to_num.
    "sysvi",
    "scanvi_single",
    "scanvi_multi",
    # "gaushvi",  # superseded by scvi_normal (which is itself superseded by sysvi)
    # "gaushanvi",  # superseded: scANVI + Gaussian likelihood, same stability issues
    "scpoli",
    # "scpoli_pca", # performs quite a bit worse than normal scpoli
    "sphering",
    "seurat_cca_v4",
    "seurat_rpca_v4",
    "seurat_cca_v5",
    "seurat_rpca_v5",
    # "cpDistiller_B",
    # "cpDistiller_S",
    # "cpDistiller_SBP",
]

# DESC has a hardcoded 200K cell limit that causes crashes on larger datasets.
# Skip it automatically when the preprocessed input exceeds that threshold.
DESC_MAX_CELLS = 200_000
_preproc_path = f"outputs/{config.get('scenario', '')}/{config.get('preproc', 'mad_int_featselect')}.parquet"
try:
    import pyarrow.parquet as pq
    _n_cells = pq.read_metadata(_preproc_path).num_rows
    if _n_cells > DESC_MAX_CELLS:
        METHODS = [m for m in METHODS if m != "desc"]
        print(f"NOTE: Skipping DESC (input has {_n_cells:,} cells, limit is {DESC_MAX_CELLS:,})")
except (FileNotFoundError, Exception):
    pass  # file doesn't exist yet (preprocessing not done); keep desc in the list

# R-based methods (fastMNN, Seurat CCA/RPCA) fail on large datasets due to R's
# 2^31-1 element vector limit in intermediate computations (e.g. batchelor's
# .average_correction smoothing kernel) and memory allocation patterns in Seurat's
# IntegrateLayers. See tasks/r_methods_scalability.md for full analysis.
R_METHODS_MAX_CELLS = 100_000
_r_methods = ("fastMNN", "seurat_cca_v4", "seurat_rpca_v4", "seurat_cca_v5", "seurat_rpca_v5")
try:
    _n_cells_r = _n_cells  # reuse from DESC check above
except NameError:
    try:
        import pyarrow.parquet as pq
        _n_cells_r = pq.read_metadata(_preproc_path).num_rows
    except (FileNotFoundError, Exception):
        _n_cells_r = 0
if _n_cells_r > R_METHODS_MAX_CELLS:
    METHODS = [m for m in METHODS if m not in _r_methods]
    print(f"NOTE: Skipping R methods {list(_r_methods)} (input has {_n_cells_r:,} cells, limit is {R_METHODS_MAX_CELLS:,})")

# scANVI's loss function uses broadcast_labels which is O(batch × n_labels × latent).
# With COMPOUND plates, even after coarsening, the label count can be thousands,
# causing OOM in broadcast_labels (e.g. 555 GiB allocation on scenario_3).
# Skip scANVI for scenarios that include COMPOUND plates.
if "COMPOUND" in config.get("plate_types", []):
    METHODS = [m for m in METHODS if m not in ("scanvi_single", "scanvi_multi")]
    print("NOTE: Skipping scANVI (COMPOUND plates → too many labels for broadcast_labels)")

# Skip methods whose HPO produced zero COMPLETE trials (e.g. Seurat CCA failing
# on certain data configurations). Without valid params, correction will crash.
import csv as _csv
_scenario = config.get("scenario", "")
_methods_to_skip = []
for _m in list(METHODS):
    _optuna_name = _m.lower().replace("mnn", "mnn")  # fastMNN -> optuna_fastmnn
    _optuna_path = f"outputs/{_scenario}/optimization/optuna_{_optuna_name}.csv"
    try:
        with open(_optuna_path) as _f:
            _reader = _csv.DictReader(_f)
            if not any(r.get("state") == "COMPLETE" for r in _reader):
                _methods_to_skip.append(_m)
    except (FileNotFoundError, Exception):
        pass  # HPO not run yet; keep in list
if _methods_to_skip:
    METHODS = [m for m in METHODS if m not in _methods_to_skip]
    print(f"NOTE: Skipping {_methods_to_skip} (0 COMPLETE HPO trials)")

# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"
include: "rules/plots.smk"
include: "rules/tune.smk"


rule all:
    input:
        expand(plots_pattern, plot=PDF_PLOTS, ext=["pdf"]),
