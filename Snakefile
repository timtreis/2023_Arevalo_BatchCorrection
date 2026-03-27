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
    "sysvi",
    "scanvi_single",
    "scanvi_multi",
    # "gaushvi",  # deferred: custom scvi-tools fork needs investigation
    # "gaushanvi",  # deferred: custom scvi-tools fork needs investigation
    "scpoli",
    # "scpoli_pca", # performs quite a bit worse than normal scpoli
    "sphering",
    "seurat_cca",
    "seurat_rpca",
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

# NOTE: scANVI and scPoli use label coarsening (utils.coarsen_labels) to handle
# high-cardinality label columns.  Compounds appearing in <5 batches are marked
# "Unknown" during training.  No Snakefile-level skip is needed.


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
