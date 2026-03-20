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
