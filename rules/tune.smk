
# --- RAM estimation for scheduling ---
# Compute the base RAM footprint from the input parquet dimensions.
# Each method gets a multiplier reflecting its peak memory pattern.
# Snakemake uses mem_mb resources to avoid oversubscribing node RAM.
import os as _os

def _estimate_base_mem_mb():
    """Estimate base RAM (MB) from input data: n_cells × n_features × 8 bytes."""
    _path = f"outputs/{config.get('scenario', '')}/{config.get('preproc', 'mad_int_featselect')}.parquet"
    try:
        import pyarrow.parquet as _pq
        _meta = _pq.read_metadata(_path)
        _schema = _pq.read_schema(_path)
        n_cells = _meta.num_rows
        n_features = sum(1 for f in _schema.names if not f.startswith("Metadata_"))
        return int(n_cells * n_features * 8 / 1e6)  # raw matrix in MB
    except Exception:
        return 3000  # conservative 3 GB default

_BASE_MEM_MB = _estimate_base_mem_mb()

def _get_available_mem_mb():
    """Get available RAM in MB, respecting SLURM allocation if present."""
    # SLURM: use allocated memory (most reliable on shared nodes)
    mem_per_cpu = _os.environ.get("SLURM_MEM_PER_CPU")
    cpus = _os.environ.get("SLURM_CPUS_ON_NODE")
    mem_per_node = _os.environ.get("SLURM_MEM_PER_NODE")
    if mem_per_cpu and cpus:
        return int(mem_per_cpu) * int(cpus)
    if mem_per_node:
        return int(mem_per_node)
    # Fallback: use system available memory
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) // 1024  # kB → MB
    except Exception:
        pass
    return 256_000  # conservative 256 GB default

_AVAILABLE_MEM_MB = _get_available_mem_mb()

# RAM multipliers per method (peak RAM / raw matrix size):
#   R methods (Seurat/fastMNN): counts + data + scaled + PCA + integration copies → ~12×
#   GPU methods (scVI/scANVI/sysVI/scPoli/DESC): AnnData + model overhead → ~6×
#   Lightweight CPU (harmony/scanorama): data + PCA + working copies → ~5×
_MEM_MULTIPLIERS = {
    "seurat_cca": 12, "seurat_rpca": 12, "fastmnn": 12,
    "scvi_single": 6, "scvi_multi": 6, "scanvi_single": 6, "scanvi_multi": 6,
    "sysvi": 6, "scpoli": 6, "desc": 6,
    "harmony_v1": 5, "harmony_v2": 5, "scanorama": 5,
}

def _method_mem_mb(method):
    return int(_BASE_MEM_MB * _MEM_MULTIPLIERS.get(method, 6) * 1.2)  # 20% safety buffer

# When use_defaults is set, skip HPO and use pre-made default parameter files.
# This rule uses a {method} wildcard that matches all optimization outputs,
# so snakemake will prefer it over the method-specific rules below.
if config.get("use_defaults", False):
    ruleorder: use_default_params > optimize_scpoli > optimize_scvi_single > optimize_scvi_multi > optimize_scanvi_single > optimize_scanvi_multi > optimize_sysvi > optimize_harmony_v1 > optimize_harmony_v2 > optimize_scanorama > optimize_desc > optimize_fastmnn > optimize_seurat_cca > optimize_seurat_rpca

    rule use_default_params:
        input:
            "inputs/defaults/optuna_{method}.csv"
        output:
            "outputs/{scenario}/optimization/optuna_{method}.csv"
        shell:
            "cp '{input}' '{output}'"


rule optimize_scpoli:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scpoli.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scpoli.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scpoli.log"
    container:
        "containers/scpoli.sif"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("scpoli")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule optimize_scvi_single:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scvi.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scvi_single.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scvi_single.log"
    container:
        "containers/scvi.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("scvi_single")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule optimize_scvi_multi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scvi.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scvi_multi.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scvi_multi.log"
    container:
        "containers/scvi.sif"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("scvi_multi")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --multi \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule optimize_scanvi_single:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        scvi_params_path="outputs/{scenario}/optimization/optuna_scvi_single.csv",
        script="scripts/optimise_scanvi.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scanvi_single.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scanvi_single.log"
    container:
        "containers/scvi.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("scanvi_single")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --scvi_params_path '{input.scvi_params_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_scanvi_multi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        scvi_params_path="outputs/{scenario}/optimization/optuna_scvi_multi.csv",
        script="scripts/optimise_scanvi.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scanvi_multi.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scanvi_multi.log"
    container:
        "containers/scvi.sif"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("scanvi_multi")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --scvi_params_path '{input.scvi_params_path}' \
            --multi \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule optimize_sysvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_sysvi.py"
    output:
        path="outputs/{scenario}/optimization/optuna_sysvi.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_sysvi.log"
    container:
        "containers/scvi.sif"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("sysvi")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_harmony_v1:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_harmony.py"
    output:
        path="outputs/{scenario}/optimization/optuna_harmony_v1.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_harmony_v1.log"
    container:
        "containers/harmony_v1.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("harmony_v1")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_harmony_v2:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_harmony.py"
    output:
        path="outputs/{scenario}/optimization/optuna_harmony_v2.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_harmony_v2.log"
    container:
        "containers/lightweight.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("harmony_v2")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_scanorama:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scanorama.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scanorama.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scanorama.log"
    container:
        "containers/lightweight.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("scanorama")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_desc:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_desc.py"
    output:
        path="outputs/{scenario}/optimization/optuna_desc.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_desc.log"
    container:
        "containers/desc.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1,
        mem_mb=_method_mem_mb("desc")
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_fastmnn:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_fastmnn.py",
        trial_runner="scripts/run_fastmnn_trial.R"
    output:
        path="outputs/{scenario}/optimization/optuna_fastmnn.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_fastmnn.log"
    container:
        "containers/r.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("fastmnn")
    shell:
        """
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_seurat_cca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_seurat.py",
        trial_runner="scripts/run_seurat_trial.R"
    output:
        path="outputs/{scenario}/optimization/optuna_seurat_cca.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_seurat_cca.log"
    container:
        "containers/r.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("seurat_cca")
    shell:
        """
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --method 'cca' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule optimize_seurat_rpca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_seurat.py",
        trial_runner="scripts/run_seurat_trial.R"
    output:
        path="outputs/{scenario}/optimization/optuna_seurat_rpca.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_seurat_rpca.log"
    container:
        "containers/r.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        mem_mb=_method_mem_mb("seurat_rpca")
    shell:
        """
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --method 'rpca' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """
