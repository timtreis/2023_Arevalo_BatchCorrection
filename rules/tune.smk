
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
        nvidia_gpu=1
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
        nvidia_gpu=1
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
        nvidia_gpu=1
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
        nvidia_gpu=1
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
        nvidia_gpu=1
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
        nvidia_gpu=1
    shell:
        """
        export JAX_PLATFORMS=cpu && \
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
        nvidia_gpu=1
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
        nvidia_gpu=1
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
    shell:
        """
        export JAX_PLATFORMS=cpu && \
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
        nvidia_gpu=1
    shell:
        """
        export JAX_PLATFORMS=cpu && \
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
        script="scripts/optimise_fastmnn.R"
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
    shell:
        """
        Rscript '{input.script}' \
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
        script="scripts/optimise_seurat.R"
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
    shell:
        """
        Rscript '{input.script}' \
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
        script="scripts/optimise_seurat.R"
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
    shell:
        """
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --method 'rpca' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """
