
rule optimize_scpoli:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scpoli.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scpoli.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_scpoli.log"
    conda:
        "../envs/scpoli.yaml"
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
    conda:
        "../envs/scvi.yaml"
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
    conda:
        "../envs/scvi.yaml"
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
    conda:
        "../envs/scvi.yaml"
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
    conda:
        "../envs/scvi.yaml"
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
    conda:
        "../envs/sysvi.yaml"
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

rule optimize_harmony:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_harmony.py"
    output:
        path="outputs/{scenario}/optimization/optuna_harmony.csv"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_optimize_harmony.log"
    conda:
        "../envs/harmony.yaml"
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
