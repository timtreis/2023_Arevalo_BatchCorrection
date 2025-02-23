rule aggregate_method_outputs_into_adata:
    input:
        unintegrated="outputs/{scenario}/" + config["preproc"] + ".parquet",
        integrated=expand(
            "outputs/{{scenario}}/" + config["preproc"] + "_{method}.parquet",
            method=METHODS,
        )
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_all_methods.h5ad"
    run:
        metrics.scib.aggregate_method_outputs_into_adata(input.unintegrated, input.integrated, output.path)

rule methods_combat:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_combat.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_combat.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_combat.log"
    conda:
        "../envs/harmony.yaml"  # we only need scanpy so this will do
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_sphering:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_sphering.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_sphering.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_sphering.log"
    conda:
        "../envs/sphering.yaml"
    params:
        method="ZCA-cor",
        column_norm="Metadata_JCP2022",
        values_norm="DMSO",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --method '{params.method}' \
            --column_norm '{params.column_norm}' \
            --values_norm '{params.values_norm}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_harmony:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_harmony.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_harmony.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_harmony.log"
    conda:
        "../envs/harmony.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode 'harmony' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_harmony_pca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_harmony.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_harmony_pca.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_harmony_pca.log"
    conda:
        "../envs/harmony.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode 'harmony_pca' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scanorama:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scanorama.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanorama.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scanorama.log"
    conda:
        "../envs/scanorama.yaml"
    params:
        method="scanorama",
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode '{params.method}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_scanorama_pca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scanorama.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanorama_pca.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scanorama_pca.log"
    conda:
        "../envs/scanorama.yaml"
    params:
        method="scanorama_pca",
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode '{params.method}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_mnn:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_mnn.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_mnn.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_mnn.log"
    conda:
        "../envs/mnn.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_desc:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_desc.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_desc.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_desc.log"
    conda:
        "../envs/desc.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scvi.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scvi.log"
    conda:
        "../envs/scvi.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scanvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scanvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanvi.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scanvi.log"
    conda:
        "../envs/scvi.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_gaushvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_gaushvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_gaushvi.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_gaushvi.log"
    conda:
        "../envs/gaushvi.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_gaushanvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_gaushanvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_gaushanvi.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_gaushanvi.log"
    conda:
        "../envs/gaushvi.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_sysvi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_sysvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_sysvi.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_sysvi.log"
    conda:
        "../envs/sysvi.yaml"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scpoli:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scpoli.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scpoli.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scpoli.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scpoli.log"
    conda:
        "../envs/scpoli.yaml"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
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
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule methods_scpoli_pca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scpoli.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scpoli_pca.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_scpoli_pca.log"
    conda:
        "../envs/scpoli.yaml"
    params:
        batch_key=','.join(config["batch_key"]),
        label_key=config["label_key"],
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
            --output_path '{output.path}' \
            --preproc 'pca' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_fastMNN:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_fastmnn.R"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_fastMNN.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_fastmnn.log"
    conda:
        "../envs/fastmnn.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_cca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat.R"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_cca.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_seurat_cca.log"
    conda:
        "../envs/seurat.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="cca",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_rpca:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat.R"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_rpca.parquet"
    log:
        "logs/{scenario}/" + config["preproc"] + "_correct_seurat_rpca.log"
    conda:
        "../envs/seurat.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="rpca",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --output_path '{output.path}' \
            &> '{log}'
        """
