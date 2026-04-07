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
        script="scripts/correct_with_combat.py",
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_combat.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_combat.log"
    container:
        "containers/base.sif"  # combat only needs scanpy, which is in base
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_sphering.log"
    container:
        "containers/base.sif"
    params:
        method="ZCA-cor",
        column_norm="Metadata_JCP2022",
        values_norm="DMSO",
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

rule methods_harmony_v1:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_harmony.py",
        parameter_path="outputs/{scenario}/optimization/optuna_harmony_v1.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_harmony_v1.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_harmony_v1.log"
    container:
        "containers/harmony_v1.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode 'harmony' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_harmony_v2:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_harmony.py",
        parameter_path="outputs/{scenario}/optimization/optuna_harmony_v2.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_harmony_v2.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_harmony_v2.log"
    container:
        "containers/lightweight.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode 'harmony' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --parameter_path '{input.parameter_path}' \
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_harmony_pca.log"
    container:
        "containers/lightweight.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest="--smoketest" if config["smoketest"] else "",
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
        script="scripts/correct_with_scanorama.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scanorama.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanorama.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scanorama.log"
    container:
        "containers/lightweight.sif"
    params:
        method="scanorama",
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --mode '{params.method}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --parameter_path '{input.parameter_path}' \
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scanorama_pca.log"
    container:
        "containers/lightweight.sif"
    params:
        method="scanorama_pca",
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_mnn.log"
    conda:
        "../envs/mnn.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
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
        script="scripts/correct_with_desc.py",
        parameter_path="outputs/{scenario}/optimization/optuna_desc.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_desc.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_desc.log"
    container:
        "containers/desc.sif"
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
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scvi_single:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scvi.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scvi_single.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scvi_single.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scvi_single.log"
    container:
        "containers/scvi.sif"
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
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scvi_multi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scvi.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scvi_multi.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scvi_multi.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scvi_multi.log"
    container:
        "containers/scvi.sif"
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
            --multi \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scvi_normal:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scvi.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scvi_normal.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scvi_normal.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scvi_normal.log"
    container:
        "containers/scvi.sif"
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
            --parameter_path '{input.parameter_path}' \
            --gene_likelihood normal \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """

rule methods_scanvi_single:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        scvi_parameter_path="outputs/{scenario}/optimization/optuna_scvi_single.csv",
        scanvi_parameter_path="outputs/{scenario}/optimization/optuna_scanvi_single.csv",
        script="scripts/correct_with_scanvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanvi_single.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scanvi_single.log"
    container:
        "containers/scvi.sif"
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
            --scvi_parameter_path '{input.scvi_parameter_path}' \
            --scanvi_parameter_path '{input.scanvi_parameter_path}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """


rule methods_scanvi_multi:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        scvi_parameter_path="outputs/{scenario}/optimization/optuna_scvi_multi.csv",
        scanvi_parameter_path="outputs/{scenario}/optimization/optuna_scanvi_multi.csv",
        script="scripts/correct_with_scanvi.py"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scanvi_multi.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scanvi_multi.log"
    container:
        "containers/scvi.sif"
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
            --multi \
            --scvi_parameter_path '{input.scvi_parameter_path}' \
            --scanvi_parameter_path '{input.scanvi_parameter_path}' \
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_gaushvi.log"
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_gaushanvi.log"
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
        script="scripts/correct_with_sysvi.py",
        parameter_path="outputs/{scenario}/optimization/optuna_sysvi.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_sysvi.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_sysvi.log"
    container:
        "containers/scvi.sif"
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

rule methods_scpoli:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_scpoli.py",
        parameter_path="outputs/{scenario}/optimization/optuna_scpoli.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_scpoli.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scpoli.log"
    container:
        "containers/scpoli.sif"
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
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_scpoli_pca.log"
    container:
        "containers/scpoli.sif"
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
        script="scripts/correct_with_fastmnn.R",
        parameter_path="outputs/{scenario}/optimization/optuna_fastmnn.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_fastMNN.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_fastmnn.log"
    container:
        "containers/r_v4.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_cca_v4:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat_v4.R",
        parameter_path="outputs/{scenario}/optimization/optuna_seurat_cca_v4.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_cca_v4.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_seurat_cca_v4.log"
    container:
        "containers/r_v4.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="cca",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_rpca_v4:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat_v4.R",
        parameter_path="outputs/{scenario}/optimization/optuna_seurat_rpca_v4.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_rpca_v4.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_seurat_rpca_v4.log"
    container:
        "containers/r_v4.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="rpca",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_cca_v5:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat_v5.R",
        parameter_path="outputs/{scenario}/optimization/optuna_seurat_cca_v5.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_cca_v5.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_seurat_cca_v5.log"
    container:
        "containers/r_v5.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="cca",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            &> '{log}'
        """

rule methods_seurat_rpca_v5:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/correct_with_seurat_v5.R",
        parameter_path="outputs/{scenario}/optimization/optuna_seurat_rpca_v5.csv"
    output:
        path="outputs/{scenario}/" + config["preproc"] + "_seurat_rpca_v5.parquet"
    log:
        "outputs/{scenario}/logs/" + config["preproc"] + "_correct_seurat_rpca_v5.log"
    container:
        "containers/r_v5.sif"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        method="rpca",
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        Rscript '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --method '{params.method}' \
            --parameter_path '{input.parameter_path}' \
            --output_path '{output.path}' \
            &> '{log}'
        """
