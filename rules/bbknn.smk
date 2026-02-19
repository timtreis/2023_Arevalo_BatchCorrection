BBKNN_METRICS = ["graph_conn", "kbet", "lisi_label", "lisi_batch", "ari", "nmi"]

rule bbknn_clustering:
    input:
        data="outputs/{prefix}/{pipeline}.parquet",
        script="scripts/correct_with_bbknn.py",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_bbknn_clusters.h5ad",
    log:
        "logs/{prefix}/{criteria}_{pipeline}_bbknn_clustering.log",
    params:
        batch_key=config["batch_key"],
    conda:
        "../envs/bbknn.yaml"
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --output_path '{output}' \
            &> '{log}'
        """


rule bbknn_all:
    input:
        expand(
            "outputs/{{scenario}}/metrics/{{criteria}}/scib/{{pipeline}}_{metric}.npy",
            metric=BBKNN_METRICS,
        ),
    output:
        output_path="outputs/{scenario}/metrics/{criteria}/{pipeline}_bbknn_scib.parquet",
    log:
        "logs/{scenario}/{criteria}_{pipeline}_bbknn_all.log",
    run:
        metrics.scib.concat(*input, **output)


ruleorder: bbknn_clustering > clustering
ruleorder: bbknn_all > scib_all
