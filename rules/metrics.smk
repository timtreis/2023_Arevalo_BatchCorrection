include: "scib.smk"
include: "consistency.smk"
include: "bbknn.smk"

import logging


metrics_baseline_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_all_metrics.parquet"
)
metrics_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_all_metrics.parquet"
)


rule metrics_concat_and_average_all_metrics:
    input:
        scib_path=f"outputs/{scenario}/metrics/{config['preproc']}_scibmetrics_benchmarker.parquet",
        map_path=f"outputs/{scenario}/metrics/{config['preproc']}_map.parquet",
    output:
        path=f"outputs/{scenario}/metrics/{config['preproc']}_all_metrics.parquet",
    params:
        metrics_redlist=[
            # "pcr",
            # "pcr_batch",
            # "il_f1",
            # "il_asw",
            "nonrep_frac_p",
            "nonrep_frac_q",
            "negcon_frac_p",
            "negcon_frac_q",
        ],
        methods_redlist=[],
    run:
        metrics.concat_and_average_all_metrics(
            input.scib_path,
            input.map_path,
            params.metrics_redlist,
            params.methods_redlist,
            output.path,
        )

rule metrics_run_scibmetrics_benchmarker:
    input:
        script="scripts/run_scibmetrics_benchmarker.py",
        adata_path=f"outputs/{scenario}/{config['preproc']}_all_methods.h5ad",
    output:
        path=f"outputs/{scenario}/metrics/{config['preproc']}_scibmetrics_benchmarker.parquet",
    params:
        eval_keys=config["eval_key"],
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        methods=METHODS
    log:
        f"outputs/{scenario}/logs/{config['preproc']}_scibmetrics_benchmarker.log"
    container:
        "containers/scibmetrics-gpu.sif"
    resources:
        nvidia_gpu=1
    threads: 6
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) && \
        case "$GPU_NAME" in *H100*|*H200*|*B100*|*B200*) export CUDA_VISIBLE_DEVICES=''; echo "GPU $GPU_NAME: FAISS not compiled for sm_90+, forcing CPU" ;; esac && \
        python {input.script} '{input.adata_path}' '{output.path}' '{params.batch_key}' '{params.eval_keys}' '{params.methods}' &> {log}
        """
        
rule metrics_mean_average_precision:
    input:
        adata=f"outputs/{scenario}/{config['preproc']}_all_methods.h5ad",
    output:
        path=f"outputs/{scenario}/metrics/{config['preproc']}_map.parquet",
    params:
        plate_type=config["plate_types"],
        eval_keys=config["eval_key"],
    log:
        path=f"outputs/{scenario}/logs/{config['preproc']}_map_metrics.log"
    run:
        logging.basicConfig(
            filename=log.path,
            filemode="w",
            format="%(asctime)s - %(levelname)s - %(message)s",
            level=logging.INFO
        )
        
        logger = logging.getLogger()
        metrics.mean_average_precision(input.adata, output.path, params.plate_type, params.eval_keys)
