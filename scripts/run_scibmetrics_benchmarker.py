import sys
import logging
from scvi.external import SysVI
from preprocessing import io

logger = logging.getLogger(__name__)

def run_scibmetrics_benchmarker(adata_path, output_path, eval_keys):
    print(eval_keys)

if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
    eval_keys = sys.argv[3]
    run_scibmetrics_benchmarker(adata_path, output_path, eval_keys)
