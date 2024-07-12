# scripts/harmony.py

import sys
import logging
from harmonypy import run_harmony
from preprocessing import io
import scanpy as sc

logger = logging.getLogger(__name__)

def harmony(dframe_path, batch_key, output_path):
    '''Harmony correction'''
    meta, feats, features = io.split_parquet(dframe_path)
    harmony_out = run_harmony(feats,
                              meta,
                              batch_key,
                              max_iter_harmony=20,
                              nclust=300)  # Number of compounds

    feats = harmony_out.Z_corr.T
    features = [f'harmony_{i}' for i in range(feats.shape[1])]
    io.merge_parquet(meta, feats, features, output_path)

def pca_harmony(dframe_path, batch_key, output_path):
    '''Harmony correction with PCA'''
    meta, feats, features = io.split_parquet(dframe_path)
    n_latent = min(feats.shape) - 1  # required for arpack
    logger.info('Computing PCA...')
    sc.tl.pca(feats, n_comps=n_latent)  # Generates X_pca
    logger.info('Computing PCA Done.')
    harmony_out = run_harmony(feats,
                              meta,
                              batch_key,
                              max_iter_harmony=20,
                              nclust=300)  # Number of compounds

    feats = harmony_out.Z_corr.T
    features = [f'harmony_{i}' for i in range(feats.shape[1])]
    io.merge_parquet(meta, feats, features, output_path)

if __name__ == "__main__":
    mode = sys.argv[1]
    dframe_path = sys.argv[2]
    batch_key = sys.argv[3]
    output_path = sys.argv[4]

    if mode == "harmony":
        harmony(dframe_path, batch_key, output_path)
    elif mode == "pca_harmony":
        pca_harmony(dframe_path, batch_key, output_path)
    else:
        raise ValueError("Invalid mode. Choose either 'harmony' or 'pca_harmony'")
