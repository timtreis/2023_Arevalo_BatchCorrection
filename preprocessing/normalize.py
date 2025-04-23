import numpy as np
import pandas as pd
from preprocessing.io import merge_parquet, split_parquet
import anndata as ad
import statsmodels.formula.api as smf
from tqdm.auto import tqdm


def mad(variant_feats_path, neg_stats_path, normalized_path):
    meta, vals, features = split_parquet(variant_feats_path)
    neg_stats = pd.read_parquet(neg_stats_path)
    neg_stats = neg_stats.query('feature in @features')

    # get counts and sort by plate
    plates, counts = np.unique(meta['Metadata_Plate'], return_counts=True)
    ix = np.argsort(meta['Metadata_Plate'])
    meta = meta.iloc[ix]
    vals = vals[ix]

    # get mad and median matrices for MAD normalization
    mads = neg_stats.pivot(index='Metadata_Plate',
                           columns='feature',
                           values='mad')
    mads = mads.loc[plates, features].values
    medians = neg_stats.pivot(index='Metadata_Plate',
                              columns='feature',
                              values='median')
    medians = medians.loc[plates, features].values

    # Get normalized features (epsilon = 0) for all plates that have MAD stats
    # -= and /= are inplace operations. i.e save memory
    vals -= np.repeat(medians, counts, axis=0)
    vals /= np.repeat(mads, counts, axis=0)

    merge_parquet(meta, vals, features, normalized_path)


def normalize_plate_effect(
    meta: pd.DataFrame,
    vals: np.ndarray,
    plate_col: str = "Metadata_Plate",
    drug_id_col: str = "Metadata_JCP2022",
    drug_id_for_replicates: str = "DMSO",
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Normalize out plate row/column effects by regressing on control replicates, separately per plate.

    Parameters
    ----------
    meta: pd.DataFrame
        Must contain columns:
          • plate_col: plate identifier
          • "row": plate row identifier
          • "column": plate column identifier
          • "well_type": with value "replicate" for control wells
    vals: np.ndarray[n_wells, n_features]
        Raw measurements.
    plate_col: str
        Name of the column in `meta` that identifies each plate.

    Returns
    -------
    meta_corrected: pd.DataFrame
        Same as input `meta`.
    vals_corrected: np.ndarray[n_wells, n_features]
        Bias‐corrected measurements.
    """
    # Wrap in AnnData
    adata = ad.AnnData(X=vals.copy(), obs=meta.copy())

    # Ensure categorical dtype
    adata.obs[plate_col] = adata.obs[plate_col].astype("category")
    adata.obs["row"]    = adata.obs["row"].astype("category")
    adata.obs["column"] = adata.obs["column"].astype("category")

    # Which wells are replicates?
    is_rep = adata.obs[drug_id_col] == drug_id_for_replicates

    # Prepare output array
    corrected = adata.X.copy()
    n_features = corrected.shape[1]

    # Loop over plates
    for plate in tqdm(adata.obs[plate_col].cat.categories, desc="Plates"):
        mask_plate = adata.obs[plate_col] == plate
        df_plate   = adata.obs.loc[mask_plate].copy()
        vals_plate = adata.X[mask_plate, :]

        # Local replicate mask
        rep_local = is_rep[mask_plate]
        if not rep_local.any():
            continue  # skip plates with no replicates

        # Regress out row+column for each feature
        for feat_idx in tqdm(range(n_features), desc=f"Features [{plate}]", leave=False):
            df_plate["value"] = vals_plate[:, feat_idx]
            model = smf.ols("value ~ C(row) + C(column)", data=df_plate[rep_local]).fit()
            bias  = model.predict(df_plate)
            corrected[mask_plate, feat_idx] = vals_plate[:, feat_idx] - bias

    return adata.obs, corrected