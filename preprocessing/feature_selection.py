import logging

import pandas as pd
from pycytominer.operations import correlation_threshold, variance_threshold
import numpy as np
from scipy.interpolate import SmoothBivariateSpline


from .metadata import find_feat_cols

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def normalize_plate_effect_no_dask(
    meta: pd.DataFrame,
    vals: np.ndarray,
    plate_col: str = "Metadata_Plate",
    drug_id_col: str = "Metadata_JCP2022",
    drug_id_for_replicates: str = "DMSO",
    spline_s: float = 1.0,
    exclude_extremes: int = 0
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Sequential plate‐effect normalizer using 2D spline smoothing,
    dropping the n most extreme replicates per feature before fitting.

    Returns
    -------
    meta_cp : pd.DataFrame
        Copy of input `meta`, unchanged except indexed the same.
    corrected : np.ndarray
        Bias‐corrected matrix, same shape as `vals`.
    """
    meta_cp = meta.copy()
    corrected = np.empty_like(vals)
    
    plates = meta_cp[plate_col].astype(str).unique().tolist()
    for plate in plates:
        idx_plate = np.where(meta_cp[plate_col] == plate)[0]
        Vp = vals[idx_plate, :]
        rep_mask = (meta_cp.iloc[idx_plate][drug_id_col] == drug_id_for_replicates).values

        # skip plates with too few replicates
        if rep_mask.sum() < max(3, 2*exclude_extremes+1):
            corrected[idx_plate, :] = Vp
            continue

        sub_meta   = meta_cp.iloc[idx_plate]
        rows_all   = sub_meta["Metadata_Row"].map(lambda x: ord(x.upper())-65+1).astype(float).values
        cols_all   = sub_meta["Metadata_Column"].astype(int).astype(float).values
        rows_rep0  = rows_all[rep_mask]
        cols_rep0  = cols_all[rep_mask]
        corrected_plate = np.zeros_like(Vp)

        for j in range(Vp.shape[1]):
            y_rep0 = Vp[rep_mask, j]
            # drop extremes
            if exclude_extremes > 0 and y_rep0.size > 2*exclude_extremes:
                order = np.argsort(y_rep0)
                keep  = order[exclude_extremes:-exclude_extremes]
                rows_rep = rows_rep0[keep]
                cols_rep = cols_rep0[keep]
                y_rep    = y_rep0[keep]
            else:
                rows_rep, cols_rep, y_rep = rows_rep0, cols_rep0, y_rep0

            med     = np.median(y_rep)
            s_adapt = max(spline_s, len(y_rep) * np.var(y_rep))
            spline  = SmoothBivariateSpline(rows_rep, cols_rep, y_rep, s=s_adapt)
            bias    = spline(rows_all, cols_all, grid=False)

            corrected_plate[:, j] = Vp[:, j] - bias + med

        corrected[idx_plate, :] = corrected_plate

    return meta_cp, corrected

# def select_features(dframe_path, feat_selected_path):
#     '''Run feature selection'''
#     dframe = pd.read_parquet(dframe_path)
#     features = find_feat_cols(dframe.columns)
#     low_variance = variance_threshold(dframe, features)
#     features = [f for f in features if f not in low_variance]
#     logger.info(f'{len(low_variance)} features removed by variance_threshold')
#     high_corr = correlation_threshold(dframe, features)
#     features = [f for f in features if f not in high_corr]
#     logger.info(f'{len(high_corr)} features removed by correlation_threshold')

#     dframe.drop(columns=low_variance + high_corr, inplace=True)


def select_features(
    dframe_path: str,
    output_path: str,
    plate_col: str = "Metadata_Plate",
    drug_id_col: str = "Metadata_JCP2022",
    drug_id_for_replicates: str = "DMSO",
    spline_s: float = 1.0,
    exclude_extremes: int = 0
) -> None:
    """
    1) Load raw data.
    2) Select features by variance & correlation.
    3) Normalize plate effects on selected features.
    4) Save the final DataFrame.
    """
    # Load and split
    dframe = pd.read_parquet(dframe_path)
    features = find_feat_cols(dframe.columns)
    meta_cols = [c for c in dframe.columns if c not in features]

    # Raw feature selection
    low_var = variance_threshold(dframe, features)
    logger.info(f"{len(low_var)} features removed by variance_threshold")
    feats_kept = [f for f in features if f not in low_var]

    high_corr = correlation_threshold(dframe, feats_kept)
    logger.info(f"{len(high_corr)} features removed by correlation_threshold")
    final_feats = [f for f in feats_kept if f not in high_corr]

    # Subset to metadata + selected features
    df_selected = dframe[meta_cols + final_feats].copy()

    # Plate‐effect normalization
    meta = df_selected[meta_cols].reset_index(drop=True)
    vals = df_selected[final_feats].to_numpy()
    meta_norm, vals_norm = normalize_plate_effect_no_dask(
        meta, vals,
        plate_col=plate_col,
        drug_id_col=drug_id_col,
        drug_id_for_replicates=drug_id_for_replicates,
        spline_s=spline_s,
        exclude_extremes=exclude_extremes
    )

    # Reconstruct and save
    df_norm = pd.concat([
        meta_norm.reset_index(drop=True),
        pd.DataFrame(vals_norm, columns=final_feats)
    ], axis=1)
    df_norm.reset_index(drop=True).to_parquet(output_path)
    logger.info(f"Saved normalized selected features to {output_path}")