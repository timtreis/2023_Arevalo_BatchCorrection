"""Thin wrappers around scib-metrics for the reference-mapping notebooks.

We keep the wrappers small so the metrics notebook stays readable. Any
heavy preprocessing (neighbours, joint embeddings) lives in the notebook,
not here.

NOTE: The scib-metrics API has churned across versions. The wrappers below
target the >=0.5 layout (`scib_metrics.mean_average_precision`, `kbet`,
`ilisi_graph`). If install resolves an older version, adjust here.
"""
from __future__ import annotations
import numpy as np
import scib_metrics


def compound_map(embedding: np.ndarray, compound_ids) -> float:
    """mAP over compound identities — high = same-compound cells cluster."""
    return float(
        scib_metrics.mean_average_precision(embedding, np.asarray(compound_ids))
    )


def moa_map(embedding: np.ndarray, moa_ids) -> float:
    """mAP over mechanism-of-action labels — high = MOA-coherent embedding."""
    return float(
        scib_metrics.mean_average_precision(embedding, np.asarray(moa_ids))
    )


def kbet_score(embedding: np.ndarray, batch_labels) -> float:
    """kBET acceptance rate on `embedding` against `batch_labels`."""
    return float(
        scib_metrics.kbet(embedding, np.asarray(batch_labels))
    )


def ilisi_score(embedding: np.ndarray, batch_labels) -> float:
    """iLISI on `embedding` — higher = better mixing across batches."""
    return float(
        scib_metrics.ilisi_graph(embedding, np.asarray(batch_labels))
    )
