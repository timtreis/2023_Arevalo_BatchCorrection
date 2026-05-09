"""Thin wrappers around scib-metrics for the reference-mapping notebooks.

Uses scib_metrics >=0.5 API. Functions that require pre-computed kNN
accept a NeighborsResults object from scib_metrics.nearest_neighbors.
"""
from __future__ import annotations
import numpy as np
import scib_metrics
from scib_metrics.nearest_neighbors import NeighborsResults


def compound_asw(embedding: np.ndarray, compound_ids) -> float:
    """ASW over compound identities — higher = same-compound cells cluster."""
    return float(
        scib_metrics.silhouette_label(embedding, np.asarray(compound_ids))
    )


def batch_asw(embedding: np.ndarray, compound_ids, batch_ids) -> float:
    """Batch ASW — higher = better batch mixing within compound groups."""
    return float(
        scib_metrics.silhouette_batch(
            embedding, np.asarray(compound_ids), np.asarray(batch_ids)
        )
    )


def kbet_score(nn: NeighborsResults, batch_labels) -> float:
    """kBET acceptance rate — higher = better batch mixing."""
    return float(scib_metrics.kbet(nn, np.asarray(batch_labels))[0])


def ilisi_score(nn: NeighborsResults, batch_labels) -> float:
    """iLISI — higher = better batch mixing."""
    return float(scib_metrics.ilisi_knn(nn, np.asarray(batch_labels)))
