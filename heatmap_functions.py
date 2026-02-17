# heatmap_functions.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
except Exception:
    sp = None


Metric = Literal["mean", "pct"]
Scale = Literal["raw", "zscore"]


@dataclass
class HeatmapResult:
    mean_expr: pd.DataFrame  # index=clusters, columns=genes
    pct_expr: pd.DataFrame   # index=clusters, columns=genes
    clusters: list[str]
    genes: list[str]


def _to_dense_2d(X) -> np.ndarray:
    """
    Convert sparse/dense matrix to a dense np.ndarray (cells x genes).
    Use only after subsetting to a small set of genes.
    """
    if sp is not None and sp.issparse(X):
        return X.toarray()
    return np.asarray(X)


def _get_expr_matrix(
    adata,
    genes: list[str],
    *,
    use_layer: Optional[str] = None,
    use_raw: bool = False,
) -> tuple[np.ndarray, list[str]]:
    """
    Returns:
      X: (n_cells, n_genes) dense float array
      genes_kept: genes actually found
    """
    genes = [str(g) for g in genes if str(g).strip()]
    if not genes:
        return np.zeros((adata.n_obs, 0), dtype=float), []

    # choose source
    src = adata
    if use_raw and getattr(adata, "raw", None) is not None:
        src = adata.raw

    # filter to genes that exist
    var_names = getattr(src, "var_names", None)
    if var_names is None:
        raise ValueError("adata.var_names is missing")

    genes_kept = [g for g in genes if g in var_names]
    if not genes_kept:
        return np.zeros((adata.n_obs, 0), dtype=float), []

    # subset expression
    if use_layer is not None:
        if not hasattr(adata, "layers") or use_layer not in adata.layers:
            raise KeyError(f'Layer "{use_layer}" not found in adata.layers')
        X = adata[:, genes_kept].layers[use_layer]
    else:
        X = src[:, genes_kept].X

    X = _to_dense_2d(X).astype(float, copy=False)
    return X, genes_kept


def compute_cluster_gene_heatmap(
    adata,
    genes: list[str],
    *,
    cluster_key: str = "leiden",
    use_layer: Optional[str] = None,
    use_raw: bool = False,
    expr_threshold: float = 0.0,
) -> HeatmapResult:
    """
    Heavy compute:
      - extracts per-cell expression for selected genes
      - aggregates mean expression and percent expressing by cluster
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f'Cluster key "{cluster_key}" not found in adata.obs')

    clusters_series = adata.obs[cluster_key].astype(str)
    clusters = sorted(clusters_series.unique().tolist(), key=lambda x: (0, int(x)) if x.isdigit() else (1, x))

    X, genes_kept = _get_expr_matrix(adata, genes, use_layer=use_layer, use_raw=use_raw)

    if X.shape[1] == 0:
        # no genes found
        mean_df = pd.DataFrame(index=clusters, columns=[], dtype=float)
        pct_df = pd.DataFrame(index=clusters, columns=[], dtype=float)
        return HeatmapResult(mean_df, pct_df, clusters, genes_kept)

    # build outputs
    mean_df = pd.DataFrame(index=clusters, columns=genes_kept, dtype=float)
    pct_df = pd.DataFrame(index=clusters, columns=genes_kept, dtype=float)

    # vectorized per cluster
    # (this is “heavy” but still efficient for typical gene panel sizes)
    for c in clusters:
        mask = (clusters_series.values == c)
        if not np.any(mask):
            mean_df.loc[c, :] = np.nan
            pct_df.loc[c, :] = np.nan
            continue

        sub = X[mask, :]  # cells_in_cluster x genes
        mean_df.loc[c, :] = np.mean(sub, axis=0)

        if expr_threshold <= 0.0:
            pct_df.loc[c, :] = np.mean(sub > 0.0, axis=0) * 100.0
        else:
            pct_df.loc[c, :] = np.mean(sub > float(expr_threshold), axis=0) * 100.0

    return HeatmapResult(mean_df, pct_df, clusters, genes_kept)


def prepare_heatmap_matrix(
    result: HeatmapResult,
    *,
    metric: Metric = "mean",
    scale: Scale = "raw",
) -> pd.DataFrame:
    """
    Light-ish transform:
      metric='mean' or 'pct'
      scale='raw' or 'zscore' (zscore per gene across clusters)
    """
    df = result.mean_expr if metric == "mean" else result.pct_expr
    df = df.copy()

    if scale == "zscore":
        # z-score each gene column across clusters
        vals = df.to_numpy(dtype=float)
        mu = np.nanmean(vals, axis=0)
        sigma = np.nanstd(vals, axis=0)
        sigma[sigma == 0] = 1.0
        z = (vals - mu) / sigma
        df.iloc[:, :] = z

    return df
