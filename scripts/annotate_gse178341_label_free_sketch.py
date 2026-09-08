#!/usr/bin/env python3
"""Select a label-free Leiden solution and annotate it with canonical markers.

Publisher annotations are used only after marker-based labels have been fixed,
to quantify post-hoc concordance.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/single_cell/GSE178341/label_free_atlas"
SEED = 20260723
MIN_CLUSTERS = 8
MAX_CLUSTERS = 22
SILHOUETTE_CAP_PER_CLUSTER = 500

MARKERS = {
    "Epithelial": [
        "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "KRT17",
        "CEACAM5", "CEACAM6",
    ],
    "T/NK/ILC": [
        "CD3D", "CD3E", "TRBC1", "NKG7", "GNLY", "KLRD1", "IL7R",
    ],
    "Myeloid": [
        "LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "CD68", "CSF1R",
    ],
    "B": [
        "MS4A1", "CD79A", "CD79B", "CD37", "CD74", "CD19", "CD22",
    ],
    "Plasma": [
        "MZB1", "JCHAIN", "SDC1", "IGKC", "XBP1", "PRDM1", "TNFRSF17",
    ],
    "Stromal": [
        "COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "COL6A2",
        "VWF", "EMCN", "PECAM1", "RGS5", "ACTA2", "PDGFRB",
    ],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "MS4A2", "CPA3", "HDC"],
}


def balanced_silhouette_indices(
    labels: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    chosen: list[np.ndarray] = []
    for label in np.unique(labels):
        idx = np.flatnonzero(labels == label)
        if len(idx) > SILHOUETTE_CAP_PER_CLUSTER:
            idx = rng.choice(
                idx, SILHOUETTE_CAP_PER_CLUSTER, replace=False
            )
        chosen.append(np.sort(idx))
    return np.sort(np.concatenate(chosen))


meta = pd.read_csv(OUT / "label_free_sketch_metadata.tsv.gz", sep="\t")
memberships = pd.read_csv(
    OUT / "label_free_sketch_leiden_memberships.tsv.gz", sep="\t"
)
pcs = np.load(OUT / "label_free_sketch_pca.npy", mmap_mode="r")
resolution_audit = pd.read_csv(
    OUT / "label_free_sketch_leiden_resolution_audit.tsv", sep="\t"
)
rng = np.random.default_rng(SEED)

silhouette_rows = []
for column in memberships.columns:
    if column == "sketch_row":
        continue
    labels = memberships[column].to_numpy()
    n_clusters = np.unique(labels).size
    if not MIN_CLUSTERS <= n_clusters <= MAX_CLUSTERS:
        continue
    idx = balanced_silhouette_indices(labels, rng)
    value = silhouette_score(
        np.asarray(pcs[idx, :30]),
        labels[idx],
        metric="euclidean",
    )
    silhouette_rows.append(
        {
            "membership_column": column,
            "clusters": n_clusters,
            "silhouette_cells": len(idx),
            "balanced_silhouette": value,
        }
    )

silhouette_audit = pd.DataFrame(silhouette_rows).sort_values(
    ["balanced_silhouette", "clusters"], ascending=[False, True]
)
if silhouette_audit.empty:
    raise RuntimeError("No Leiden solution met the prespecified cluster range")
selected_column = silhouette_audit.iloc[0]["membership_column"]
selected_labels = memberships[selected_column].to_numpy(dtype=int)
resolution = float(selected_column.replace("leiden_r", ""))
silhouette_audit["selected"] = (
    silhouette_audit["membership_column"] == selected_column
)
silhouette_audit = silhouette_audit.merge(
    resolution_audit, on="clusters", how="left", suffixes=("", "_grid")
)
silhouette_audit.to_csv(
    OUT / "label_free_leiden_silhouette_selection.tsv",
    sep="\t",
    index=False,
)

features = pd.read_csv(OUT / "features.tsv.gz", sep="\t")
arrays = np.load(OUT / "label_free_sketch_model_arrays.npz")
hvg_indices = arrays["hvg_feature_indices"].astype(int)
hvg_features = features.set_index("feature_index").loc[hvg_indices].reset_index()
expression = sparse.load_npz(
    OUT / "label_free_sketch_log_normalized_cells_by_hvgs.npz"
).tocsr()

gene_to_columns: dict[str, list[int]] = {}
for column, symbol in enumerate(hvg_features["gene_symbol"].astype(str)):
    gene_to_columns.setdefault(symbol, []).append(column)

available_markers = []
for cell_type, symbols in MARKERS.items():
    for symbol in symbols:
        if symbol in gene_to_columns:
            available_markers.append(
                {
                    "cell_type": cell_type,
                    "gene_symbol": symbol,
                    "hvg_columns": ",".join(
                        str(x) for x in gene_to_columns[symbol]
                    ),
                }
            )
pd.DataFrame(available_markers).to_csv(
    OUT / "label_free_marker_availability.tsv", sep="\t", index=False
)

clusters = np.sort(np.unique(selected_labels))
marker_rows = []
for cluster in clusters:
    rows = np.flatnonzero(selected_labels == cluster)
    local = expression[rows]
    for cell_type, symbols in MARKERS.items():
        for symbol in symbols:
            columns = gene_to_columns.get(symbol, [])
            if not columns:
                continue
            values = local[:, columns].mean(axis=1)
            values = np.asarray(values).ravel()
            marker_rows.append(
                {
                    "cluster": cluster,
                    "cell_type_marker_set": cell_type,
                    "gene_symbol": symbol,
                    "cells": len(rows),
                    "mean_log_normalized_expression": values.mean(),
                    "detection_rate": np.mean(values > 0),
                }
            )
marker_table = pd.DataFrame(marker_rows)

# Standardize each marker across clusters, then average the expression and
# detection evidence within each canonical marker set.
for value_column in [
    "mean_log_normalized_expression",
    "detection_rate",
]:
    marker_table[f"z_{value_column}"] = marker_table.groupby(
        "gene_symbol", observed=True
    )[value_column].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0)
        if x.std(ddof=0) > 0
        else 0.0
    )
marker_table["marker_evidence"] = (
    marker_table["z_mean_log_normalized_expression"]
    + marker_table["z_detection_rate"]
) / 2
score_table = (
    marker_table.groupby(
        ["cluster", "cell_type_marker_set"], observed=True
    )["marker_evidence"]
    .mean()
    .unstack(fill_value=-np.inf)
)
cluster_annotation = score_table.idxmax(axis=1).rename(
    "marker_assigned_cell_type"
).reset_index()
cluster_annotation["runner_up_cell_type"] = score_table.apply(
    lambda row: row.sort_values(ascending=False).index[1], axis=1
).to_numpy()
cluster_annotation["marker_score"] = score_table.max(axis=1).to_numpy()
cluster_annotation["marker_score_margin"] = score_table.apply(
    lambda row: row.sort_values(ascending=False).iloc[0]
    - row.sort_values(ascending=False).iloc[1],
    axis=1,
).to_numpy()
cluster_annotation["cells"] = cluster_annotation["cluster"].map(
    pd.Series(selected_labels).value_counts()
)

marker_table.to_csv(
    OUT / "label_free_cluster_marker_statistics.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
score_table.reset_index().to_csv(
    OUT / "label_free_cluster_cell_type_scores.tsv",
    sep="\t",
    index=False,
)

cluster_to_type = cluster_annotation.set_index("cluster")[
    "marker_assigned_cell_type"
]
meta["label_free_leiden_cluster"] = selected_labels
meta["label_free_marker_cell_type"] = pd.Series(selected_labels).map(
    cluster_to_type
).to_numpy()

# Post-hoc validation begins here; publisher labels played no role above.
publisher_map = {
    "Epi": "Epithelial",
    "TNKILC": "T/NK/ILC",
    "Myeloid": "Myeloid",
    "B": "B",
    "Plasma": "Plasma",
    "Strom": "Stromal",
    "Mast": "Mast",
}
meta["publisher_cell_type_standardized"] = meta[
    "publisher_top_level"
].map(publisher_map)
confusion = pd.crosstab(
    meta["label_free_marker_cell_type"],
    meta["publisher_cell_type_standardized"],
    margins=True,
)
confusion.to_csv(
    OUT / "label_free_marker_vs_publisher_confusion.tsv", sep="\t"
)
concordance = pd.DataFrame(
    {
        "metric": [
            "selected_resolution",
            "selected_clusters",
            "balanced_silhouette",
            "cell_level_exact_concordance",
            "adjusted_rand_index",
            "normalized_mutual_information",
        ],
        "value": [
            resolution,
            np.unique(selected_labels).size,
            silhouette_audit.loc[
                silhouette_audit["selected"], "balanced_silhouette"
            ].iloc[0],
            np.mean(
                meta["label_free_marker_cell_type"]
                == meta["publisher_cell_type_standardized"]
            ),
            adjusted_rand_score(
                meta["publisher_cell_type_standardized"],
                meta["label_free_marker_cell_type"],
            ),
            normalized_mutual_info_score(
                meta["publisher_cell_type_standardized"],
                meta["label_free_marker_cell_type"],
            ),
        ],
    }
)
concordance.to_csv(
    OUT / "label_free_annotation_posthoc_concordance.tsv",
    sep="\t",
    index=False,
)
cluster_annotation.to_csv(
    OUT / "label_free_cluster_annotation.tsv", sep="\t", index=False
)
meta.to_csv(
    OUT / "label_free_sketch_annotated_metadata.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

print(
    f"Selected {selected_column} with "
    f"{np.unique(selected_labels).size} clusters; marker annotation complete"
)
