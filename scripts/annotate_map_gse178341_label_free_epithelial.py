#!/usr/bin/env python3
"""Select, annotate and map label-free GSE178341 epithelial subclusters."""

from pathlib import Path

import numpy as np
import pandas as pd
from pynndescent import NNDescent
from scipy import sparse
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/omics/crc_single_cell/GSE178341"
OUT = ROOT / "results/single_cell/GSE178341/label_free_epithelial"
SEED = 20260723
MARKERS = {
    "Stem/TA-like": ["OLFM4", "LGR5", "SMOC2", "SOX9", "ASCL2", "LEFTY1"],
    "Proliferating": ["MKI67", "TOP2A", "UBE2C", "CENPF", "STMN1", "TUBA1B"],
    "Immature goblet": ["SPINK4", "AGR2", "TFF3", "CLCA1", "MUC13"],
    "Enterocyte": ["CA1", "CA2", "APOA1", "FABP1", "KRT20", "ALDOB"],
    "Goblet": ["MUC2", "FCGBP", "ZG16", "CLCA1", "SPINK4", "TFF3"],
    "BEST4+": ["BEST4", "OTOP2", "CA7", "GUCA2A"],
    "Tuft": ["POU2F3", "TRPM5", "IL17RB", "AVIL", "SH2D6"],
    "Enteroendocrine": ["CHGA", "CHGB", "NEUROD1", "PAX6", "TPH1"],
}
MIN_CLUSTERS = 8
MAX_CLUSTERS = 15


meta = pd.read_csv(
    OUT / "epithelial_label_free_sketch_metadata.tsv.gz", sep="\t"
)
memberships = pd.read_csv(
    OUT / "epithelial_leiden_memberships.tsv.gz", sep="\t"
)
pcs = np.load(OUT / "epithelial_sketch_pca.npy")
rng = np.random.default_rng(SEED)
selection_rows = []
for column in memberships.columns[1:]:
    labels = memberships[column].to_numpy(dtype=int)
    clusters = np.unique(labels)
    if not MIN_CLUSTERS <= len(clusters) <= MAX_CLUSTERS:
        continue
    chosen = []
    for cluster in clusters:
        idx = np.flatnonzero(labels == cluster)
        if len(idx) > 600:
            idx = rng.choice(idx, 600, replace=False)
        chosen.append(idx)
    chosen = np.sort(np.concatenate(chosen))
    selection_rows.append(
        {
            "membership_column": column,
            "clusters": len(clusters),
            "silhouette_cells": len(chosen),
            "balanced_silhouette": silhouette_score(
                pcs[chosen, :30], labels[chosen]
            ),
        }
    )
selection = pd.DataFrame(selection_rows).sort_values(
    ["balanced_silhouette", "clusters"], ascending=[False, True]
)
selected_column = selection.iloc[0]["membership_column"]
selection["selected"] = selection["membership_column"].eq(selected_column)
selection.to_csv(
    OUT / "epithelial_leiden_silhouette_selection.tsv",
    sep="\t",
    index=False,
)
labels = memberships[selected_column].to_numpy(dtype=int)

features = pd.read_csv(
    ROOT / "results/single_cell/GSE178341/label_free_atlas/features.tsv.gz",
    sep="\t",
)
arrays = np.load(OUT / "epithelial_sketch_model_arrays.npz")
hvg_indices = arrays["hvg_indices"].astype(int)
hvg = features.set_index("feature_index").loc[hvg_indices].reset_index()
expression = sparse.load_npz(
    OUT / "epithelial_sketch_log_normalized.npz"
).tocsr()
gene_columns: dict[str, list[int]] = {}
for column, gene in enumerate(hvg["gene_symbol"].astype(str)):
    gene_columns.setdefault(gene, []).append(column)

marker_rows = []
for cluster in np.sort(np.unique(labels)):
    rows = np.flatnonzero(labels == cluster)
    local = expression[rows]
    for state, genes in MARKERS.items():
        for gene in genes:
            columns = gene_columns.get(gene, [])
            if not columns:
                continue
            values = np.asarray(local[:, columns].mean(axis=1)).ravel()
            marker_rows.append(
                {
                    "cluster": cluster,
                    "state_marker_set": state,
                    "gene": gene,
                    "cells": len(rows),
                    "mean_expression": values.mean(),
                    "detection_rate": np.mean(values > 0),
                }
            )
marker = pd.DataFrame(marker_rows)
for value in ["mean_expression", "detection_rate"]:
    marker[f"z_{value}"] = marker.groupby("gene", observed=True)[
        value
    ].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0)
        if x.std(ddof=0) > 0
        else 0
    )
marker["marker_evidence"] = (
    marker["z_mean_expression"] + marker["z_detection_rate"]
) / 2
score = (
    marker.groupby(["cluster", "state_marker_set"], observed=True)[
        "marker_evidence"
    ]
    .mean()
    .unstack(fill_value=-np.inf)
)
annotation = score.idxmax(axis=1).rename("marker_state").reset_index()
annotation["marker_score"] = score.max(axis=1).to_numpy()
annotation["marker_score_margin"] = score.apply(
    lambda x: x.sort_values(ascending=False).iloc[0]
    - x.sort_values(ascending=False).iloc[1],
    axis=1,
).to_numpy()
annotation["cells"] = annotation["cluster"].map(
    pd.Series(labels).value_counts()
)
annotation["subcluster_label"] = annotation.apply(
    lambda x: f"E{int(x['cluster']):02d} ({x['marker_state']})", axis=1
)
cluster_to_state = annotation.set_index("cluster")["marker_state"]
cluster_to_label = annotation.set_index("cluster")["subcluster_label"]
meta["epithelial_leiden_cluster"] = labels
meta["epithelial_marker_state"] = pd.Series(labels).map(
    cluster_to_state
).to_numpy()
meta["epithelial_subcluster_label"] = pd.Series(labels).map(
    cluster_to_label
).to_numpy()

# Publisher subtypes are appended only for post-hoc validation.
publisher = pd.read_csv(
    RAW / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
)
global_cells = meta["full_matrix_column"].to_numpy(dtype=int)
meta["publisher_epithelial_subcluster"] = publisher.iloc[global_cells][
    "cl295v11SubFull"
].to_numpy()
posthoc = pd.crosstab(
    meta["epithelial_subcluster_label"],
    meta["publisher_epithelial_subcluster"],
)
posthoc.to_csv(
    OUT / "epithelial_label_free_vs_publisher_subtype_confusion.tsv",
    sep="\t",
)
pd.DataFrame(
    {
        "metric": [
            "selected_membership",
            "selected_clusters",
            "balanced_silhouette",
            "posthoc_adjusted_rand_index",
            "posthoc_normalized_mutual_information",
        ],
        "value": [
            selected_column,
            len(np.unique(labels)),
            selection.iloc[0]["balanced_silhouette"],
            adjusted_rand_score(
                meta["publisher_epithelial_subcluster"],
                meta["epithelial_leiden_cluster"],
            ),
            normalized_mutual_info_score(
                meta["publisher_epithelial_subcluster"],
                meta["epithelial_leiden_cluster"],
            ),
        ],
    }
).to_csv(
    OUT / "epithelial_subcluster_posthoc_audit.tsv",
    sep="\t",
    index=False,
)

# Project all newly mapped epithelial cells and transfer the de novo labels.
full_expression = sparse.load_npz(
    OUT / "full_epithelial_hvg_log_normalized.npz"
).tocsr()
index = NNDescent(
    pcs[:, :30],
    n_neighbors=30,
    metric="cosine",
    random_state=SEED,
    n_jobs=-1,
    low_memory=True,
    compressed=True,
)
index.prepare()
n_epi = full_expression.shape[0]
mapped = np.empty(n_epi, dtype=np.int16)
confidence = np.empty(n_epi, dtype=np.float32)
unique_labels = np.sort(np.unique(labels))
for start in range(0, n_epi, 5000):
    end = min(start + 5000, n_epi)
    dense = full_expression[start:end].toarray().astype(np.float32, copy=False)
    dense -= arrays["gene_mean"]
    dense /= arrays["gene_sd"]
    np.clip(dense, -10, 10, out=dense)
    projected = (
        dense - arrays["pca_mean"]
    ) @ arrays["pca_components"][:30].T
    neighbor, distance = index.query(projected, k=15, epsilon=0.1)
    weight = 1 / np.maximum(distance, 1e-4)
    votes = np.column_stack(
        [
            np.sum(weight * (labels[neighbor] == value), axis=1)
            for value in unique_labels
        ]
    )
    best = np.argmax(votes, axis=1)
    mapped[start:end] = unique_labels[best]
    confidence[start:end] = votes[
        np.arange(len(best)), best
    ] / np.maximum(votes.sum(axis=1), 1e-12)

sketch_local = arrays["sketch_local"].astype(int)
mapped[sketch_local] = labels
confidence[sketch_local] = 1
full_meta = pd.read_csv(
    ROOT
    / "results/single_cell/GSE178341/label_free_atlas/"
    "full_cell_label_free_mapped_metadata.tsv.gz",
    sep="\t",
)
epi_global = arrays["epi_global"].astype(int)
epi = full_meta.iloc[epi_global].copy().reset_index(drop=True)
epi["epithelial_leiden_cluster"] = mapped
epi["epithelial_marker_state"] = pd.Series(mapped).map(
    cluster_to_state
).to_numpy()
epi["epithelial_subcluster_label"] = pd.Series(mapped).map(
    cluster_to_label
).to_numpy()
epi["epithelial_label_transfer_confidence"] = confidence
epi["publisher_epithelial_subcluster"] = publisher.iloc[epi_global][
    "cl295v11SubFull"
].to_numpy()
epi.to_csv(
    OUT / "full_epithelial_label_free_subclusters.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
annotation.to_csv(
    OUT / "epithelial_label_free_cluster_annotation.tsv",
    sep="\t",
    index=False,
)
marker.to_csv(
    OUT / "epithelial_label_free_marker_statistics.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
score.reset_index().to_csv(
    OUT / "epithelial_label_free_marker_scores.tsv",
    sep="\t",
    index=False,
)
meta.to_csv(
    OUT / "epithelial_label_free_sketch_annotated.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
print(
    f"Selected {selected_column}; mapped {n_epi:,} epithelial cells to "
    f"{len(unique_labels)} de novo subclusters"
)
