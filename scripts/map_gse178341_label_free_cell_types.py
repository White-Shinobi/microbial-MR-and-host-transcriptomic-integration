#!/usr/bin/env python3
"""Project all GSE178341 cells to the label-free sketch and transfer labels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pynndescent import NNDescent
from scipy import sparse
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/omics/crc_single_cell/GSE178341"
OUT = ROOT / "results/single_cell/GSE178341/label_free_atlas"
META = RAW / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
PUBLISHER = RAW / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
SEED = 20260723
CHUNK = 5000
K_NEIGHBORS = 15


def weighted_vote(
    neighbor_labels: np.ndarray,
    distances: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    unique_labels = np.sort(np.unique(neighbor_labels))
    weights = 1 / np.maximum(distances, 1e-4)
    vote = np.zeros((len(neighbor_labels), len(unique_labels)), dtype=np.float32)
    for j, label in enumerate(unique_labels):
        vote[:, j] = np.sum(weights * (neighbor_labels == label), axis=1)
    best = np.argmax(vote, axis=1)
    confidence = vote[np.arange(len(vote)), best] / np.maximum(
        vote.sum(axis=1), 1e-12
    )
    return unique_labels[best], confidence


metadata = pd.read_csv(META)
publisher = pd.read_csv(PUBLISHER)
if not np.array_equal(
    metadata["cellID"].to_numpy(dtype=object),
    publisher["sampleID"].to_numpy(dtype=object),
):
    raise RuntimeError("Metadata and publisher cell order differ")

sketch_meta = pd.read_csv(
    OUT / "label_free_sketch_annotated_metadata.tsv.gz", sep="\t"
)
arrays = np.load(OUT / "label_free_sketch_model_arrays.npz")
full_expression = sparse.load_npz(
    OUT / "full_cell_hvg_log_normalized_cells_by_genes.npz"
).tocsr()
sketch_pcs = np.load(OUT / "label_free_sketch_pca.npy").astype(
    np.float32, copy=False
)

sketch_cells = arrays["sketch_cells"].astype(int)
gene_mean = arrays["gene_mean"].astype(np.float32)
gene_sd = arrays["gene_sd"].astype(np.float32)
pca_components = arrays["pca_components"].astype(np.float32)
pca_mean = arrays["pca_mean"].astype(np.float32)
if not np.array_equal(
    sketch_cells, sketch_meta["full_matrix_column"].to_numpy(dtype=int)
):
    raise RuntimeError("Sketch model and annotated metadata are misaligned")

index = NNDescent(
    sketch_pcs[:, :30],
    n_neighbors=30,
    metric="cosine",
    random_state=SEED,
    n_jobs=-1,
    low_memory=True,
    compressed=True,
)
index.prepare()

n_cells = len(metadata)
mapped_cluster = np.empty(n_cells, dtype=np.int16)
mapping_confidence = np.empty(n_cells, dtype=np.float32)
sketch_cluster = sketch_meta["label_free_leiden_cluster"].to_numpy(dtype=int)

for start in range(0, n_cells, CHUNK):
    end = min(start + CHUNK, n_cells)
    dense = full_expression[start:end].toarray().astype(np.float32, copy=False)
    dense -= gene_mean
    dense /= gene_sd
    np.clip(dense, -10, 10, out=dense)
    projected = (dense - pca_mean) @ pca_components[:30].T
    neighbor_indices, distances = index.query(
        projected, k=K_NEIGHBORS, epsilon=0.1
    )
    labels, confidence = weighted_vote(
        sketch_cluster[neighbor_indices], distances
    )
    mapped_cluster[start:end] = labels.astype(np.int16)
    mapping_confidence[start:end] = confidence
    if end % 50_000 == 0 or end == n_cells:
        print(f"  mapped {end:,}/{n_cells:,} cells", flush=True)

# Preserve exact labels for the sketch cells themselves.
mapped_cluster[sketch_cells] = sketch_cluster.astype(np.int16)
mapping_confidence[sketch_cells] = 1.0

cluster_annotation = pd.read_csv(
    OUT / "label_free_cluster_annotation.tsv", sep="\t"
)
cluster_to_type = cluster_annotation.set_index("cluster")[
    "marker_assigned_cell_type"
]
mapped_type = pd.Series(mapped_cluster).map(cluster_to_type).to_numpy()
if pd.isna(mapped_type).any():
    raise RuntimeError("Some mapped clusters have no marker annotation")

publisher_map = {
    "Epi": "Epithelial",
    "TNKILC": "T/NK/ILC",
    "Myeloid": "Myeloid",
    "B": "B",
    "Plasma": "Plasma",
    "Strom": "Stromal",
    "Mast": "Mast",
}
publisher_standardized = publisher["clTopLevel"].map(publisher_map)

full = metadata.copy()
full.insert(0, "full_matrix_column", np.arange(n_cells, dtype=int))
full["label_free_leiden_cluster"] = mapped_cluster
full["label_free_marker_cell_type"] = mapped_type
full["label_transfer_confidence"] = mapping_confidence
full["selected_for_label_free_sketch"] = False
full.loc[sketch_cells, "selected_for_label_free_sketch"] = True
# Post-hoc fields: these were not supplied to the mapping model.
full["publisher_top_level"] = publisher["clTopLevel"].to_numpy()
full["publisher_cell_type_standardized"] = publisher_standardized.to_numpy()
full.to_csv(
    OUT / "full_cell_label_free_mapped_metadata.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

pd.crosstab(
    full["label_free_marker_cell_type"],
    full["publisher_cell_type_standardized"],
    margins=True,
).to_csv(
    OUT / "full_cell_label_free_vs_publisher_confusion.tsv", sep="\t"
)
pd.DataFrame(
    {
        "metric": [
            "full_cells",
            "sketch_cells",
            "nearest_neighbors",
            "median_label_transfer_confidence",
            "cells_with_confidence_below_0.5",
            "cell_level_exact_concordance",
            "adjusted_rand_index",
            "normalized_mutual_information",
        ],
        "value": [
            n_cells,
            len(sketch_cells),
            K_NEIGHBORS,
            np.median(mapping_confidence),
            np.sum(mapping_confidence < 0.5),
            np.mean(mapped_type == publisher_standardized.to_numpy()),
            adjusted_rand_score(publisher_standardized, mapped_type),
            normalized_mutual_info_score(publisher_standardized, mapped_type),
        ],
    }
).to_csv(
    OUT / "full_cell_label_transfer_audit.tsv",
    sep="\t",
    index=False,
)

print("Completed full-cell label-free projection and label transfer")
