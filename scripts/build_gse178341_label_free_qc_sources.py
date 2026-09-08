#!/usr/bin/env python3
"""Create Figure 6A-F source tables from the label-free sketch."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
CACHE = ATLAS / "raw_csc_cache"
TARGET_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)


meta = pd.read_csv(
    ATLAS / "label_free_sketch_annotated_metadata.tsv.gz", sep="\t"
)
pcs = np.load(ATLAS / "label_free_sketch_pca.npy")
umap = np.load(ATLAS / "label_free_sketch_umap.npy")
arrays = np.load(ATLAS / "label_free_sketch_model_arrays.npz")
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
indptr = np.load(CACHE / "cell_indptr_int64.npy", mmap_mode="r")
gene_index = np.load(CACHE / "gene_indices_int32.npy", mmap_mode="r")
data = np.load(CACHE / "data_float32.npy", mmap_mode="r")
library_size = np.load(TARGET_CACHE)["library_size"].astype(float)

cells = meta["full_matrix_column"].to_numpy(dtype=int)
detected_features = np.diff(indptr)[cells].astype(float)
total_umi = library_size[cells]
mt_lookup = np.zeros(len(features), dtype=bool)
mt_lookup[
    features.loc[
        features["gene_symbol"].astype(str).str.startswith("MT-"),
        "feature_index",
    ].to_numpy(dtype=int)
] = True
mt_count = np.zeros(len(cells), dtype=float)
for row, cell in enumerate(cells):
    left, right = int(indptr[cell]), int(indptr[cell + 1])
    genes = gene_index[left:right]
    values = data[left:right]
    mt_count[row] = values[mt_lookup[genes]].sum()
mitochondrial_percent = 100 * mt_count / np.maximum(total_umi, 1)

pc_rows = []
for j in range(pcs.shape[1]):
    pc_rows.append(
        {
            "PC": j + 1,
            "standard_deviation": np.sqrt(
                arrays["pca_explained_variance"][j]
            ),
            "explained_variance": arrays["pca_explained_variance"][j],
            "explained_variance_ratio": arrays[
                "pca_explained_variance_ratio"
            ][j],
            "cumulative_explained_variance_ratio": arrays[
                "pca_explained_variance_ratio"
            ][: j + 1].sum(),
            "spearman_rho_log_total_UMI": spearmanr(
                pcs[:, j], np.log1p(total_umi)
            ).statistic,
            "spearman_rho_detected_features": spearmanr(
                pcs[:, j], detected_features
            ).statistic,
            "spearman_rho_mitochondrial_percent": spearmanr(
                pcs[:, j], mitochondrial_percent
            ).statistic,
        }
    )
pd.DataFrame(pc_rows).to_csv(
    ATLAS / "label_free_principal_component_statistics.tsv",
    sep="\t",
    index=False,
)

embedding = meta.copy()
embedding["PC1"] = pcs[:, 0]
embedding["PC2"] = pcs[:, 1]
embedding["UMAP1"] = umap[:, 0]
embedding["UMAP2"] = umap[:, 1]
embedding["total_UMI"] = total_umi
embedding["detected_features"] = detected_features
embedding["mitochondrial_percent"] = mitochondrial_percent
embedding.to_csv(
    ATLAS / "label_free_sketch_embedding_and_clusters.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
print("Completed label-free Figure 6A-F QC source tables")
