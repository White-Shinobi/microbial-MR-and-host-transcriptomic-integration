#!/usr/bin/env python3
"""Prepare coherent Figure 7 inputs from label-free epithelial subclusters."""

from pathlib import Path
import gzip

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from scipy.sparse.csgraph import dijkstra
from scipy.stats import spearmanr, wilcoxon
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
EPI = ROOT / "results/single_cell/GSE178341/label_free_epithelial"
OUT = ROOT / "results/single_cell/GSE178341/label_free_figure7"
SOURCE = ROOT / "manuscript/source_data"
CACHE = ATLAS / "raw_csc_cache"
TARGET_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)
BIOMARKERS = ["GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"]
MATURE_STATES = {"Enterocyte", "BEST4+", "Goblet"}
CELLCHAT_CAP = 5
SEED = 20260723


def safe_wilcoxon(tumor: np.ndarray, normal: np.ndarray) -> float:
    diff = tumor - normal
    if len(diff) < 3:
        return np.nan
    if np.allclose(diff, 0):
        return 1.0
    return float(
        wilcoxon(
            tumor, normal, alternative="two-sided", method="auto"
        ).pvalue
    )


OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)
sketch = pd.read_csv(
    EPI / "epithelial_label_free_sketch_annotated.tsv.gz", sep="\t"
)
full_epi = pd.read_csv(
    EPI / "full_epithelial_label_free_subclusters.tsv.gz", sep="\t"
)
pcs = np.load(EPI / "epithelial_sketch_pca.npy")
umap = np.load(EPI / "epithelial_sketch_umap.npy")
sketch["UMAP1"] = umap[:, 0]
sketch["UMAP2"] = umap[:, 1]
sketch["PC1"] = pcs[:, 0]
sketch["PC2"] = pcs[:, 1]

# Patient-level composition of the de novo epithelial subclusters.
subclusters = sorted(full_epi["epithelial_subcluster_label"].unique())
patient_tissue = full_epi[["PID", "SPECIMEN_TYPE"]].drop_duplicates()
grid = (
    patient_tissue.assign(_key=1)
    .merge(
        pd.DataFrame(
            {"epithelial_subcluster_label": subclusters, "_key": 1}
        ),
        on="_key",
    )
    .drop(columns="_key")
)
observed = (
    full_epi.groupby(
        ["PID", "SPECIMEN_TYPE", "epithelial_subcluster_label"],
        observed=True,
    )
    .size()
    .rename("n_cells")
    .reset_index()
)
composition = grid.merge(
    observed,
    on=["PID", "SPECIMEN_TYPE", "epithelial_subcluster_label"],
    how="left",
)
composition["n_cells"] = composition["n_cells"].fillna(0).astype(int)
composition["total_epithelial_cells"] = composition.groupby(
    ["PID", "SPECIMEN_TYPE"]
)["n_cells"].transform("sum")
composition["fraction"] = (
    composition["n_cells"] / composition["total_epithelial_cells"]
)
stats = []
for label in subclusters:
    wide = composition.loc[
        composition["epithelial_subcluster_label"].eq(label)
    ].pivot(index="PID", columns="SPECIMEN_TYPE", values="fraction")
    wide = wide.dropna(subset=["N", "T"])
    stats.append(
        {
            "epithelial_subcluster_label": label,
            "n_pairs": len(wide),
            "median_normal": wide["N"].median(),
            "median_CRC": wide["T"].median(),
            "median_paired_difference_CRC_minus_normal": (
                wide["T"] - wide["N"]
            ).median(),
            "paired_wilcoxon_P": safe_wilcoxon(
                wide["T"].to_numpy(), wide["N"].to_numpy()
            ),
        }
    )
stats = pd.DataFrame(stats)
valid = stats["paired_wilcoxon_P"].notna()
stats["paired_wilcoxon_BH_FDR"] = np.nan
stats.loc[valid, "paired_wilcoxon_BH_FDR"] = multipletests(
    stats.loc[valid, "paired_wilcoxon_P"], method="fdr_bh"
)[1]
composition.to_csv(
    OUT / "label_free_epithelial_subcluster_patient_composition.tsv",
    sep="\t",
    index=False,
)
stats.to_csv(
    OUT / "label_free_epithelial_subcluster_paired_statistics.tsv",
    sep="\t",
    index=False,
)

# Exploratory graph pseudotime on the label-free sketch. The root state is
# selected among mature marker states by the largest normal enrichment.
annotation = pd.read_csv(
    EPI / "epithelial_label_free_cluster_annotation.tsv", sep="\t"
)
cluster_state = annotation.set_index("cluster")["marker_state"]
counts_by_tissue = pd.crosstab(
    sketch["epithelial_leiden_cluster"], sketch["SPECIMEN_TYPE"]
)
counts_by_tissue = counts_by_tissue.reindex(
    columns=["N", "T"], fill_value=0
)
candidate_clusters = [
    cluster
    for cluster, state in cluster_state.items()
    if state in MATURE_STATES
]
if not candidate_clusters:
    candidate_clusters = list(cluster_state.index)
normal_enrichment = (
    (counts_by_tissue.loc[candidate_clusters, "N"] + 0.5)
    / (counts_by_tissue.loc[candidate_clusters, "T"] + 0.5)
)
root_cluster = int(normal_enrichment.idxmax())
root_candidates = np.flatnonzero(
    sketch["epithelial_leiden_cluster"].to_numpy() == root_cluster
)
centroid = pcs[root_candidates, :30].mean(axis=0)
root_index = int(
    root_candidates[
        np.argmin(
            np.sum(
                (pcs[root_candidates, :30] - centroid) ** 2, axis=1
            )
        )
    ]
)
neighbors = NearestNeighbors(
    n_neighbors=31, metric="euclidean", n_jobs=-1
).fit(pcs[:, :30])
distance, index = neighbors.kneighbors(pcs[:, :30])
rows = np.repeat(np.arange(len(sketch)), 30)
graph = sparse.csr_matrix(
    (
        distance[:, 1:].reshape(-1),
        (rows, index[:, 1:].reshape(-1)),
    ),
    shape=(len(sketch), len(sketch)),
)
graph = graph.maximum(graph.T)
pt = dijkstra(graph, directed=False, indices=root_index)
if not np.isfinite(pt).all():
    raise RuntimeError("Epithelial sketch graph is disconnected")
pt = 100 * (pt - pt.min()) / max(pt.max() - pt.min(), 1e-8)
sketch["graph_pseudotime"] = pt
sketch["root_epithelial_cluster"] = root_cluster
sketch.to_csv(
    OUT / "label_free_epithelial_embedding_pseudotime.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

patient_pt = (
    sketch.groupby(["PID", "SPECIMEN_TYPE"], observed=True)
    .agg(
        median_pseudotime=("graph_pseudotime", "median"),
        mean_pseudotime=("graph_pseudotime", "mean"),
        n_cells=("cellID", "size"),
    )
    .reset_index()
)
wide = patient_pt.pivot(
    index="PID", columns="SPECIMEN_TYPE", values="median_pseudotime"
).dropna(subset=["N", "T"])
pt_stats = pd.DataFrame(
    [
        {
            "n_pairs": len(wide),
            "median_normal": patient_pt.loc[
                patient_pt["SPECIMEN_TYPE"].eq("N"),
                "median_pseudotime",
            ].median(),
            "median_CRC": patient_pt.loc[
                patient_pt["SPECIMEN_TYPE"].eq("T"),
                "median_pseudotime",
            ].median(),
            "median_paired_difference_CRC_minus_normal": (
                wide["T"] - wide["N"]
            ).median(),
            "paired_wilcoxon_P": safe_wilcoxon(
                wide["T"].to_numpy(), wide["N"].to_numpy()
            ),
            "root_cluster": root_cluster,
            "root_state": cluster_state.loc[root_cluster],
        }
    ]
)
patient_pt.to_csv(
    OUT / "label_free_epithelial_patient_pseudotime.tsv",
    sep="\t",
    index=False,
)
pt_stats.to_csv(
    OUT / "label_free_epithelial_patient_pseudotime_statistics.tsv",
    sep="\t",
    index=False,
)

target = np.load(TARGET_CACHE)
library = target["library_size"]
sketch_global = sketch["full_matrix_column"].to_numpy(dtype=int)
biomarker_rows = []
for i, gene in enumerate(BIOMARKERS):
    expression = np.log1p(
        target["target_counts"][i, sketch_global]
        / np.maximum(library[sketch_global], 1)
        * 10_000
    )
    rho, p_value = spearmanr(expression, pt)
    biomarker_rows.append(
        {
            "gene": gene,
            "sketch_cells": len(expression),
            "detected_cells": np.count_nonzero(expression),
            "spearman_rho": rho,
            "spearman_P": p_value,
        }
    )
biomarker_pt = pd.DataFrame(biomarker_rows)
biomarker_pt["spearman_BH_FDR"] = multipletests(
    biomarker_pt["spearman_P"], method="fdr_bh"
)[1]
biomarker_pt.to_csv(
    OUT / "label_free_biomarker_epithelial_pseudotime.tsv",
    sep="\t",
    index=False,
)

# CellChat: epithelial identities are the de novo subclusters, while other
# cells use the newly mapped label-free broad cell types.
all_cells = pd.read_csv(
    ATLAS / "full_cell_label_free_mapped_metadata.tsv.gz", sep="\t"
)
epi_identity = full_epi.set_index("full_matrix_column")[
    "epithelial_subcluster_label"
]
all_cells["cellchat_identity"] = all_cells[
    "label_free_marker_cell_type"
]
is_epi = all_cells["label_free_marker_cell_type"].eq("Epithelial")
all_cells.loc[is_epi, "cellchat_identity"] = all_cells.loc[
    is_epi, "full_matrix_column"
].map(epi_identity)
selected_parts = []
for indices in all_cells.groupby(
    ["PID", "SPECIMEN_TYPE", "cellchat_identity"],
    observed=True,
    sort=True,
).indices.values():
    indices = np.asarray(indices, dtype=int)
    if len(indices) > CELLCHAT_CAP:
        indices = rng.choice(indices, CELLCHAT_CAP, replace=False)
    selected_parts.append(np.sort(indices))
selected = np.sort(np.concatenate(selected_parts))
selected_meta = all_cells.iloc[selected].copy().reset_index(drop=True)

data = np.load(CACHE / "data_float32.npy", mmap_mode="r")
genes = np.load(CACHE / "gene_indices_int32.npy", mmap_mode="r")
indptr = np.load(CACHE / "cell_indptr_int64.npy", mmap_mode="r")
data_parts = []
gene_parts = []
selected_indptr = np.zeros(len(selected) + 1, dtype=np.int64)
for i, cell in enumerate(selected):
    left, right = int(indptr[cell]), int(indptr[cell + 1])
    data_parts.append(np.asarray(data[left:right], dtype=np.float32))
    gene_parts.append(np.asarray(genes[left:right], dtype=np.int32))
    selected_indptr[i + 1] = selected_indptr[i] + (right - left)
selected_matrix = sparse.csc_matrix(
    (
        np.concatenate(data_parts),
        np.concatenate(gene_parts),
        selected_indptr,
    ),
    shape=(len(pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")), len(selected)),
)
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
valid_symbol = features["gene_symbol"].notna() & features[
    "gene_symbol"
].ne("")
symbols = features.loc[valid_symbol, "gene_symbol"].astype(str)
unique_symbols = np.sort(symbols.unique())
symbol_index = {symbol: i for i, symbol in enumerate(unique_symbols)}
feature_rows = np.flatnonzero(valid_symbol.to_numpy())
collapse = sparse.csr_matrix(
    (
        np.ones(len(feature_rows), dtype=np.float32),
        (
            np.asarray(
                [symbol_index[symbol] for symbol in symbols], dtype=int
            ),
            feature_rows,
        ),
    ),
    shape=(len(unique_symbols), len(features)),
)
collapsed = (collapse @ selected_matrix).tocoo()
with gzip.open(
    OUT / "label_free_cellchat_counts_genes_by_cells.mtx.gz", "wb"
) as handle:
    mmwrite(handle, collapsed, field="integer")
pd.DataFrame({"gene_symbol": unique_symbols}).to_csv(
    OUT / "label_free_cellchat_gene_symbols.tsv",
    sep="\t",
    index=False,
)
selected_meta[
    [
        "cellID",
        "PID",
        "SPECIMEN_TYPE",
        "label_free_marker_cell_type",
        "cellchat_identity",
    ]
].to_csv(
    OUT / "label_free_cellchat_cell_metadata.tsv",
    sep="\t",
    index=False,
)
selected_meta.groupby(
    ["SPECIMEN_TYPE", "cellchat_identity"], observed=True
).size().rename("sampled_cells").reset_index().to_csv(
    OUT / "label_free_cellchat_sampling_audit.tsv",
    sep="\t",
    index=False,
)
pd.DataFrame(
    {
        "metric": [
            "mapped_epithelial_cells",
            "epithelial_sketch_cells",
            "label_free_epithelial_subclusters",
            "matched_patients",
            "root_cluster",
            "root_state",
            "cellchat_cells",
            "cellchat_cap_per_patient_tissue_identity",
        ],
        "value": [
            len(full_epi),
            len(sketch),
            len(subclusters),
            len(wide),
            root_cluster,
            cluster_state.loc[root_cluster],
            len(selected),
            CELLCHAT_CAP,
        ],
    }
).to_csv(OUT / "label_free_figure7_audit.tsv", sep="\t", index=False)
print(
    f"Prepared Figure 7: {len(subclusters)} epithelial subclusters, "
    f"{len(wide)} matched patients, {len(selected):,} CellChat cells"
)
