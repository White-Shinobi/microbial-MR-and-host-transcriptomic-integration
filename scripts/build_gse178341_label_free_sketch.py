#!/usr/bin/env python3
"""Build a label-free, patient-tissue-balanced GSE178341 cell sketch.

No publisher cell-type label is used for HVG selection, leverage-score
calculation, sampling, PCA, neighborhood construction, or UMAP. Publisher
annotations are copied only into the final audit table for post-hoc validation.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
import umap


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/omics/crc_single_cell/GSE178341"
OUT = ROOT / "results/single_cell/GSE178341/label_free_atlas"
H5 = RAW / "GSE178341_crc10x_full_c295v4_submit.h5"
META = RAW / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
PUBLISHER = RAW / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
LIB_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)
CSC_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/label_free_atlas/"
    "raw_csc_cache"
)

SEED = 20260723
HVG_ESTIMATION_CAP_PER_PATIENT_TISSUE = 200
SKETCH_CAP_PER_PATIENT_TISSUE = 500
N_HVG = 2000
N_LEVERAGE_DIMS = 50
N_PCS = 50
RETAINED_PCS = 30
H5_CELL_BLOCK = 1000


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            value.decode("utf-8").rstrip("\x00")
            if isinstance(value, bytes)
            else str(value)
            for value in values
        ],
        dtype=object,
    )


def sample_by_group(
    groups: pd.Series,
    cap: int,
    rng: np.random.Generator,
    probabilities: np.ndarray | None = None,
) -> np.ndarray:
    selected: list[np.ndarray] = []
    group_indices = groups.groupby(groups, sort=True).indices
    for _, indices in group_indices.items():
        indices = np.asarray(indices, dtype=np.int64)
        if len(indices) <= cap:
            chosen = indices
        else:
            prob = None
            if probabilities is not None:
                prob = np.maximum(probabilities[indices].astype(float), 0)
                if np.isfinite(prob).all() and prob.sum() > 0:
                    prob = prob / prob.sum()
                else:
                    prob = None
            chosen = rng.choice(indices, cap, replace=False, p=prob)
        selected.append(np.sort(chosen))
    return np.sort(np.concatenate(selected))


def estimate_hvg_moments_streaming(
    data_all: np.ndarray,
    gene_all: np.ndarray,
    indptr_all: np.ndarray,
    selected_cells: np.ndarray,
    library_size: np.ndarray,
    n_genes: int,
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate gene moments without materializing cells-by-all-genes.

    The HDF5 matrix is scanned in contiguous cell blocks. Only entries belonging
    to the prespecified label-free HVG-estimation sample contribute to the
    moments. This preserves the calculation while keeping peak memory bounded.
    """
    selected_mask = np.zeros(n_cells, dtype=bool)
    selected_mask[selected_cells] = True
    sum_log = np.zeros(n_genes, dtype=np.float64)
    sumsq_log = np.zeros(n_genes, dtype=np.float64)
    detected = np.zeros(n_genes, dtype=np.int64)
    processed = 0

    for cell_start in range(0, n_cells, H5_CELL_BLOCK):
        cell_end = min(cell_start + H5_CELL_BLOCK, n_cells)
        local_selected = selected_mask[cell_start:cell_end]
        if not local_selected.any():
            continue
        data_start = int(indptr_all[cell_start])
        data_end = int(indptr_all[cell_end])
        raw_gene = gene_all[data_start:data_end]
        raw_data = data_all[data_start:data_end]
        counts_per_cell = np.diff(indptr_all[cell_start : cell_end + 1])
        raw_cell = np.repeat(
            np.arange(cell_end - cell_start, dtype=np.int32), counts_per_cell
        )
        keep = local_selected[raw_cell]
        kept_cell = raw_cell[keep] + cell_start
        values = np.log1p(
            raw_data[keep]
            * (10_000 / np.maximum(library_size[kept_cell], 1))
        )
        genes = raw_gene[keep]
        sum_log += np.bincount(genes, weights=values, minlength=n_genes)
        sumsq_log += np.bincount(
            genes, weights=values * values, minlength=n_genes
        )
        detected += np.bincount(genes, minlength=n_genes)
        processed += int(local_selected.sum())
        if processed and processed % 5000 < int(local_selected.sum()):
            print(
                f"  accumulated HVG moments for {processed:,}/"
                f"{len(selected_cells):,} selected cells",
                flush=True,
            )

    if processed != len(selected_cells):
        raise RuntimeError(
            f"Processed {processed} HVG-estimation cells; expected "
            f"{len(selected_cells)}"
        )
    mean = sum_log / processed
    variance = np.maximum(sumsq_log / processed - mean**2, 0)
    detection = detected.astype(np.float64) / processed
    return mean, variance, detection


def extract_hvg_all_cells(
    data_all: np.ndarray,
    gene_all: np.ndarray,
    indptr_all: np.ndarray,
    hvg_indices: np.ndarray,
    n_genes: int,
    n_cells: int,
) -> sparse.csr_matrix:
    """Stream all cells while retaining only selected gene rows."""
    lookup = np.full(n_genes, -1, dtype=np.int32)
    lookup[hvg_indices] = np.arange(len(hvg_indices), dtype=np.int32)
    blocks: list[sparse.csr_matrix] = []
    for cell_start in range(0, n_cells, H5_CELL_BLOCK):
        cell_end = min(cell_start + H5_CELL_BLOCK, n_cells)
        data_start = int(indptr_all[cell_start])
        data_end = int(indptr_all[cell_end])
        raw_gene = gene_all[data_start:data_end]
        raw_data = data_all[data_start:data_end]
        counts_per_cell = np.diff(indptr_all[cell_start : cell_end + 1])
        raw_cell = np.repeat(
            np.arange(cell_end - cell_start, dtype=np.int32), counts_per_cell
        )
        mapped_gene = lookup[raw_gene]
        keep = mapped_gene >= 0
        block = sparse.coo_matrix(
            (raw_data[keep], (raw_cell[keep], mapped_gene[keep])),
            shape=(cell_end - cell_start, len(hvg_indices)),
            dtype=np.float32,
        ).tocsr()
        blocks.append(block)
        if (len(blocks) % 50) == 0:
            print(
                f"  retained HVGs for {cell_end:,}/{n_cells:,} cells",
                flush=True,
            )
    return sparse.vstack(blocks, format="csr")


def ensure_raw_csc_cache(
    matrix: h5py.Group,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decompress the unusually large HDF5 chunks once into reusable memmaps."""
    CSC_CACHE.mkdir(parents=True, exist_ok=True)
    data_path = CSC_CACHE / "data_float32.npy"
    gene_path = CSC_CACHE / "gene_indices_int32.npy"
    indptr_path = CSC_CACHE / "cell_indptr_int64.npy"
    nnz = int(matrix["data"].shape[0])
    n_cells = int(matrix["shape"][1])

    valid = False
    if data_path.exists() and gene_path.exists() and indptr_path.exists():
        try:
            cached_data = np.load(data_path, mmap_mode="r")
            cached_gene = np.load(gene_path, mmap_mode="r")
            cached_indptr = np.load(indptr_path, mmap_mode="r")
            valid = (
                cached_data.shape == (nnz,)
                and cached_data.dtype == np.float32
                and cached_gene.shape == (nnz,)
                and cached_gene.dtype == np.int32
                and cached_indptr.shape == (n_cells + 1,)
                and cached_indptr.dtype == np.int64
                and int(cached_indptr[-1]) == nnz
            )
        except (OSError, ValueError):
            valid = False
    if valid:
        print("Using the reusable decompressed CSC cache", flush=True)
        return cached_data, cached_gene, cached_indptr

    print("Creating one-time decompressed CSC cache", flush=True)
    data_cache = np.lib.format.open_memmap(
        data_path, mode="w+", dtype=np.float32, shape=(nnz,)
    )
    matrix["data"].read_direct(data_cache)
    data_cache.flush()
    del data_cache

    gene_cache = np.lib.format.open_memmap(
        gene_path, mode="w+", dtype=np.int32, shape=(nnz,)
    )
    matrix["indices"].read_direct(gene_cache)
    gene_cache.flush()
    del gene_cache

    indptr_cache = np.lib.format.open_memmap(
        indptr_path, mode="w+", dtype=np.int64, shape=(n_cells + 1,)
    )
    matrix["indptr"].read_direct(indptr_cache)
    indptr_cache.flush()
    del indptr_cache

    return (
        np.load(data_path, mmap_mode="r"),
        np.load(gene_path, mmap_mode="r"),
        np.load(indptr_path, mmap_mode="r"),
    )


def log_normalize(
    counts: sparse.csr_matrix, library_size: np.ndarray
) -> sparse.csr_matrix:
    normalized = counts.astype(np.float32, copy=True)
    normalized = normalized.multiply(
        (10_000 / np.maximum(library_size.astype(float), 1))[:, None]
    ).tocsr()
    np.log1p(normalized.data, out=normalized.data)
    return normalized


def build_hvg_statistics(
    mean: np.ndarray,
    variance: np.ndarray,
    detection: np.ndarray,
    features: pd.DataFrame,
) -> pd.DataFrame:
    dispersion = np.zeros_like(mean)
    eligible = (detection >= 0.005) & (mean > 0) & (variance > 0)
    dispersion[eligible] = variance[eligible] / np.maximum(mean[eligible], 1e-8)

    symbol = features["gene_symbol"].fillna("").astype(str)
    technical = symbol.str.startswith("MT-") | symbol.str.match(r"^RP[SL][0-9]")
    eligible &= ~technical.to_numpy()

    standardized = np.full(len(mean), -np.inf)
    eligible_index = np.flatnonzero(eligible)
    if len(eligible_index):
        bins = pd.qcut(
            mean[eligible_index],
            q=min(20, len(eligible_index)),
            labels=False,
            duplicates="drop",
        )
        bins = np.asarray(bins)
        for value in np.unique(bins):
            idx = eligible_index[bins == value]
            local = dispersion[idx]
            sd = local.std(ddof=1)
            standardized[idx] = (
                (local - local.mean()) / sd if sd > 0 else local - local.mean()
            )

    table = features.copy()
    table["hvg_estimation_mean_log_normalized"] = mean
    table["hvg_estimation_variance_log_normalized"] = variance
    table["hvg_estimation_detection_rate"] = detection
    table["hvg_estimation_dispersion"] = dispersion
    table["hvg_estimation_standardized_dispersion"] = standardized
    table["technical_gene_excluded"] = technical
    table["eligible_for_hvg"] = eligible
    order = np.argsort(standardized)[::-1]
    selected = order[np.isfinite(standardized[order])][:N_HVG]
    table["selected_hvg"] = False
    table.loc[selected, "selected_hvg"] = True
    table["hvg_rank"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    table.loc[selected, "hvg_rank"] = np.arange(1, len(selected) + 1)
    return table


OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)

metadata = pd.read_csv(META)
publisher = pd.read_csv(PUBLISHER)
if not np.array_equal(
    metadata["cellID"].to_numpy(dtype=object),
    publisher["sampleID"].to_numpy(dtype=object),
):
    raise RuntimeError("GSE178341 metadata and publisher annotations are not aligned")
metadata = metadata.copy()
metadata["patient_tissue_id"] = (
    metadata["PID"].astype(str) + "__" + metadata["SPECIMEN_TYPE"].astype(str)
)
library_cache = np.load(LIB_CACHE)
library_size = library_cache["library_size"].astype(np.float64)
if len(library_size) != len(metadata):
    raise RuntimeError("Cached full-cell library sizes do not match GSE178341")

with h5py.File(H5, "r") as h5:
    matrix = h5["matrix"]
    n_genes, n_cells = (int(x) for x in matrix["shape"][:])
    barcodes = decode(matrix["barcodes"][:])
    if n_cells != len(metadata) or not np.array_equal(
        barcodes, metadata["cellID"].to_numpy(dtype=object)
    ):
        raise RuntimeError("HDF5 cell order does not match metadata")
    features = pd.DataFrame(
        {
            "feature_index": np.arange(n_genes, dtype=int),
            "ensembl_id": decode(matrix["features/id"][:]),
            "gene_symbol": decode(matrix["features/name"][:]),
        }
    )
    data_all, gene_all, indptr_all = ensure_raw_csc_cache(matrix)

hvg_estimation_cells = sample_by_group(
    metadata["patient_tissue_id"],
    HVG_ESTIMATION_CAP_PER_PATIENT_TISSUE,
    rng,
)
print(
    f"Extracting {len(hvg_estimation_cells):,} label-free HVG-estimation cells",
    flush=True,
)
hvg_mean, hvg_variance, hvg_detection = estimate_hvg_moments_streaming(
    data_all,
    gene_all,
    indptr_all,
    hvg_estimation_cells,
    library_size,
    n_genes,
    n_cells,
)
hvg_table = build_hvg_statistics(
    hvg_mean,
    hvg_variance,
    hvg_detection,
    features,
)
hvg_indices = (
    hvg_table.loc[hvg_table["selected_hvg"]]
    .sort_values("hvg_rank")["feature_index"]
    .to_numpy(dtype=np.int64)
)
if len(hvg_indices) != N_HVG:
    raise RuntimeError(f"Expected {N_HVG} HVGs; obtained {len(hvg_indices)}")
print("Extracting the selected HVGs for all cells", flush=True)
full_hvg_counts = extract_hvg_all_cells(
    data_all,
    gene_all,
    indptr_all,
    hvg_indices,
    n_genes,
    n_cells,
)

full_hvg_normalized = log_normalize(full_hvg_counts, library_size)
sparse.save_npz(
    OUT / "full_cell_hvg_log_normalized_cells_by_genes.npz",
    full_hvg_normalized,
    compressed=True,
)
hvg_table.to_csv(
    OUT / "label_free_hvg_statistics.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
features.to_csv(
    OUT / "features.tsv.gz", sep="\t", index=False, compression="gzip"
)

print("Calculating global randomized-SVD leverage scores", flush=True)
svd = TruncatedSVD(
    n_components=N_LEVERAGE_DIMS,
    algorithm="randomized",
    n_iter=7,
    random_state=SEED,
)
svd_scores = svd.fit_transform(full_hvg_normalized).astype(np.float32)
singular = np.maximum(svd.singular_values_.astype(np.float64), 1e-12)
left_vectors = svd_scores.astype(np.float64) / singular[None, :]
leverage = np.sum(left_vectors**2, axis=1)
leverage = np.maximum(leverage, np.finfo(float).eps)
del left_vectors

sketch_cells = sample_by_group(
    metadata["patient_tissue_id"],
    SKETCH_CAP_PER_PATIENT_TISSUE,
    rng,
    probabilities=leverage,
)
sketch_meta = metadata.iloc[sketch_cells].copy().reset_index(drop=True)
sketch_meta.insert(0, "full_matrix_column", sketch_cells)
# Publisher labels are appended only after label-free sampling for validation.
sketch_meta["publisher_top_level"] = publisher.iloc[sketch_cells][
    "clTopLevel"
].to_numpy(dtype=object)
sketch_meta["publisher_mid_level"] = publisher.iloc[sketch_cells][
    "clMidwayPr"
].to_numpy(dtype=object)
sketch_meta["publisher_subcluster"] = publisher.iloc[sketch_cells][
    "cl295v11SubFull"
].to_numpy(dtype=object)
sketch_meta["global_leverage_score"] = leverage[sketch_cells]
sketch_meta.to_csv(
    OUT / "label_free_sketch_metadata.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
pd.DataFrame(
    {
        "full_matrix_column": np.arange(n_cells, dtype=int),
        "cellID": metadata["cellID"],
        "patient_tissue_id": metadata["patient_tissue_id"],
        "global_leverage_score": leverage,
        "selected_for_sketch": np.isin(np.arange(n_cells), sketch_cells),
    }
).to_csv(
    OUT / "full_cell_leverage_scores.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

print(f"Running PCA and UMAP on {len(sketch_cells):,} sketch cells", flush=True)
sketch_normalized = full_hvg_normalized[sketch_cells].tocsr()
sketch_dense = sketch_normalized.toarray().astype(np.float32, copy=False)
gene_mean = sketch_dense.mean(axis=0, dtype=np.float64).astype(np.float32)
gene_sd = sketch_dense.std(axis=0, ddof=1, dtype=np.float64).astype(np.float32)
gene_sd = np.where(gene_sd > 0, gene_sd, 1).astype(np.float32)
sketch_dense -= gene_mean
sketch_dense /= gene_sd
np.clip(sketch_dense, -10, 10, out=sketch_dense)

pca = PCA(
    n_components=N_PCS,
    svd_solver="randomized",
    random_state=SEED,
)
pcs = pca.fit_transform(sketch_dense).astype(np.float32)
del sketch_dense

umap_model = umap.UMAP(
    n_neighbors=30,
    min_dist=0.30,
    n_components=2,
    metric="cosine",
    random_state=SEED,
    n_jobs=1,
    low_memory=True,
)
embedding = umap_model.fit_transform(pcs[:, :RETAINED_PCS]).astype(np.float32)
graph = sparse.triu(umap_model.graph_, k=1).tocoo()

np.savez_compressed(
    OUT / "label_free_sketch_model_arrays.npz",
    sketch_cells=sketch_cells,
    hvg_feature_indices=hvg_indices,
    gene_mean=gene_mean,
    gene_sd=gene_sd,
    pca_components=pca.components_.astype(np.float32),
    pca_mean=pca.mean_.astype(np.float32),
    pca_explained_variance=pca.explained_variance_.astype(np.float64),
    pca_explained_variance_ratio=pca.explained_variance_ratio_.astype(np.float64),
    leverage_singular_values=svd.singular_values_.astype(np.float64),
)
np.save(OUT / "label_free_sketch_pca.npy", pcs)
np.save(OUT / "label_free_sketch_umap.npy", embedding)
sparse.save_npz(
    OUT / "label_free_sketch_log_normalized_cells_by_hvgs.npz",
    sketch_normalized,
    compressed=True,
)
pd.DataFrame(
    {
        "source": graph.row.astype(np.int32),
        "target": graph.col.astype(np.int32),
        "weight": graph.data.astype(np.float32),
    }
).to_csv(
    OUT / "label_free_sketch_umap_graph_edges.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

patient_tissue_counts = (
    metadata.groupby(["PID", "SPECIMEN_TYPE", "patient_tissue_id"], observed=True)
    .size()
    .rename("full_cells")
    .reset_index()
)
sketch_counts = (
    sketch_meta.groupby(
        ["PID", "SPECIMEN_TYPE", "patient_tissue_id"], observed=True
    )
    .size()
    .rename("sketch_cells")
    .reset_index()
)
patient_tissue_counts.merge(
    sketch_counts,
    on=["PID", "SPECIMEN_TYPE", "patient_tissue_id"],
    how="left",
).to_csv(
    OUT / "label_free_sketch_patient_tissue_audit.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame(
    {
        "field": [
            "full_cells",
            "patient_tissue_strata",
            "hvg_estimation_cap_per_patient_tissue",
            "hvg_estimation_cells",
            "sketch_cap_per_patient_tissue",
            "sketch_cells",
            "HVGs",
            "leverage_dimensions",
            "PCA_components",
            "PCA_components_for_neighbors",
            "publisher_labels_used_before_clustering",
            "random_seed",
            "python",
            "platform",
        ],
        "value": [
            n_cells,
            metadata["patient_tissue_id"].nunique(),
            HVG_ESTIMATION_CAP_PER_PATIENT_TISSUE,
            len(hvg_estimation_cells),
            SKETCH_CAP_PER_PATIENT_TISSUE,
            len(sketch_cells),
            N_HVG,
            N_LEVERAGE_DIMS,
            N_PCS,
            RETAINED_PCS,
            False,
            SEED,
            sys.version.replace("\n", " "),
            platform.platform(),
        ],
    }
).to_csv(OUT / "label_free_sketch_analysis_metadata.tsv", sep="\t", index=False)

print(
    f"Completed label-free sketch: {len(sketch_cells):,} cells from "
    f"{metadata['patient_tissue_id'].nunique()} patient-tissue strata",
    flush=True,
)
