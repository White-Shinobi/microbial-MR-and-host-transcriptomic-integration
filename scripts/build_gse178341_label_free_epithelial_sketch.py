#!/usr/bin/env python3
"""Build an epithelial-specific label-free sketch from newly mapped cells."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
import umap


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
OUT = (
    ROOT
    / "results/single_cell/GSE178341/label_free_epithelial"
)
CACHE = ATLAS / "raw_csc_cache"
TARGET_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)
SEED = 20260723
HVG_CAP = 200
SKETCH_CAP = 500
N_HVG = 2000
BLOCK = 5000


def sample_by_group(
    groups: pd.Series,
    cap: int,
    rng: np.random.Generator,
    probability: np.ndarray | None = None,
) -> np.ndarray:
    chosen = []
    for indices in groups.groupby(groups, sort=True).indices.values():
        indices = np.asarray(indices, dtype=int)
        if len(indices) > cap:
            p = None
            if probability is not None:
                p = np.maximum(probability[indices], 0)
                p = p / p.sum() if np.isfinite(p).all() and p.sum() > 0 else None
            indices = rng.choice(indices, cap, replace=False, p=p)
        chosen.append(np.sort(indices))
    return np.sort(np.concatenate(chosen))


def hvg_moments(
    selected_global: np.ndarray,
    data: np.ndarray,
    genes: np.ndarray,
    indptr: np.ndarray,
    library: np.ndarray,
    n_genes: int,
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = np.zeros(n_cells, dtype=bool)
    selected[selected_global] = True
    total = len(selected_global)
    sum_x = np.zeros(n_genes, dtype=float)
    sumsq_x = np.zeros(n_genes, dtype=float)
    detected = np.zeros(n_genes, dtype=np.int64)
    for start in range(0, n_cells, BLOCK):
        end = min(start + BLOCK, n_cells)
        keep_cells = selected[start:end]
        if not keep_cells.any():
            continue
        left, right = int(indptr[start]), int(indptr[end])
        local = np.repeat(
            np.arange(end - start, dtype=np.int32),
            np.diff(indptr[start : end + 1]),
        )
        keep = keep_cells[local]
        global_cell = local[keep] + start
        gene = genes[left:right][keep]
        value = np.log1p(
            data[left:right][keep]
            * (10_000 / np.maximum(library[global_cell], 1))
        )
        sum_x += np.bincount(gene, weights=value, minlength=n_genes)
        sumsq_x += np.bincount(
            gene, weights=value * value, minlength=n_genes
        )
        detected += np.bincount(gene, minlength=n_genes)
    mean = sum_x / total
    variance = np.maximum(sumsq_x / total - mean**2, 0)
    return mean, variance, detected / total


def select_hvgs(
    features: pd.DataFrame,
    mean: np.ndarray,
    variance: np.ndarray,
    detection: np.ndarray,
) -> pd.DataFrame:
    dispersion = np.zeros_like(mean)
    eligible = (detection >= 0.005) & (mean > 0) & (variance > 0)
    symbols = features["gene_symbol"].fillna("").astype(str)
    technical = symbols.str.startswith("MT-") | symbols.str.match(
        r"^RP[SL][0-9]"
    )
    eligible &= ~technical.to_numpy()
    dispersion[eligible] = variance[eligible] / np.maximum(
        mean[eligible], 1e-8
    )
    standardized = np.full(len(mean), -np.inf)
    idx = np.flatnonzero(eligible)
    bins = pd.qcut(
        mean[idx], q=min(20, len(idx)), labels=False, duplicates="drop"
    )
    for value in np.unique(bins):
        local_idx = idx[np.asarray(bins) == value]
        local = dispersion[local_idx]
        sd = local.std(ddof=1)
        standardized[local_idx] = (
            (local - local.mean()) / sd if sd > 0 else local - local.mean()
        )
    order = np.argsort(standardized)[::-1]
    selected = order[np.isfinite(standardized[order])][:N_HVG]
    table = features.copy()
    table["mean_log_normalized"] = mean
    table["variance_log_normalized"] = variance
    table["detection_rate"] = detection
    table["dispersion"] = dispersion
    table["standardized_dispersion"] = standardized
    table["technical_excluded"] = technical
    table["selected_hvg"] = False
    table.loc[selected, "selected_hvg"] = True
    table["hvg_rank"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    table.loc[selected, "hvg_rank"] = np.arange(1, len(selected) + 1)
    return table


def extract_epi_hvgs(
    epi_global: np.ndarray,
    hvg_indices: np.ndarray,
    data: np.ndarray,
    genes: np.ndarray,
    indptr: np.ndarray,
    n_genes: int,
    n_cells: int,
) -> sparse.csr_matrix:
    hvg_lookup = np.full(n_genes, -1, dtype=np.int32)
    hvg_lookup[hvg_indices] = np.arange(len(hvg_indices), dtype=np.int32)
    epi_lookup = np.full(n_cells, -1, dtype=np.int32)
    epi_lookup[epi_global] = np.arange(len(epi_global), dtype=np.int32)
    blocks = []
    for start in range(0, n_cells, BLOCK):
        end = min(start + BLOCK, n_cells)
        epi_rows_global = epi_lookup[start:end]
        local_epi = np.flatnonzero(epi_rows_global >= 0)
        if not len(local_epi):
            continue
        local_to_block = np.full(end - start, -1, dtype=np.int32)
        local_to_block[local_epi] = np.arange(len(local_epi), dtype=np.int32)
        left, right = int(indptr[start]), int(indptr[end])
        local_cell = np.repeat(
            np.arange(end - start, dtype=np.int32),
            np.diff(indptr[start : end + 1]),
        )
        mapped_cell = local_to_block[local_cell]
        mapped_gene = hvg_lookup[genes[left:right]]
        keep = (mapped_cell >= 0) & (mapped_gene >= 0)
        block = sparse.coo_matrix(
            (
                data[left:right][keep],
                (mapped_cell[keep], mapped_gene[keep]),
            ),
            shape=(len(local_epi), len(hvg_indices)),
            dtype=np.float32,
        ).tocsr()
        blocks.append(block)
    return sparse.vstack(blocks, format="csr")


OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)
full_meta = pd.read_csv(
    ATLAS / "full_cell_label_free_mapped_metadata.tsv.gz", sep="\t"
)
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
epi_meta = full_meta.loc[
    full_meta["label_free_marker_cell_type"].eq("Epithelial")
].copy()
epi_meta = epi_meta.sort_values("full_matrix_column").reset_index(drop=True)
epi_meta["patient_tissue_id"] = (
    epi_meta["PID"].astype(str)
    + "__"
    + epi_meta["SPECIMEN_TYPE"].astype(str)
)
epi_global = epi_meta["full_matrix_column"].to_numpy(dtype=int)
n_cells = len(full_meta)
n_genes = len(features)
data = np.load(CACHE / "data_float32.npy", mmap_mode="r")
genes = np.load(CACHE / "gene_indices_int32.npy", mmap_mode="r")
indptr = np.load(CACHE / "cell_indptr_int64.npy", mmap_mode="r")
library = np.load(TARGET_CACHE)["library_size"].astype(float)

hvg_local = sample_by_group(
    epi_meta["patient_tissue_id"], HVG_CAP, rng
)
hvg_global = epi_global[hvg_local]
mean, variance, detection = hvg_moments(
    hvg_global, data, genes, indptr, library, n_genes, n_cells
)
hvg_table = select_hvgs(features, mean, variance, detection)
hvg_indices = (
    hvg_table.loc[hvg_table["selected_hvg"]]
    .sort_values("hvg_rank")["feature_index"]
    .to_numpy(dtype=int)
)
epi_counts = extract_epi_hvgs(
    epi_global, hvg_indices, data, genes, indptr, n_genes, n_cells
)
epi_norm = epi_counts.astype(np.float32).multiply(
    (10_000 / np.maximum(library[epi_global], 1))[:, None]
).tocsr()
np.log1p(epi_norm.data, out=epi_norm.data)
sparse.save_npz(
    OUT / "full_epithelial_hvg_log_normalized.npz",
    epi_norm,
    compressed=True,
)

svd = TruncatedSVD(
    n_components=50, n_iter=7, random_state=SEED
)
svd_score = svd.fit_transform(epi_norm).astype(np.float32)
left = svd_score / np.maximum(svd.singular_values_, 1e-12)[None, :]
leverage = np.sum(left * left, axis=1)
sketch_local = sample_by_group(
    epi_meta["patient_tissue_id"], SKETCH_CAP, rng, leverage
)
sketch_meta = epi_meta.iloc[sketch_local].copy().reset_index(drop=True)
sketch_meta.insert(0, "epithelial_row", sketch_local)
sketch_meta["epithelial_leverage_score"] = leverage[sketch_local]

sketch_sparse = epi_norm[sketch_local].tocsr()
dense = sketch_sparse.toarray().astype(np.float32, copy=False)
gene_mean = dense.mean(axis=0, dtype=np.float64).astype(np.float32)
gene_sd = dense.std(axis=0, ddof=1, dtype=np.float64).astype(np.float32)
gene_sd = np.where(gene_sd > 0, gene_sd, 1).astype(np.float32)
dense -= gene_mean
dense /= gene_sd
np.clip(dense, -10, 10, out=dense)
pca = PCA(n_components=50, svd_solver="randomized", random_state=SEED)
pcs = pca.fit_transform(dense).astype(np.float32)
embedding_model = umap.UMAP(
    n_neighbors=30,
    min_dist=0.25,
    n_components=2,
    metric="cosine",
    random_state=SEED,
    n_jobs=1,
    low_memory=True,
)
embedding = embedding_model.fit_transform(pcs[:, :30]).astype(np.float32)
graph = sparse.triu(embedding_model.graph_, k=1).tocoo()

hvg_table.to_csv(
    OUT / "epithelial_hvg_statistics.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
sketch_meta.to_csv(
    OUT / "epithelial_label_free_sketch_metadata.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
sparse.save_npz(
    OUT / "epithelial_sketch_log_normalized.npz",
    sketch_sparse,
    compressed=True,
)
np.save(OUT / "epithelial_sketch_pca.npy", pcs)
np.save(OUT / "epithelial_sketch_umap.npy", embedding)
np.savez_compressed(
    OUT / "epithelial_sketch_model_arrays.npz",
    epi_global=epi_global,
    sketch_local=sketch_local,
    hvg_indices=hvg_indices,
    gene_mean=gene_mean,
    gene_sd=gene_sd,
    pca_components=pca.components_.astype(np.float32),
    pca_mean=pca.mean_.astype(np.float32),
    explained_variance=pca.explained_variance_,
    explained_variance_ratio=pca.explained_variance_ratio_,
)
pd.DataFrame(
    {
        "source": graph.row.astype(np.int32),
        "target": graph.col.astype(np.int32),
        "weight": graph.data.astype(np.float32),
    }
).to_csv(
    OUT / "epithelial_sketch_graph_edges.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
pd.DataFrame(
    {
        "metric": [
            "mapped_epithelial_cells",
            "patient_tissue_strata",
            "HVG_estimation_cells",
            "HVGs",
            "epithelial_sketch_cells",
            "publisher_subtype_used_before_clustering",
            "seed",
        ],
        "value": [
            len(epi_meta),
            epi_meta["patient_tissue_id"].nunique(),
            len(hvg_local),
            len(hvg_indices),
            len(sketch_local),
            False,
            SEED,
        ],
    }
).to_csv(
    OUT / "epithelial_sketch_audit.tsv", sep="\t", index=False
)
print(
    f"Completed epithelial label-free sketch: {len(sketch_local):,}/"
    f"{len(epi_meta):,} cells"
)
