#!/usr/bin/env python3
"""Aggregate GSE178341 using the newly mapped label-free broad cell types."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
OUT = ATLAS / "pseudobulk"
CACHE = ATLAS / "raw_csc_cache"
TARGET_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)
BIOMARKERS = ["GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"]
CELL_TYPES = [
    "Epithelial", "T/NK/ILC", "Myeloid", "B", "Plasma", "Stromal", "Mast"
]
MIN_CELLS = 20


def safe_wilcoxon(tumor: np.ndarray, normal: np.ndarray) -> float:
    difference = tumor - normal
    if len(difference) < 3:
        return np.nan
    if np.allclose(difference, 0):
        return 1.0
    return float(
        wilcoxon(
            difference,
            zero_method="wilcox",
            alternative="two-sided",
            method="auto",
        ).pvalue
    )


OUT.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(
    ATLAS / "full_cell_label_free_mapped_metadata.tsv.gz", sep="\t"
)
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
target = np.load(TARGET_CACHE)
library_size = target["library_size"].astype(np.int64)
target_counts = target["target_counts"].astype(np.int32)
if target_counts.shape != (len(BIOMARKERS), len(metadata)):
    raise RuntimeError("Five-gene count cache and mapped metadata differ")

metadata["cell_type"] = metadata["label_free_marker_cell_type"]
metadata["pseudobulk_id"] = (
    metadata["PID"].astype(str)
    + "__"
    + metadata["SPECIMEN_TYPE"].astype(str)
    + "__"
    + metadata["cell_type"].astype(str)
)
group_meta = (
    metadata.groupby(
        ["pseudobulk_id", "PID", "SPECIMEN_TYPE", "cell_type"],
        observed=True,
        sort=True,
    )
    .size()
    .rename("n_cells")
    .reset_index()
)
group_meta.insert(0, "pseudobulk_index", np.arange(len(group_meta), dtype=int))
lookup = dict(zip(group_meta["pseudobulk_id"], group_meta["pseudobulk_index"]))
cell_group = metadata["pseudobulk_id"].map(lookup).to_numpy(dtype=np.int32)
group_meta.to_csv(
    OUT / "label_free_pseudobulk_sample_metadata.tsv",
    sep="\t",
    index=False,
)

# Patient-balanced biomarker localization: first summarize each patient,
# tissue and mapped broad cell type, then average those patient-level values.
patient_rows = []
base = metadata[["PID", "SPECIMEN_TYPE", "cell_type"]].copy()
for gene_index, gene in enumerate(BIOMARKERS):
    counts = target_counts[gene_index].astype(np.int64, copy=False)
    log_expression = np.log1p(
        10_000 * counts / np.maximum(library_size, 1)
    )
    work = base.copy()
    work["detected"] = counts > 0
    work["log_expression"] = log_expression
    work["gene_umi"] = counts
    work["total_umi"] = library_size
    grouped = (
        work.groupby(
            ["PID", "SPECIMEN_TYPE", "cell_type"],
            observed=True,
            sort=False,
        )
        .agg(
            n_cells=("detected", "size"),
            detection_rate=("detected", "mean"),
            mean_log_normalized_expression=("log_expression", "mean"),
            gene_umi=("gene_umi", "sum"),
            total_umi=("total_umi", "sum"),
        )
        .reset_index()
    )
    grouped.insert(0, "gene", gene)
    grouped["patient_gene_detected"] = grouped["gene_umi"] > 0
    grouped["eligible_min_20_cells"] = grouped["n_cells"] >= MIN_CELLS
    grouped["log2_cpm"] = np.log2(
        grouped["gene_umi"] / np.maximum(grouped["total_umi"], 1) * 1e6
        + 0.5
    )
    patient_rows.append(grouped)

patient = pd.concat(patient_rows, ignore_index=True)
patient.to_csv(
    OUT / "label_free_patient_tissue_cell_type_biomarker_values.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
eligible = patient.loc[patient["eligible_min_20_cells"]].copy()
localization = (
    eligible.groupby(
        ["gene", "cell_type", "SPECIMEN_TYPE"],
        observed=True,
        sort=False,
    )
    .agg(
        eligible_patients=("PID", "nunique"),
        patient_balanced_detection_rate=("detection_rate", "mean"),
        patient_balanced_mean_expression=(
            "mean_log_normalized_expression", "mean"
        ),
        patient_detection_coverage=("patient_gene_detected", "mean"),
        median_patient_log2_cpm=("log2_cpm", "median"),
    )
    .reset_index()
)
localization.to_csv(
    OUT / "label_free_patient_balanced_biomarker_localization.tsv",
    sep="\t",
    index=False,
)

paired_rows = []
for (gene, cell_type), frame in eligible.groupby(
    ["gene", "cell_type"], observed=True, sort=False
):
    wide = frame.pivot_table(
        index="PID", columns="SPECIMEN_TYPE", values="log2_cpm", aggfunc="first"
    )
    if {"N", "T"}.issubset(wide.columns):
        wide = wide.dropna(subset=["N", "T"])
        normal = wide["N"].to_numpy(dtype=float)
        tumor = wide["T"].to_numpy(dtype=float)
    else:
        normal = np.array([], dtype=float)
        tumor = np.array([], dtype=float)
    paired_rows.append(
        {
            "gene": gene,
            "cell_type": cell_type,
            "n_paired_patients": len(normal),
            "median_normal_log2_cpm": (
                np.median(normal) if len(normal) else np.nan
            ),
            "median_CRC_log2_cpm": (
                np.median(tumor) if len(tumor) else np.nan
            ),
            "median_paired_difference_CRC_minus_normal": (
                np.median(tumor - normal) if len(normal) else np.nan
            ),
            "paired_wilcoxon_P": safe_wilcoxon(tumor, normal),
        }
    )
paired = pd.DataFrame(paired_rows)
valid = paired["paired_wilcoxon_P"].notna()
paired["paired_wilcoxon_BH_FDR_35_tests"] = np.nan
paired.loc[valid, "paired_wilcoxon_BH_FDR_35_tests"] = multipletests(
    paired.loc[valid, "paired_wilcoxon_P"], method="fdr_bh"
)[1]
paired.to_csv(
    OUT / "label_free_paired_wilcoxon_biomarker_results_35_tests.tsv",
    sep="\t",
    index=False,
)

print("Aggregating the full raw matrix to label-free pseudobulk samples")
data = np.load(CACHE / "data_float32.npy", mmap_mode="r")
gene_index = np.load(CACHE / "gene_indices_int32.npy", mmap_mode="r")
indptr = np.load(CACHE / "cell_indptr_int64.npy", mmap_mode="r")
full = sparse.csc_matrix(
    (data, gene_index, indptr),
    shape=(len(features), len(metadata)),
    copy=False,
)
cell_to_group = sparse.csr_matrix(
    (
        np.ones(len(metadata), dtype=np.float32),
        (np.arange(len(metadata), dtype=np.int64), cell_group),
    ),
    shape=(len(metadata), len(group_meta)),
)
pseudobulk = (full @ cell_to_group).tocsc()
if pseudobulk.nnz and pseudobulk.data.max() >= np.iinfo(np.int32).max:
    raise RuntimeError("Pseudobulk count exceeded int32 range")
pseudobulk.data = np.rint(pseudobulk.data).astype(np.int32, copy=False)
sparse.save_npz(
    OUT / "label_free_pseudobulk_raw_counts_genes_by_samples.npz",
    pseudobulk,
    compressed=True,
)

pd.DataFrame(
    {
        "metric": [
            "full_cells",
            "genes",
            "mapped_broad_cell_types",
            "pseudobulk_samples",
            "minimum_cells_for_localization_and_pairing",
            "formal_biomarker_tests",
        ],
        "value": [
            len(metadata),
            len(features),
            metadata["cell_type"].nunique(),
            len(group_meta),
            MIN_CELLS,
            len(paired),
        ],
    }
).to_csv(
    OUT / "label_free_pseudobulk_preparation_audit.tsv",
    sep="\t",
    index=False,
)
print(
    f"Completed {pseudobulk.shape[1]} label-free pseudobulk samples and "
    f"{len(paired)} biomarker-cell-type tests"
)
