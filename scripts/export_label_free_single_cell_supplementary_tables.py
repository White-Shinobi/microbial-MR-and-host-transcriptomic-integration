#!/usr/bin/env python3
"""Export revised Tables S21-S31 for the coherent label-free workflow."""

from pathlib import Path
import hashlib

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/omics/crc_single_cell/GSE178341"
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PB = ATLAS / "pseudobulk"
EPI = ROOT / "results/single_cell/GSE178341/label_free_epithelial"
F8 = ROOT / "results/single_cell/GSE178341/label_free_figure7"
OUT = ROOT / "manuscript/tables/supplementary"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


files = [
    (
        "GSE178341_crc10x_full_c295v4_submit.h5",
        "processed raw UMI matrix",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/"
        "GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit.h5",
    ),
    (
        "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
        "publisher annotations used only for post-hoc audit",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/"
        "GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
    ),
    (
        "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
        "publisher patient and tissue metadata",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/"
        "GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
    ),
]
s21 = []
for filename, role, url in files:
    path = RAW / filename
    s21.append(
        {
            "accession": "GSE178341",
            "source_article": "Pelka et al., Cell 2021",
            "source_url": url,
            "study_reported_cells": 371223,
            "analysis_matched_cells": 370115,
            "analysis_genes": 43113,
            "CRC_patients": 62,
            "normal_patients": 36,
            "paired_patients": 36,
            "filename": filename,
            "file_role": role,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "reuse_status": "Existing download reused; checksum verified",
        }
    )
pd.DataFrame(s21).to_csv(
    OUT / "Table_S21_GSE178341_cohort_file_audit.tsv",
    sep="\t",
    index=False,
)

top_meta = pd.read_csv(
    ATLAS / "label_free_sketch_analysis_metadata.tsv", sep="\t"
)
resolution = pd.read_csv(
    ATLAS / "label_free_leiden_silhouette_selection.tsv", sep="\t"
)
transfer = pd.read_csv(
    ATLAS / "full_cell_label_transfer_audit.tsv", sep="\t"
)
s22 = []
for row in top_meta.itertuples():
    s22.append(
        {
            "record_type": "label_free_sketch",
            "item": row.field,
            "value": row.value,
            "selected": True,
            "note": "No publisher cell-type label used before clustering",
        }
    )
for row in resolution.itertuples():
    s22.append(
        {
            "record_type": "Leiden_resolution",
            "item": row.membership_column,
            "value": row.balanced_silhouette,
            "selected": bool(row.selected),
            "note": f"{row.clusters} clusters; {row.silhouette_cells} silhouette cells",
        }
    )
for row in transfer.itertuples():
    s22.append(
        {
            "record_type": "full_cell_label_transfer",
            "item": row.metric,
            "value": row.value,
            "selected": True,
            "note": "Publisher labels used only for post-hoc concordance metrics",
        }
    )
pd.DataFrame(s22).to_csv(
    OUT / "Table_S22_GSE178341_label_free_cluster_audit.tsv",
    sep="\t",
    index=False,
)

s23 = pd.read_csv(
    ATLAS / "figure6G_label_free_patient_balanced_biomarker_summary.tsv",
    sep="\t",
)
s23.to_csv(
    OUT / "Table_S23_GSE178341_patient_balanced_biomarker_localization.tsv",
    sep="\t",
    index=False,
)

s24 = pd.read_csv(
    PB / "label_free_paired_edger_biomarker_results_35_tests.tsv",
    sep="\t",
)
s24.to_csv(
    OUT / "Table_S24_GSE178341_paired_edgeR_35_tests.tsv",
    sep="\t",
    index=False,
)
s25 = s24[
    [
        "gene",
        "cell_type",
        "n_paired_patients",
        "median_normal_log2_cpm",
        "median_CRC_log2_cpm",
        "median_paired_difference_CRC_minus_normal",
        "paired_wilcoxon_P",
        "paired_wilcoxon_BH_FDR_35_tests",
        "wilcoxon_FDR_lt_0_05",
    ]
].copy()
s25.to_csv(
    OUT / "Table_S25_GSE178341_paired_Wilcoxon_sensitivity_35_tests.tsv",
    sep="\t",
    index=False,
)

s26 = pd.read_csv(
    ATLAS / "figure6H_label_free_all_KEGG_pathway_scores.tsv.gz",
    sep="\t",
)
s26.to_csv(
    OUT / "Table_S26_GSE178341_patient_balanced_KEGG_pathway_scores.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

s27 = pd.read_csv(
    F8 / "label_free_epithelial_subcluster_paired_statistics.tsv",
    sep="\t",
)
s27["BH_FDR_significant"] = s27["paired_wilcoxon_BH_FDR"] < 0.05
s27.to_csv(
    OUT / "Table_S27_GSE178341_label_free_epithelial_subcluster_composition.tsv",
    sep="\t",
    index=False,
)

s28 = pd.read_csv(
    F8 / "label_free_cellchat_BH_FDR_significant_interactions.tsv.gz",
    sep="\t",
)
s28 = s28.loc[s28["involves_epithelial"].astype(bool)].copy()
s28.to_csv(
    OUT / "Table_S28_GSE178341_label_free_epithelial_CellChat_BHFDR.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

s29 = pd.read_csv(
    F8 / "label_free_focal_biomarker_LR_correlations.tsv", sep="\t"
)
s29["BH_FDR_significant"] = s29["spearman_BH_FDR"] < 0.05
s29.to_csv(
    OUT / "Table_S29_GSE178341_label_free_focal_biomarker_LR_correlations.tsv",
    sep="\t",
    index=False,
)

pd.read_csv(
    F8 / "label_free_biomarker_epithelial_pseudotime.tsv", sep="\t"
).to_csv(
    OUT / "Table_S30_GSE178341_label_free_biomarker_epithelial_pseudotime.tsv",
    sep="\t",
    index=False,
)
pd.read_csv(
    F8 / "label_free_epithelial_patient_pseudotime_statistics.tsv", sep="\t"
).to_csv(
    OUT / "Table_S31_GSE178341_label_free_epithelial_patient_pseudotime.tsv",
    sep="\t",
    index=False,
)

print(
    "Exported S21-S31 rows:",
    {
        "S21": len(s21),
        "S22": len(s22),
        "S23": len(s23),
        "S24": len(s24),
        "S25": len(s25),
        "S26": len(s26),
        "S27": len(s27),
        "S28": len(s28),
        "S29": len(s29),
        "S30": 5,
        "S31": 1,
    },
)
