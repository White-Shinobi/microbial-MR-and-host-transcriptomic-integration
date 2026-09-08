#!/usr/bin/env python3
"""Freeze V2 single-cell outputs as submission supplementary tables S24-S34."""

from __future__ import annotations

import csv
import gzip
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/omics/crc_single_cell/GSE132465"
RESULT = ROOT / "results/single_cell/GSE132465"
SOURCE = ROOT / "manuscript/source_data"
OUT = ROOT / "manuscript/tables/supplementary"


def copy(source: Path, name: str) -> None:
    destination = OUT / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cohort_audit() -> None:
    audit_path = RESULT / "core/QC_cell_audit.tsv.gz"
    with gzip.open(audit_path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    retained = [row for row in rows if row["retained"].lower() == "true"]
    patients = {row["Patient"] for row in retained}
    tumor_patients = {row["Patient"] for row in retained if row["Class"] == "Tumor"}
    normal_patients = {row["Patient"] for row in retained if row["Class"] == "Normal"}
    tumor_samples = {row["Sample"] for row in retained if row["Class"] == "Tumor"}
    normal_samples = {row["Sample"] for row in retained if row["Class"] == "Normal"}
    matrix = RAW / "10x_sparse/matrix.mtx.gz"
    annotation = RAW / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
    record = {
        "accession": "GSE132465",
        "source": "GEO-processed raw UMI matrix and cell annotations",
        "processed_cells": len(rows),
        "qc_retained_cells": len(retained),
        "qc_retained_tumor_cells": sum(row["Class"] == "Tumor" for row in retained),
        "qc_retained_normal_cells": sum(row["Class"] == "Normal" for row in retained),
        "patients": len(patients),
        "tumor_samples": len(tumor_samples),
        "normal_samples": len(normal_samples),
        "matched_patients": len(tumor_patients & normal_patients),
        "matrix_path": str(matrix),
        "matrix_sha256": sha256(matrix),
        "annotation_path": str(annotation),
        "annotation_sha256": sha256(annotation),
        "reuse_mode": "V1 raw-data symbolic link; no repeat download",
    }
    destination = OUT / "Table_S24_single_cell_cohort_audit.tsv"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record), delimiter="\t")
        writer.writeheader()
        writer.writerow(record)


def biomarker_lr_table() -> None:
    source = (
        RESULT
        / "cellchat_nboot1000_fdr/"
        "stromal_biomarker_ligand_receptor_correlations.tsv"
    )
    destination = OUT / "Table_S34_stromal_biomarker_ligand_receptor_correlations.tsv"
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = [
        "biomarker",
        "ligand_receptor_gene",
        "n_patient_tissue_profiles",
        "spearman_rho",
        "spearman_p",
        "spearman_fdr",
        "bh_fdr_supported",
        "reference_display_rule_supported",
    ]
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **{key: row[key] for key in fields[:6]},
                    "bh_fdr_supported": float(row["spearman_fdr"]) < 0.05,
                    "reference_display_rule_supported": row["ref1_supported"],
                }
            )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cohort_audit()
    copy(RESULT / "core/QC_summary.tsv", "Table_S25_single_cell_QC_summary.tsv")
    copy(
        RESULT / "core/biomarker_tumor_normal_by_cell_type.tsv",
        "Table_S26_single_cell_biomarker_celltype_comparisons.tsv",
    )
    copy(
        RESULT / "core/key_cell_type_ranking.tsv",
        "Table_S27_single_cell_key_cell_ranking.tsv",
    )
    copy(
        RESULT / "core/cell_composition_tumor_normal.tsv",
        "Table_S28_single_cell_composition.tsv",
    )
    copy(
        RESULT / "reactomegsa/ReactomeGSA_cell_type_scores.tsv.gz",
        "Table_S29_ReactomeGSA_cell_type_scores.tsv.gz",
    )
    copy(
        RESULT / "stromal_pseudotime/biomarker_pseudotime_spearman.tsv",
        "Table_S30_stromal_biomarker_pseudotime.tsv",
    )
    if not (
        OUT / "Table_S31_stromal_pseudotime_tissue_comparison.tsv"
    ).exists():
        raise FileNotFoundError("Pseudotime tissue-comparison table S31 is missing")
    copy(
        SOURCE / "Figure_7C_stromal_subcluster_statistics.tsv",
        "Table_S32_stromal_subcluster_composition.tsv",
    )
    copy(
        RESULT
        / "cellchat_nboot1000_fdr/"
        "cellchat_significant_interactions_by_tissue.tsv",
        "Table_S33_CellChat_interactions_nboot1000_BHFDR.tsv",
    )
    biomarker_lr_table()
    print("Frozen supplementary tables S24-S34")


if __name__ == "__main__":
    main()
