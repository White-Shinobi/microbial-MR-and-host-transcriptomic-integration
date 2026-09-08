#!/usr/bin/env python3
"""Reuse revised-MR regional files for final-biomarker three-layer coloc."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GATE = ROOT / "results/bulk/crc/machine_learning/biomarker_gate_summary.tsv"
MAPPING = ROOT / "results/gmrg/ref1_candidate_feature_snp_gene.tsv"
REGIONS = ROOT / "data/processed/coloc/v2_fdr_mr_loci/region_manifest.tsv"
OUT = ROOT / "data/processed/coloc/v2_biomarker_loci"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


selected = {
    row["gene_symbol"] for row in rows(GATE)
    if row.get("gene_symbol") and truthy(row.get("final_biomarker", ""))
}
mapped = [
    row for row in rows(MAPPING)
    if row.get("gene_symbol") in selected and row.get("ensembl_gene_id")
]
regions = rows(REGIONS)
OUT.mkdir(parents=True, exist_ok=True)
for old in OUT.glob("*"):
    if old.is_file() or old.is_symlink():
        old.unlink()

manifest: list[dict[str, str]] = []
seen: set[tuple[str, str, str]] = set()
for row in mapped:
    key = (row["gene_symbol"], row["accession"], row["SNP"])
    if key in seen:
        continue
    seen.add(key)
    exposure = next(
        item for item in regions
        if item["dataset"] == row["accession"]
        and item["accession"] == row["accession"]
        and item["snp"] == row["SNP"]
    )
    outcome = next(
        item for item in regions
        if item["dataset"] == "GCST90013862"
        and item["accession"] == row["accession"]
        and item["snp"] == row["SNP"]
    )
    start = max(1, int(exposure["start"]))
    end = int(exposure["end"])
    chromosome = exposure["chr"]
    gene = row["gene_symbol"]
    accession = row["accession"]
    exposure_name = f"{gene}_{accession}_chr{chromosome}_{start}_{end}.tsv"
    outcome_name = f"{gene}_GCST90013862_chr{chromosome}_{start}_{end}.tsv"
    shutil.copy2(ROOT / exposure["output"], OUT / exposure_name)
    shutil.copy2(ROOT / outcome["output"], OUT / outcome_name)
    manifest.append({
        "gene": gene,
        "gene_id": row["ensembl_gene_id"].split(".")[0],
        "accession": accession,
        "feature": row["feature_name"],
        "model": exposure["model"],
        "snp": row["SNP"],
        "chr": chromosome,
        "center": exposure["center"],
        "start": str(start),
        "end": str(end),
    })

if {row["gene"] for row in manifest} != selected:
    missing = sorted(selected - {row["gene"] for row in manifest})
    raise RuntimeError(f"Missing mapped locus for biomarkers: {missing}")

fields = ["gene", "gene_id", "accession", "feature", "model", "snp",
          "chr", "center", "start", "end"]
with (OUT / "biomarker_locus_manifest.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest)
print(f"Prepared {len(manifest)} mapped biomarker loci for {len(selected)} genes.")
