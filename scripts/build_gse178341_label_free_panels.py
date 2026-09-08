#!/usr/bin/env python3
"""Build patient-balanced Figure 6G-H data using label-free mapped types."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PB = ATLAS / "pseudobulk"
SOURCE = ROOT / "manuscript/source_data"
GMT = Path(os.environ.get(
    "MSIGDB_KEGG_GMT",
    ROOT / "data" / "reference" / "msigdb" / "c2.cp.kegg.v7.5.1.symbols.gmt",
))
BIOMARKERS = ["GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"]
CELL_ORDER = [
    "T/NK/ILC", "Epithelial", "Myeloid", "B", "Mast", "Stromal", "Plasma"
]
CELL_LABELS = {
    "T/NK/ILC": "T/NK/ILC cells",
    "Epithelial": "Epithelial cells",
    "Myeloid": "Myeloid cells",
    "B": "B cells",
    "Mast": "Mast cells",
    "Stromal": "Stromal cells",
    "Plasma": "Plasma cells",
}
MIN_CELLS = 20


SOURCE.mkdir(parents=True, exist_ok=True)
patient = pd.read_csv(
    PB / "label_free_patient_tissue_cell_type_biomarker_values.tsv.gz",
    sep="\t",
)
eligible = patient.loc[patient["eligible_min_20_cells"]].copy()
g_summary = (
    eligible.groupby(["gene", "cell_type"], observed=True, sort=False)
    .agg(
        patient_tissue_units=("PID", "size"),
        unique_patients=("PID", "nunique"),
        patient_balanced_detection_rate=("detection_rate", "mean"),
        patient_balanced_mean_expression=(
            "mean_log_normalized_expression", "mean"
        ),
        patient_coverage=("patient_gene_detected", "mean"),
        detection_rate_SEM=(
            "detection_rate",
            lambda x: x.std(ddof=1) / np.sqrt(len(x)),
        ),
        mean_expression_SEM=(
            "mean_log_normalized_expression",
            lambda x: x.std(ddof=1) / np.sqrt(len(x)),
        ),
    )
    .reset_index()
)
g_summary["cell_type_label"] = g_summary["cell_type"].map(CELL_LABELS)
g_summary["cell_order"] = g_summary["cell_type"].map(
    {x: i for i, x in enumerate(CELL_ORDER)}
)
g_summary["gene_order"] = g_summary["gene"].map(
    {x: i for i, x in enumerate(BIOMARKERS)}
)
g_summary = g_summary.sort_values(["cell_order", "gene_order"]).drop(
    columns=["cell_order", "gene_order"]
)
g_summary.to_csv(
    ATLAS / "figure6G_label_free_patient_balanced_biomarker_summary.tsv",
    sep="\t",
    index=False,
)
g_summary.to_csv(
    SOURCE / "Figure_6G_GSE178341_label_free_biomarker_source.tsv",
    sep="\t",
    index=False,
)

features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
pb_meta = pd.read_csv(
    PB / "label_free_pseudobulk_sample_metadata.tsv", sep="\t"
)
pseudobulk = sparse.load_npz(
    PB / "label_free_pseudobulk_raw_counts_genes_by_samples.npz"
).tocsc()
if pseudobulk.shape != (len(features), len(pb_meta)):
    raise RuntimeError("Pseudobulk matrix dimensions do not match metadata")

valid = features["gene_symbol"].notna() & features["gene_symbol"].ne("")
symbol_values = features.loc[valid, "gene_symbol"].astype(str)
unique_symbols = np.sort(symbol_values.unique())
symbol_index = {symbol: i for i, symbol in enumerate(unique_symbols)}
feature_rows = np.flatnonzero(valid.to_numpy())
collapse = sparse.csr_matrix(
    (
        np.ones(len(feature_rows), dtype=np.float32),
        (
            np.asarray(
                [symbol_index[symbol] for symbol in symbol_values],
                dtype=np.int64,
            ),
            feature_rows,
        ),
    ),
    shape=(len(unique_symbols), len(features)),
)
symbol_counts = (collapse @ pseudobulk).tocsc()
eligible_samples = pb_meta["n_cells"].ge(MIN_CELLS).to_numpy()
symbol_counts = symbol_counts[:, eligible_samples]
eligible_meta = pb_meta.loc[eligible_samples].reset_index(drop=True)
library = np.asarray(symbol_counts.sum(axis=0)).ravel()
log_cpm = np.log2(
    symbol_counts.toarray().astype(np.float32)
    / np.maximum(library, 1)[None, :]
    * 1e6
    + 0.5
)

celltype_expression = np.empty(
    (len(unique_symbols), len(CELL_ORDER)), dtype=np.float32
)
celltype_units = []
for j, cell_type in enumerate(CELL_ORDER):
    mask = eligible_meta["cell_type"].eq(cell_type).to_numpy()
    if not mask.any():
        raise RuntimeError(f"No eligible pseudobulk samples for {cell_type}")
    celltype_expression[:, j] = log_cpm[:, mask].mean(axis=1)
    celltype_units.append(int(mask.sum()))

gene_to_row = {gene: i for i, gene in enumerate(unique_symbols)}
pathways: list[tuple[str, list[int]]] = []
with GMT.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        name = (
            fields[0]
            .replace("KEGG_", "")
            .replace("_", " ")
            .title()
        )
        rows = sorted(
            {
                gene_to_row[gene]
                for gene in fields[2:]
                if gene in gene_to_row
            }
        )
        if 15 <= len(rows) <= 500:
            pathways.append((name, rows))

rank_matrix = np.column_stack(
    [
        rankdata(celltype_expression[:, j], method="average")
        / len(unique_symbols)
        for j in range(len(CELL_ORDER))
    ]
).astype(np.float32)
pathway_rows = []
for pathway, rows in pathways:
    scores = rank_matrix[rows].mean(axis=0)
    for j, cell_type in enumerate(CELL_ORDER):
        pathway_rows.append(
            {
                "pathway": pathway,
                "cell_type": cell_type,
                "cell_type_label": CELL_LABELS[cell_type],
                "rank_based_pathway_score": float(scores[j]),
                "mapped_gene_count": len(rows),
                "patient_tissue_units": celltype_units[j],
            }
        )
pathway_scores = pd.DataFrame(pathway_rows)
pathway_scores["between_cell_type_variance"] = pathway_scores.groupby(
    "pathway"
)["rank_based_pathway_score"].transform("var")
pathway_scores.to_csv(
    ATLAS / "figure6H_label_free_all_KEGG_pathway_scores.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
top_pathways = (
    pathway_scores[["pathway", "between_cell_type_variance"]]
    .drop_duplicates()
    .sort_values(
        ["between_cell_type_variance", "pathway"],
        ascending=[False, True],
    )
    .head(10)["pathway"]
    .tolist()
)
heatmap = pathway_scores.loc[
    pathway_scores["pathway"].isin(top_pathways)
].copy()
heatmap["pathway_z_score"] = heatmap.groupby("pathway")[
    "rank_based_pathway_score"
].transform(
    lambda x: (x - x.mean()) / x.std(ddof=1)
    if x.std(ddof=1) > 0
    else 0
)
heatmap["pathway_rank_by_variance"] = heatmap["pathway"].map(
    {x: i + 1 for i, x in enumerate(top_pathways)}
)
heatmap = heatmap.sort_values(
    ["pathway_rank_by_variance", "cell_type"]
)
heatmap.to_csv(
    ATLAS / "figure6H_label_free_top10_pathway_heatmap.tsv",
    sep="\t",
    index=False,
)
heatmap.to_csv(
    SOURCE / "Figure_6H_GSE178341_label_free_pathway_source.tsv",
    sep="\t",
    index=False,
)
print(
    f"Figure 6G: {len(g_summary)} rows; "
    f"Figure 6H: {len(pathways)} KEGG pathways scored"
)
