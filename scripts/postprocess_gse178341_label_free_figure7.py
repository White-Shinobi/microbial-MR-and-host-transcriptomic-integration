#!/usr/bin/env python3
"""Post-process label-free Figure 7 CellChat and patient correlations."""

from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PB = ATLAS / "pseudobulk"
OUT = ROOT / "results/single_cell/GSE178341/label_free_figure7"
SOURCE = ROOT / "manuscript/source_data"
FOCAL = ["GTF2IRD1", "KIAA1671"]
BIOMARKERS = ["GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"]
TARGET_CACHE = (
    ROOT
    / "results/single_cell/GSE178341/biomarker_localization/"
    "GSE178341_five_biomarker_raw_counts_and_library_sizes.npz"
)
SEED = 20260723
PARTNERS = [
    "B cells", "Epithelial subclusters", "Mast cells", "Myeloid cells",
    "Plasma cells", "Stromal cells", "T/NK/ILC cells",
]
MAP = {
    "B": "B cells",
    "Mast": "Mast cells",
    "Myeloid": "Myeloid cells",
    "Plasma": "Plasma cells",
    "Stromal": "Stromal cells",
    "T/NK/ILC": "T/NK/ILC cells",
}


def is_epi(value: str) -> bool:
    return str(value).startswith("E")


def broad(value: str) -> str:
    return "Epithelial subclusters" if is_epi(value) else MAP.get(
        str(value), str(value)
    )


def contains_component(value: str, gene: str) -> bool:
    return gene in re.split(r"[_+]", str(value))


SOURCE.mkdir(parents=True, exist_ok=True)
sig = pd.read_csv(
    OUT / "label_free_cellchat_BH_FDR_significant_interactions.tsv.gz",
    sep="\t",
)
sig["class"] = sig["tissue"].map({"N": "Normal", "T": "CRC"})
sig["source_broad"] = sig["source"].map(broad)
sig["target_broad"] = sig["target"].map(broad)
epi = sig.loc[sig["involves_epithelial"].astype(bool)].copy()
epi["direction"] = np.select(
    [
        epi["source"].map(is_epi) & epi["target"].map(is_epi),
        epi["source"].map(is_epi),
    ],
    ["Within epithelial", "Outgoing"],
    default="Incoming",
)
epi["partner"] = np.where(
    epi["direction"].eq("Outgoing"),
    epi["target_broad"],
    np.where(
        epi["direction"].eq("Incoming"),
        epi["source_broad"],
        "Epithelial subclusters",
    ),
)
d = (
    epi.groupby(
        ["class", "direction", "partner"], observed=True
    )
    .agg(
        interaction_count=("prob", "size"),
        summed_probability=("prob", "sum"),
    )
    .reset_index()
)
grid = pd.MultiIndex.from_product(
    [
        ["Normal", "CRC"],
        ["Incoming", "Within epithelial", "Outgoing"],
        PARTNERS,
    ],
    names=["class", "direction", "partner"],
).to_frame(index=False)
d = grid.merge(d, how="left").fillna(
    {"interaction_count": 0, "summed_probability": 0}
)
d["interaction_count"] = d["interaction_count"].astype(int)
d.to_csv(
    OUT / "label_free_epithelial_communication_overview.tsv",
    sep="\t",
    index=False,
)

pair = (
    sig.groupby(
        ["class", "source_broad", "target_broad"], observed=True
    )
    .agg(
        interaction_count=("prob", "size"),
        interaction_probability=("prob", "sum"),
    )
    .reset_index()
)
pair = pair.pivot_table(
    index=["source_broad", "target_broad"],
    columns="class",
    values=["interaction_count", "interaction_probability"],
    fill_value=0,
).reset_index()
pair.columns = [
    "_".join(str(x) for x in column if str(x))
    if isinstance(column, tuple)
    else str(column)
    for column in pair.columns
]
for column in [
    "interaction_count_CRC",
    "interaction_count_Normal",
    "interaction_probability_CRC",
    "interaction_probability_Normal",
]:
    if column not in pair:
        pair[column] = 0
pair["delta_count"] = (
    pair["interaction_count_CRC"] - pair["interaction_count_Normal"]
)
pair["delta_probability"] = (
    pair["interaction_probability_CRC"]
    - pair["interaction_probability_Normal"]
)
pair.to_csv(
    OUT / "label_free_cellchat_CRC_minus_normal_broad_difference.tsv",
    sep="\t",
    index=False,
)

lr = (
    epi.groupby(
        ["class", "source", "target", "ligand", "receptor"],
        observed=True,
    )
    .agg(prob=("prob", "sum"), bh_fdr=("p_fdr_bh", "min"))
    .reset_index()
)
wide = lr.pivot_table(
    index=["source", "target", "ligand", "receptor"],
    columns="class",
    values="prob",
    fill_value=0,
).reset_index()
for column in ["Normal", "CRC"]:
    if column not in wide:
        wide[column] = 0
wide["abs_delta"] = (wide["CRC"] - wide["Normal"]).abs()

lr_genes = sorted(
    {
        gene
        for value in pd.concat([epi["ligand"], epi["receptor"]])
        .dropna()
        .astype(str)
        for gene in re.split(r"[_+]", value)
        if gene
    }
)
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
pb_meta = pd.read_csv(
    PB / "label_free_pseudobulk_sample_metadata.tsv", sep="\t"
)
pb = sparse.load_npz(
    PB / "label_free_pseudobulk_raw_counts_genes_by_samples.npz"
).tocsc()
keep_sample = pb_meta["cell_type"].eq("Epithelial") & pb_meta[
    "n_cells"
].ge(20)
pb = pb[:, keep_sample.to_numpy()]
sample_meta = pb_meta.loc[
    keep_sample, ["PID", "SPECIMEN_TYPE", "n_cells"]
].reset_index(drop=True)
symbols = features["gene_symbol"].fillna("").astype(str)
available = set(symbols)
genes_needed = [x for x in FOCAL + lr_genes if x in available]
raw = []
for gene in genes_needed:
    idx = np.flatnonzero(symbols.eq(gene).to_numpy())
    raw.append(np.asarray(pb[idx].sum(axis=0)).ravel())
raw = np.vstack(raw).astype(float)
library = np.asarray(pb.sum(axis=0)).ravel()
expression = np.log2(
    raw / np.maximum(library, 1)[None, :] * 1e6 + 0.5
)
expression = pd.DataFrame(
    expression, index=genes_needed, columns=np.arange(len(sample_meta))
)
delta = {}
for gene in genes_needed:
    table = sample_meta[["PID", "SPECIMEN_TYPE"]].copy()
    table["expression"] = expression.loc[gene].to_numpy()
    gene_wide = table.pivot(
        index="PID", columns="SPECIMEN_TYPE", values="expression"
    ).dropna(subset=["N", "T"])
    delta[gene] = gene_wide["T"] - gene_wide["N"]
rows = []
for biomarker in FOCAL:
    for lr_gene in lr_genes:
        if (
            biomarker not in expression.index
            or lr_gene not in expression.index
            or lr_gene == biomarker
        ):
            continue
        common = delta[biomarker].index.intersection(delta[lr_gene].index)
        rho, p_value = spearmanr(
            delta[biomarker].loc[common],
            delta[lr_gene].loc[common],
        )
        rows.append(
            {
                "biomarker": biomarker,
                "ligand_receptor_gene": lr_gene,
                "n_paired_patients": len(common),
                "spearman_rho_CRC_minus_normal_change": rho,
                "spearman_P": p_value,
            }
        )
correlation = pd.DataFrame(rows)
valid = correlation["spearman_P"].notna()
correlation["spearman_BH_FDR"] = np.nan
correlation.loc[valid, "spearman_BH_FDR"] = multipletests(
    correlation.loc[valid, "spearman_P"], method="fdr_bh"
)[1]
correlation.to_csv(
    OUT / "label_free_focal_biomarker_LR_correlations.tsv",
    sep="\t",
    index=False,
)

# Figure 7F: select one representative CellChat interaction for every
# FDR-significant biomarker–signaling-gene association in Figure 7G. The
# representative interaction is the epithelial-involving ligand–receptor
# route with the largest absolute CRC-minus-normal probability difference.
linked_rows = []
significant_links = correlation.loc[
    correlation["spearman_BH_FDR"].lt(0.05)
].sort_values(
    ["spearman_BH_FDR", "spearman_P", "biomarker"]
).head(6)
for link in significant_links.itertuples(index=False):
    gene = link.ligand_receptor_gene
    candidates = wide.loc[
        wide["ligand"].map(lambda value: contains_component(value, gene))
        | wide["receptor"].map(
            lambda value: contains_component(value, gene)
        )
    ].copy()
    if candidates.empty:
        continue
    best = candidates.sort_values(
        ["abs_delta", "CRC", "Normal"],
        ascending=False,
    ).iloc[0].to_dict()
    best.update(
        {
            "biomarker": link.biomarker,
            "linked_signaling_gene": gene,
            "biomarker_gene_spearman_rho": (
                link.spearman_rho_CRC_minus_normal_change
            ),
            "biomarker_gene_BH_FDR": link.spearman_BH_FDR,
        }
    )
    linked_rows.append(best)
linked = pd.DataFrame(linked_rows)
f = linked.melt(
    id_vars=[
        "biomarker", "linked_signaling_gene",
        "biomarker_gene_spearman_rho", "biomarker_gene_BH_FDR",
        "source", "target", "ligand", "receptor", "abs_delta",
    ],
    value_vars=["Normal", "CRC"],
    var_name="class",
    value_name="prob",
)
f = f.merge(
    lr[
        [
            "class", "source", "target", "ligand", "receptor",
            "bh_fdr",
        ]
    ],
    how="left",
    on=["class", "source", "target", "ligand", "receptor"],
)
f["interaction"] = (
    f["ligand"].astype(str) + "–" + f["receptor"].astype(str)
)
f["route"] = (
    f["source"].astype(str) + " → " + f["target"].astype(str)
)
f.to_csv(
    OUT / "label_free_representative_epithelial_LR_signals.tsv",
    sep="\t",
    index=False,
)

# Figure 7H: patient-level biomarker expression across eight pseudotime
# quantile bins. Raw counts and library sizes are first summed within each
# patient–tissue–bin unit; units containing fewer than 20 sketch cells are
# excluded from the plotted summaries.
embed = pd.read_csv(
    OUT / "label_free_epithelial_embedding_pseudotime.tsv.gz",
    sep="\t",
)
embed["pseudotime_bin"] = pd.qcut(
    embed["graph_pseudotime"],
    q=8,
    labels=False,
    duplicates="drop",
).astype(int)
bin_axis = (
    embed.groupby("pseudotime_bin", observed=True)["graph_pseudotime"]
    .agg(
        pseudotime_min="min",
        pseudotime_median="median",
        pseudotime_max="max",
    )
    .reset_index()
)
target = np.load(TARGET_CACHE)
global_columns = embed["full_matrix_column"].to_numpy(dtype=int)
embed["library_size"] = target["library_size"][global_columns]
trajectory_units = []
for gene in FOCAL:
    gene_index = BIOMARKERS.index(gene)
    working = embed[
        [
            "PID", "SPECIMEN_TYPE", "pseudotime_bin",
            "graph_pseudotime", "library_size",
        ]
    ].copy()
    working["gene_count"] = target["target_counts"][
        gene_index, global_columns
    ]
    unit = (
        working.groupby(
            ["PID", "SPECIMEN_TYPE", "pseudotime_bin"],
            observed=True,
        )
        .agg(
            n_cells=("gene_count", "size"),
            median_cell_pseudotime=("graph_pseudotime", "median"),
            raw_gene_count=("gene_count", "sum"),
            total_library_size=("library_size", "sum"),
        )
        .reset_index()
    )
    unit["gene"] = gene
    unit["log2_CPM_plus_0_5"] = np.log2(
        unit["raw_gene_count"]
        / np.maximum(unit["total_library_size"], 1)
        * 1e6
        + 0.5
    )
    unit["eligible_for_plot"] = unit["n_cells"].ge(20)
    trajectory_units.append(unit)
trajectory_units = pd.concat(trajectory_units, ignore_index=True)
trajectory_units = trajectory_units.merge(
    bin_axis, on="pseudotime_bin", how="left"
)

rng = np.random.default_rng(SEED)
summary_rows = []
eligible_units = trajectory_units.loc[
    trajectory_units["eligible_for_plot"]
].copy()
for keys, frame in eligible_units.groupby(
    ["gene", "SPECIMEN_TYPE", "pseudotime_bin"],
    observed=True,
):
    values = frame["log2_CPM_plus_0_5"].to_numpy()
    bootstrap_medians = np.median(
        rng.choice(
            values,
            size=(1000, len(values)),
            replace=True,
        ),
        axis=1,
    )
    summary_rows.append(
        {
            "gene": keys[0],
            "SPECIMEN_TYPE": keys[1],
            "pseudotime_bin": keys[2],
            "pseudotime_median": frame[
                "pseudotime_median"
            ].iloc[0],
            "n_patient_tissue_units": len(values),
            "median_log2_CPM_plus_0_5": np.median(values),
            "bootstrap_95CI_lower": np.quantile(
                bootstrap_medians, 0.025
            ),
            "bootstrap_95CI_upper": np.quantile(
                bootstrap_medians, 0.975
            ),
        }
    )
trajectory_summary = pd.DataFrame(summary_rows)
trajectory_units.to_csv(
    OUT / "label_free_biomarker_patient_pseudotime_bins.tsv",
    sep="\t",
    index=False,
)
trajectory_summary.to_csv(
    OUT / "label_free_biomarker_pseudotime_summary.tsv",
    sep="\t",
    index=False,
)

for name, table in {
    "Figure_7D_GSE178341_label_free_communication.tsv": d,
    "Figure_7E_GSE178341_label_free_communication_difference.tsv": pair,
    "Figure_7F_GSE178341_label_free_representative_LR.tsv": f,
    "Figure_7G_GSE178341_label_free_biomarker_LR_correlations.tsv": correlation,
    "Figure_7H_GSE178341_biomarker_pseudotime_patient_bins.tsv": (
        trajectory_units
    ),
    "Figure_7H_GSE178341_biomarker_pseudotime_summary.tsv": (
        trajectory_summary
    ),
}.items():
    table.to_csv(SOURCE / name, sep="\t", index=False)
print(
    f"Postprocessed {len(sig):,} significant CellChat rows; "
    f"{len(correlation):,} biomarker–LR tests"
)
