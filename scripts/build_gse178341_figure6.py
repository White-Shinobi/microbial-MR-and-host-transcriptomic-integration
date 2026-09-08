#!/usr/bin/env python3
"""Assemble Figure 6 from the coherent label-free GSE178341 workflow."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PSEUDOBULK = RESULT / "pseudobulk"
SOURCE = ROOT / "manuscript/source_data"
FIGURE = ROOT / "manuscript/figures/main"

BIOMARKERS = ["GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"]
DISPLAY_GENES = ["GTF2IRD1", "KIAA1671"]
CELL_ORDER = ["T/NK/ILC", "Epithelial", "Myeloid", "B", "Mast", "Stromal", "Plasma"]
CELL_LABELS = {
    "T/NK/ILC": "T/NK/ILC cells",
    "Epithelial": "Epithelial cells",
    "Myeloid": "Myeloid cells",
    "B": "B cells",
    "Mast": "Mast cells",
    "Stromal": "Stromal cells",
    "Plasma": "Plasma cells",
}
CELL_COLORS = {
    "B": "#E76F51",
    "Epithelial": "#C49A00",
    "Mast": "#00A651",
    "Myeloid": "#00A9B7",
    "Plasma": "#8E5EA2",
    "Stromal": "#5794F2",
    "T/NK/ILC": "#E64DCB",
}
TISSUE_COLORS = {"T": "#B2182B", "N": "#2166AC"}
SEED = 20260722
RNG = np.random.default_rng(SEED)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def significance(q: float) -> str:
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def draw_violin(ax: plt.Axes, values: np.ndarray, position: float, color: str, width: float = 0.34) -> None:
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or np.allclose(values, values[0]):
        ax.plot([position - width / 3, position + width / 3], [values[0], values[0]], color="black", lw=0.65)
        return
    violin = ax.violinplot(
        values,
        positions=[position],
        widths=width,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.30,
    )
    body = violin["bodies"][0]
    body.set_facecolor(color)
    body.set_edgecolor("black")
    body.set_linewidth(0.38)
    body.set_alpha(0.92)


SOURCE.mkdir(parents=True, exist_ok=True)
FIGURE.mkdir(parents=True, exist_ok=True)

hvg = pd.read_csv(RESULT / "label_free_hvg_statistics.tsv.gz", sep="\t")
pc = pd.read_csv(RESULT / "label_free_principal_component_statistics.tsv", sep="\t")
embedding = pd.read_csv(RESULT / "label_free_sketch_embedding_and_clusters.tsv.gz", sep="\t")
cluster_search = pd.read_csv(RESULT / "label_free_leiden_silhouette_selection.tsv", sep="\t")
cluster_annotation = pd.read_csv(RESULT / "label_free_cluster_annotation.tsv", sep="\t")
g_source = pd.read_csv(RESULT / "figure6G_label_free_patient_balanced_biomarker_summary.tsv", sep="\t")
h_source = pd.read_csv(RESULT / "figure6H_label_free_top10_pathway_heatmap.tsv", sep="\t")
edger = pd.read_csv(PSEUDOBULK / "label_free_paired_edger_biomarker_results_35_tests.tsv", sep="\t")
pseudobulk = pd.read_csv(
    PSEUDOBULK / "label_free_patient_tissue_cell_type_biomarker_values.tsv.gz",
    sep="\t",
)

# Freeze source data used by every displayed panel.
hvg.to_csv(SOURCE / "Figure_6A_GSE178341_HVG_source.tsv.gz", sep="\t", index=False, compression="gzip")
embedding[["cellID", "PID", "SPECIMEN_TYPE", "PC1", "PC2"]].to_csv(
    SOURCE / "Figure_6B_GSE178341_PCA_source.tsv.gz", sep="\t", index=False, compression="gzip"
)
pc.to_csv(SOURCE / "Figure_6CD_GSE178341_PC_selection_QC_source.tsv", sep="\t", index=False)
embedding.to_csv(
    SOURCE / "Figure_6EF_GSE178341_embedding_source.tsv.gz", sep="\t", index=False, compression="gzip"
)
cluster_search.to_csv(SOURCE / "Figure_6E_GSE178341_cluster_selection_source.tsv", sep="\t", index=False)
cluster_annotation.to_csv(SOURCE / "Figure_6F_GSE178341_cluster_annotation_source.tsv", sep="\t", index=False)
g_source.to_csv(SOURCE / "Figure_6G_GSE178341_patient_balanced_biomarker_source.tsv", sep="\t", index=False)
h_source.to_csv(SOURCE / "Figure_6H_GSE178341_pathway_heatmap_source.tsv", sep="\t", index=False)

# Create paired patient plotting data for I-J and freeze edgeR labels.
ij_rows = []
ij_significance = []
for gene in DISPLAY_GENES:
    for cell_type in CELL_ORDER:
        part = pseudobulk.loc[
            pseudobulk["gene"].eq(gene)
            & pseudobulk["cell_type"].eq(cell_type)
            & pseudobulk["eligible_min_20_cells"],
            ["PID", "SPECIMEN_TYPE", "log2_cpm"],
        ]
        wide = part.pivot(index="PID", columns="SPECIMEN_TYPE", values="log2_cpm")
        if {"T", "N"}.issubset(wide.columns):
            wide = wide.dropna(subset=["T", "N"])
            for patient, row in wide.iterrows():
                for tissue in ["T", "N"]:
                    ij_rows.append(
                        {
                            "gene": gene,
                            "cell_type_code": cell_type,
                            "cell_type": CELL_LABELS[cell_type],
                            "PID": patient,
                            "tissue": tissue,
                            "group": "CRC" if tissue == "T" else "Normal",
                            "log2_CPM_plus_0_5": row[tissue],
                        }
                    )
        model = edger.loc[edger["gene"].eq(gene) & edger["cell_type"].eq(cell_type)].iloc[0]
        ij_significance.append(
            {
                "gene": gene,
                "cell_type_code": cell_type,
                "cell_type": CELL_LABELS[cell_type],
                "n_paired_patients": int(model["n_paired_patients"]),
                "edgeR_log2_fold_change_CRC_vs_normal": model["edgeR_log2_fold_change_CRC_vs_normal"],
                "edgeR_QL_P": model["edgeR_QL_P"],
                "edgeR_BH_FDR_35_tests": model["edgeR_BH_FDR_35_tests"],
                "paired_wilcoxon_P_sensitivity": model["paired_wilcoxon_P"],
                "paired_wilcoxon_BH_FDR_sensitivity": model["paired_wilcoxon_BH_FDR_35_tests"],
                "figure_label": significance(model["edgeR_BH_FDR_35_tests"]),
            }
        )
ij_data = pd.DataFrame(ij_rows)
ij_sig = pd.DataFrame(ij_significance)
ij_data.to_csv(SOURCE / "Figure_6IJ_GSE178341_paired_patient_pseudobulk_source.tsv", sep="\t", index=False)
ij_sig.to_csv(SOURCE / "Figure_6IJ_GSE178341_edgeR_significance_source.tsv", sep="\t", index=False)
edger.to_csv(SOURCE / "Figure_6_GSE178341_all_35_pseudobulk_tests_source.tsv", sep="\t", index=False)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7.5,
        "axes.titlesize": 8.5,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
    }
)
fig = plt.figure(figsize=(19, 17.5), facecolor="white")
outer = GridSpec(4, 6, figure=fig, height_ratios=[0.90, 1.0, 1.15, 1.0], hspace=0.58, wspace=0.72)

# A: HVG selection.
ax = fig.add_subplot(outer[0, 0:2])
eligible = np.isfinite(hvg["hvg_estimation_standardized_dispersion"])
ax.scatter(
    np.log10(
        hvg.loc[
            eligible & ~hvg["selected_hvg"],
            "hvg_estimation_mean_log_normalized",
        ].clip(lower=1e-8)
    ),
    hvg.loc[
        eligible & ~hvg["selected_hvg"],
        "hvg_estimation_standardized_dispersion",
    ],
    s=2,
    c="#1F1F1F",
    alpha=0.34,
    linewidths=0,
    rasterized=True,
)
ax.scatter(
    np.log10(
        hvg.loc[
            hvg["selected_hvg"],
            "hvg_estimation_mean_log_normalized",
        ].clip(lower=1e-8)
    ),
    hvg.loc[
        hvg["selected_hvg"],
        "hvg_estimation_standardized_dispersion",
    ],
    s=3.5,
    c="#E41A1C",
    alpha=0.70,
    linewidths=0,
    rasterized=True,
)
for _, row in hvg.loc[hvg["selected_hvg"]].nsmallest(10, "hvg_rank").iterrows():
    ax.text(
        np.log10(max(row["hvg_estimation_mean_log_normalized"], 1e-8)),
        row["hvg_estimation_standardized_dispersion"] + 0.08,
        str(row["gene_symbol"]),
        fontsize=6,
        ha="center",
    )
ax.set_title("Top 2,000 highly variable genes", fontweight="bold")
ax.set_xlabel("Average log-normalized expression (log10)")
ax.set_ylabel("Expression-bin standardized dispersion")
ax.legend(
    handles=[
        Line2D([], [], marker="o", linestyle="", color="#1F1F1F", markersize=3, label="Other genes"),
        Line2D([], [], marker="o", linestyle="", color="#E41A1C", markersize=3, label="HVGs (n=2,000)"),
    ],
    loc="lower right",
    frameon=False,
)
panel_label(ax, "a")

# B: PCA.
ax = fig.add_subplot(outer[0, 2:4])
for tissue, label in [("N", "Normal"), ("T", "CRC")]:
    d = embedding.loc[embedding["SPECIMEN_TYPE"].eq(tissue)]
    ax.scatter(d["PC1"], d["PC2"], s=1.1, alpha=0.35, linewidths=0, color=TISSUE_COLORS[tissue], label=label, rasterized=True)
ax.set_title("PCA distribution in the label-free sketch", fontweight="bold")
ax.set_xlabel(f"PC1 ({100*pc.iloc[0]['explained_variance_ratio']:.1f}%)")
ax.set_ylabel(f"PC2 ({100*pc.iloc[1]['explained_variance_ratio']:.1f}%)")
ax.legend(frameon=False, markerscale=4)
panel_label(ax, "b")

# C: elbow / retained PCs.
ax = fig.add_subplot(outer[0, 4:6])
ax.plot(pc["PC"], pc["standard_deviation"], "o-", color="#222222", markersize=2.8, linewidth=0.8)
ax.axvline(30, color="#D73027", linestyle="--", linewidth=1)
ax.text(30.7, pc["standard_deviation"].max() * 0.84, "30 PCs retained", color="#D73027", fontsize=7)
ax.set_title("Elbow plot and PC retention", fontweight="bold")
ax.set_xlabel("Principal component")
ax.set_ylabel("Standard deviation")
panel_label(ax, "c")

# D: technical-covariate PC QC.
ax = fig.add_subplot(outer[1, 0:2])
pc30 = pc.head(30)
ax.plot(pc30["PC"], pc30["spearman_rho_log_total_UMI"].abs(), marker="o", ms=2.5, lw=0.8, label="Total UMI")
ax.plot(pc30["PC"], pc30["spearman_rho_detected_features"].abs(), marker="o", ms=2.5, lw=0.8, label="Detected genes")
ax.plot(pc30["PC"], pc30["spearman_rho_mitochondrial_percent"].abs(), marker="o", ms=2.5, lw=0.8, label="Mitochondrial %")
ax.axhline(0.30, color="#777777", linestyle="--", linewidth=0.7)
ax.set_ylim(0, max(0.52, pc30.filter(like="spearman").abs().max().max() * 1.08))
ax.set_title("PC technical-covariate quality control", fontweight="bold")
ax.set_xlabel("Retained principal component")
ax.set_ylabel("Absolute Spearman correlation")
ax.legend(frameon=False, ncol=3, loc="upper right")
panel_label(ax, "d")

# E: unsupervised clusters.
ax = fig.add_subplot(outer[1, 2:4])
cluster_values = sorted(embedding["label_free_leiden_cluster"].unique())
cluster_cmap = plt.get_cmap("tab20", len(cluster_values))
for i, cluster in enumerate(cluster_values):
    d = embedding.loc[embedding["label_free_leiden_cluster"].eq(cluster)]
    ax.scatter(d["UMAP1"], d["UMAP2"], s=1.0, alpha=0.55, linewidths=0, color=cluster_cmap(i), rasterized=True)
    ax.text(d["UMAP1"].median(), d["UMAP2"].median(), str(cluster), fontsize=6.5, fontweight="bold", ha="center", va="center")
selected_k = int(cluster_search.loc[cluster_search["selected"], "clusters"].iloc[0])
ax.set_title(f"Unsupervised clusters (k={selected_k})", fontweight="bold")
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
panel_label(ax, "e")

# F: annotated major cell types.
ax = fig.add_subplot(outer[1, 4:6])
for cell_type in CELL_ORDER:
    d = embedding.loc[
        embedding["label_free_marker_cell_type"].eq(cell_type)
    ]
    ax.scatter(d["UMAP1"], d["UMAP2"], s=1.0, alpha=0.55, linewidths=0, color=CELL_COLORS[cell_type], label=CELL_LABELS[cell_type], rasterized=True)
    ax.text(d["UMAP1"].median(), d["UMAP2"].median(), CELL_LABELS[cell_type].replace(" cells", ""), fontsize=5.8, ha="center", va="center")
ax.set_title("Annotated major cell types", fontweight="bold")
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.legend(frameon=False, bbox_to_anchor=(1.01, 0.5), loc="center left", markerscale=4)
panel_label(ax, "f")

# G: biomarker localization using equal-weight patient–tissue summaries.
ax = fig.add_subplot(outer[2, 0:3])
g = g_source.copy()
g["x"] = g["gene"].map({gene: i for i, gene in enumerate(BIOMARKERS)})
g["y"] = g["cell_type"].map({cell: i for i, cell in enumerate(CELL_ORDER)})
scatter = ax.scatter(
    g["x"],
    g["y"],
    s=18 + 850 * g["patient_balanced_detection_rate"],
    c=g["patient_balanced_mean_expression"],
    cmap="YlOrRd",
    edgecolor="#333333",
    linewidth=0.35,
)
ax.set_xticks(range(len(BIOMARKERS)), BIOMARKERS, rotation=42, ha="right", fontstyle="italic")
ax.set_yticks(range(len(CELL_ORDER)), [CELL_LABELS[cell] for cell in CELL_ORDER])
ax.invert_yaxis()
ax.set_title("Localization of the five biomarkers across cell types", fontweight="bold")
ax.set_xlabel("")
ax.set_ylabel("")
cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label("Mean log-normalized expression", fontsize=6.5)
for pct in [1, 5, 10, 25]:
    ax.scatter([], [], s=18 + 8.5 * pct, facecolor="#BDBDBD", edgecolor="#333333", label=f"{pct}%")
ax.legend(title="Detection rate", frameon=False, bbox_to_anchor=(0.5, -0.20), loc="upper center", ncol=4)
panel_label(ax, "g")

# H: clustered pathway heatmap with dendrograms.
hgrid = GridSpecFromSubplotSpec(
    2,
    3,
    subplot_spec=outer[2, 3:6],
    width_ratios=[0.62, 0.12, 0.26],
    height_ratios=[0.18, 0.82],
    wspace=0.02,
    hspace=0.02,
)
ax_top = fig.add_subplot(hgrid[0, 0])
ax_corner = fig.add_subplot(hgrid[0, 1:3]); ax_corner.axis("off")
ax_h = fig.add_subplot(hgrid[1, 0])
ax_left = fig.add_subplot(hgrid[1, 1])
ax_labels = fig.add_subplot(hgrid[1, 2])
hmat = h_source.pivot(
    index="pathway", columns="cell_type", values="pathway_z_score"
).reindex(columns=CELL_ORDER)
row_link = linkage(hmat.to_numpy(), method="complete", metric="euclidean")
col_link = linkage(hmat.to_numpy().T, method="complete", metric="euclidean")
row_order = leaves_list(row_link)
col_order = leaves_list(col_link)
dendrogram(col_link, ax=ax_top, no_labels=True, color_threshold=0, above_threshold_color="#333333")
ax_top.axis("off")
dendrogram(row_link, ax=ax_left, orientation="right", no_labels=True, color_threshold=0, above_threshold_color="#333333")
ax_left.axis("off")
hordered = hmat.iloc[row_order, col_order]
im = ax_h.imshow(hordered.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax_h.set_xticks(range(len(hordered.columns)), [CELL_LABELS[x] for x in hordered.columns], rotation=42, ha="right")
short_pathway_labels = {
    "Intestinal Immune Network For Iga Production": "Intestinal IgA network",
    "Maturity Onset Diabetes Of The Young": "Maturity-onset diabetes of the young",
    "Pentose And Glucuronate Interconversions": "Pentose/glucuronate interconversion",
    "Ascorbate And Aldarate Metabolism": "Ascorbate/aldarate metabolism",
    "Ecm Receptor Interaction": "ECM–receptor interaction",
    "Graft Versus Host Disease": "Graft-versus-host disease",
}
ax_h.set_yticks([])
ax_top.set_title("KEGG pathway scores across cell types", fontweight="bold", pad=5)
for xline in np.arange(-0.5, len(hordered.columns), 1):
    ax_h.axvline(xline, color="white", lw=0.4)
for yline in np.arange(-0.5, len(hordered.index), 1):
    ax_h.axhline(yline, color="white", lw=0.4)
ax_labels.set_xlim(0, 1)
ax_labels.set_ylim(len(hordered.index) - 0.5, -0.5)
for row_index, pathway in enumerate(hordered.index):
    ax_labels.text(
        0.0,
        row_index,
        short_pathway_labels.get(pathway, pathway),
        ha="left",
        va="center",
        fontsize=5.3,
    )
ax_labels.axis("off")
cbar_ax = ax_labels.inset_axes([0.88, 0.08, 0.08, 0.84])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Pathway-wise z score", fontsize=6.5)
ax_top.text(-0.12, 1.22, "h", transform=ax_top.transAxes, fontsize=12, fontweight="bold", va="top")

# I-J: paired patient pseudobulk distributions with edgeR FDR labels.
for panel_index, (gene, span) in enumerate(zip(DISPLAY_GENES, [(0, 3), (3, 6)])):
    ax = fig.add_subplot(outer[3, span[0]:span[1]])
    gene_data = ij_data.loc[ij_data["gene"].eq(gene)]
    offsets = {"T": -0.18, "N": 0.18}
    for cell_i, cell_type in enumerate(CELL_ORDER):
        for tissue in ["T", "N"]:
            values = gene_data.loc[
                gene_data["cell_type_code"].eq(cell_type) & gene_data["tissue"].eq(tissue),
                "log2_CPM_plus_0_5",
            ].to_numpy(dtype=float)
            position = cell_i + offsets[tissue]
            draw_violin(ax, values, position, TISSUE_COLORS[tissue])
            jitter = RNG.normal(0, 0.022, len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                s=7,
                facecolor=TISSUE_COLORS[tissue],
                edgecolor="white",
                linewidth=0.22,
                alpha=0.65,
                zorder=4,
            )
            ax.scatter(position, np.median(values), s=20, facecolor="white", edgecolor="black", linewidth=0.7, zorder=6)
    gene_min = gene_data["log2_CPM_plus_0_5"].min()
    gene_max = gene_data["log2_CPM_plus_0_5"].max()
    span_y = max(gene_max - gene_min, 1)
    for cell_i, cell_type in enumerate(CELL_ORDER):
        q = ij_sig.loc[
            ij_sig["gene"].eq(gene) & ij_sig["cell_type_code"].eq(cell_type),
            "edgeR_BH_FDR_35_tests",
        ].iloc[0]
        ax.text(cell_i, gene_max + 0.07 * span_y, significance(q), ha="center", va="bottom", fontsize=7)
    ax.set_xlim(-0.55, len(CELL_ORDER) - 0.45)
    ax.set_ylim(gene_min - 0.08 * span_y, gene_max + 0.18 * span_y)
    ax.set_xticks(range(len(CELL_ORDER)), [CELL_LABELS[cell] for cell in CELL_ORDER], rotation=42, ha="right")
    ax.set_ylabel(f"Pseudobulk expression of {gene}\nlog2(CPM + 0.5)")
    ax.set_xlabel("")
    ax.legend(
        handles=[
            Patch(facecolor=TISSUE_COLORS["T"], edgecolor="black", label="CRC"),
            Patch(facecolor=TISSUE_COLORS["N"], edgecolor="black", label="Normal"),
        ],
        title="group",
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=2,
    )
    panel_label(ax, "i" if panel_index == 0 else "j")

for ax in fig.axes:
    if hasattr(ax, "spines") and ax not in [ax_top, ax_left, ax_corner]:
        for side in ["top", "right"]:
            if side in ax.spines:
                ax.spines[side].set_visible(False)

fig.text(
    0.5,
    0.006,
    "Panels A–F use a label-free patient–tissue-stratified leverage-score sketch; labels were transferred to all cells before panels G–J. I–J stars: edgeR QL BH-FDR across 35 tests.",
    ha="center",
    fontsize=7,
)

png = FIGURE / "Figure_6_GSE178341_single_cell_key_cells_preview.png"
pdf = FIGURE / "Figure_6_GSE178341_single_cell_key_cells_preview.pdf"
fig.savefig(png, dpi=320, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(
    FIGURE / "Figure_6_GSE178341_single_cell_key_cells.png",
    dpi=320,
    bbox_inches="tight",
)
fig.savefig(
    FIGURE / "Figure_6_GSE178341_single_cell_key_cells.pdf",
    bbox_inches="tight",
)
plt.close(fig)

print(f"Built complete Figure 6: {png}")
