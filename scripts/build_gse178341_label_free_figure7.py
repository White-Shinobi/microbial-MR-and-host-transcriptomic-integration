#!/usr/bin/env python3
"""Build Figure 7 from the coherent label-free GSE178341 epithelial chain."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/single_cell/GSE178341/label_free_figure7"
FIG = ROOT / "manuscript/figures/main"
SUPP = ROOT / "manuscript/figures/supplementary"
SOURCE = ROOT / "manuscript/source_data"
NORMAL = "#2A9D8F"
CRC = "#D1495B"
NORMAL_BLUE = "#2F6BBD"
SEED = 20260723


def panel_tag(ax, tag):
    ax.text(
        -0.13, 1.07, tag, transform=ax.transAxes,
        fontsize=13, fontweight="bold", va="top"
    )


def sig(q):
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


FIG.mkdir(parents=True, exist_ok=True)
SUPP.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)
embed = pd.read_csv(
    OUT / "label_free_epithelial_embedding_pseudotime.tsv.gz", sep="\t"
)
composition = pd.read_csv(
    OUT / "label_free_epithelial_subcluster_patient_composition.tsv",
    sep="\t",
)
comp_stats = pd.read_csv(
    OUT / "label_free_epithelial_subcluster_paired_statistics.tsv",
    sep="\t",
)
d = pd.read_csv(
    OUT / "label_free_epithelial_communication_overview.tsv", sep="\t"
)
e = pd.read_csv(
    OUT / "label_free_cellchat_CRC_minus_normal_broad_difference.tsv",
    sep="\t",
)
f = pd.read_csv(
    OUT / "label_free_representative_epithelial_LR_signals.tsv",
    sep="\t",
)
g = pd.read_csv(
    OUT / "label_free_focal_biomarker_LR_correlations.tsv", sep="\t"
)
patient_pt = pd.read_csv(
    OUT / "label_free_epithelial_patient_pseudotime.tsv", sep="\t"
)
pt_stats = pd.read_csv(
    OUT / "label_free_epithelial_patient_pseudotime_statistics.tsv",
    sep="\t",
).iloc[0]
trajectory_units = pd.read_csv(
    OUT / "label_free_biomarker_pseudotime_GAMM_units.tsv",
    sep="\t",
)
trajectory_predictions = pd.read_csv(
    OUT / "label_free_biomarker_pseudotime_GAMM_predictions.tsv",
    sep="\t",
)
trajectory_statistics = pd.read_csv(
    OUT / "label_free_biomarker_pseudotime_GAMM_statistics.tsv",
    sep="\t",
)

subclusters = sorted(embed["epithelial_subcluster_label"].unique())
code = {label: label.split(" ", 1)[0] for label in subclusters}
embed["subcluster_code"] = embed["epithelial_subcluster_label"].map(code)
embed["tissue"] = embed["SPECIMEN_TYPE"].map({"N": "Normal", "T": "CRC"})
palette = dict(
    zip(
        [code[x] for x in subclusters],
        sns.color_palette("tab10", n_colors=len(subclusters)),
    )
)
rng = np.random.default_rng(SEED)
plot_embed = embed.iloc[rng.permutation(len(embed))].copy()

for name, table in {
    "Figure_7ABH_GSE178341_label_free_epithelial_embedding.tsv.gz": embed,
    "Figure_7C_GSE178341_label_free_epithelial_composition.tsv": composition,
    "Figure_7C_GSE178341_label_free_epithelial_statistics.tsv": comp_stats,
    "Figure_7H_GSE178341_label_free_patient_pseudotime.tsv": patient_pt,
}.items():
    table.to_csv(
        SOURCE / name,
        sep="\t",
        index=False,
        compression="gzip" if name.endswith(".gz") else None,
    )

fig = plt.figure(figsize=(18, 13.2))
gs = GridSpec(
    3, 12, figure=fig,
    height_ratios=[1.02, 0.96, 0.92],
    hspace=0.55, wspace=0.78,
)

# A: de novo epithelial subclusters.
ax = fig.add_subplot(gs[0, 0:3])
for label in subclusters:
    short = code[label]
    part = plot_embed.loc[
        plot_embed["epithelial_subcluster_label"].eq(label)
    ]
    ax.scatter(
        part["UMAP1"], part["UMAP2"], s=2.0, alpha=0.7,
        color=palette[short], rasterized=True
    )
    ax.text(
        part["UMAP1"].median(), part["UMAP2"].median(), short,
        fontsize=6.8, fontweight="bold", ha="center", va="center"
    )
ax.set(
    title="De novo epithelial subclusters",
    xlabel="UMAP 1", ylabel="UMAP 2"
)
panel_tag(ax, "a")

# B: tissue distribution.
sub = gs[0, 3:7].subgridspec(1, 2, wspace=0.08)
for j, tissue in enumerate(["Normal", "CRC"]):
    bx = fig.add_subplot(sub[0, j])
    other = plot_embed.loc[~plot_embed["tissue"].eq(tissue)]
    part = plot_embed.loc[plot_embed["tissue"].eq(tissue)]
    bx.scatter(
        other["UMAP1"], other["UMAP2"], s=0.8,
        color="#E7E7E7", alpha=0.22, rasterized=True
    )
    bx.scatter(
        part["UMAP1"], part["UMAP2"], s=1.8,
        color=NORMAL if tissue == "Normal" else CRC,
        alpha=0.7, rasterized=True
    )
    bx.set_title(tissue, fontsize=9, fontweight="bold")
    bx.set_xlabel("UMAP 1")
    if j == 0:
        bx.set_ylabel("UMAP 2")
        panel_tag(bx, "b")
    else:
        bx.set_ylabel("")
        bx.set_yticklabels([])
fig.text(
    0.35, 0.93, "Normal and CRC epithelial states",
    ha="center", fontsize=10, fontweight="bold"
)

# C: paired patient-level composition.
subc = gs[0, 7:12].subgridspec(2, 5, hspace=0.58, wspace=0.38)
for i, label in enumerate(subclusters):
    cx = fig.add_subplot(subc[i // 5, i % 5])
    part = composition.loc[
        composition["epithelial_subcluster_label"].eq(label)
    ]
    wide = part.pivot(
        index="PID", columns="SPECIMEN_TYPE", values="fraction"
    ).dropna(subset=["N", "T"])
    for _, row in wide.iterrows():
        cx.plot(
            [0, 1], [row["N"], row["T"]],
            color="#B7B7B7", lw=0.4, alpha=0.5
        )
    cx.scatter(
        rng.normal(0, 0.025, len(wide)), wide["N"],
        s=6, color=NORMAL, alpha=0.75
    )
    cx.scatter(
        1 + rng.normal(0, 0.025, len(wide)), wide["T"],
        s=6, color=CRC, alpha=0.75
    )
    q = comp_stats.loc[
        comp_stats["epithelial_subcluster_label"].eq(label),
        "paired_wilcoxon_BH_FDR",
    ].iloc[0]
    cx.text(
        0.5, 0.96, sig(q), transform=cx.transAxes,
        ha="center", va="top", fontsize=8, fontweight="bold"
    )
    cx.set_title(code[label], fontsize=7.4, fontweight="bold")
    cx.set_xticks([0, 1], ["N", "CRC"], fontsize=6)
    cx.tick_params(axis="y", labelsize=5.7)
    cx.set_ylim(bottom=0)
    if i == 0:
        panel_tag(cx, "c")
fig.text(
    0.81, 0.93,
    "Patient-level epithelial-subcluster composition\n"
    "Paired Wilcoxon; BH across 10 subclusters",
    ha="center", fontsize=9.2, fontweight="bold"
)

# D: epithelial-centered communication.
ax = fig.add_subplot(gs[1, 0:4])
row_order = [
    f"{tissue} · {direction}"
    for tissue in ["Normal", "CRC"]
    for direction in ["Incoming", "Within epithelial", "Outgoing"]
]
d["row"] = pd.Categorical(
    d["class"] + " · " + d["direction"], row_order
)
pivot = d.pivot(
    index="row", columns="partner", values="summed_probability"
)
pivot = pd.concat(
    [
        pivot.iloc[:3],
        pd.DataFrame(
            [[np.nan] * pivot.shape[1]],
            index=[""],
            columns=pivot.columns,
        ),
        pivot.iloc[3:],
    ]
)
sns.heatmap(
    pivot, ax=ax, cmap="magma", linewidths=0.45,
    linecolor="white", cbar=False, mask=pivot.isna()
)
for iy, row in enumerate(pivot.index):
    if row == "":
        continue
    for ix, partner in enumerate(pivot.columns):
        n = int(
            d.loc[
                d["row"].eq(row) & d["partner"].eq(partner),
                "interaction_count",
            ].iloc[0]
        )
        if n:
            ax.text(
                ix + 0.5, iy + 0.5, str(n),
                color="white", fontsize=6.5, fontweight="bold",
                ha="center", va="center"
            )
ax.set(
    title="Epithelial-centered communication",
    xlabel="", ylabel=""
)
ax.tick_params(axis="x", rotation=42, labelsize=6.5)
ax.tick_params(axis="y", labelsize=6.5, length=0)
for label in ax.get_yticklabels():
    if label.get_text().startswith("Normal"):
        label.set_color(NORMAL_BLUE)
        label.set_fontweight("bold")
    elif label.get_text().startswith("CRC"):
        label.set_color(CRC)
        label.set_fontweight("bold")
ax.axhline(3, color="#D9D9D9", linewidth=0.8)
ax.axhline(4, color="#D9D9D9", linewidth=0.8)
panel_tag(ax, "d")

# E: broad communication contrast.
ax = fig.add_subplot(gs[1, 4:7])
order = [
    "B cells", "Epithelial subclusters", "Mast cells",
    "Myeloid cells", "Plasma cells", "Stromal cells",
    "T/NK/ILC cells",
]
delta = e.pivot(
    index="source_broad", columns="target_broad",
    values="delta_probability"
).reindex(index=order, columns=order).fillna(0)
limit = max(np.abs(delta.to_numpy()).max(), 1e-8)
hm = sns.heatmap(
    delta, ax=ax, cmap="RdBu_r", center=0,
    vmin=-limit, vmax=limit, linewidths=0.4,
    linecolor="white", cbar_kws={"label": "CRC − normal probability"}
)
ax.set(
    title="CRC–normal communication change",
    xlabel="", ylabel="Source"
)
ax.tick_params(axis="x", rotation=52, labelsize=6.2)
ax.tick_params(axis="y", labelsize=6.2)
panel_tag(ax, "e")

# F: biomarker-linked ligand-receptor signals shown as a dumbbell plot.
ax = fig.add_subplot(gs[1, 8:12])
def short_cell(value):
    value = str(value)
    if value == "B":
        return "B cell"
    if value.startswith("E") and " (" in value:
        return value.split(" ", 1)[0]
    return value


f["display_route"] = (
    f["source"].map(short_cell) + " → " + f["target"].map(short_cell)
)
f["label"] = (
    f["biomarker"] + ": " + f["interaction"]
    + "  |  " + f["display_route"]
)
dumbbell = f.pivot_table(
    index="label", columns="class", values="prob", fill_value=0
)
for column in ["Normal", "CRC"]:
    if column not in dumbbell:
        dumbbell[column] = 0
dumbbell = dumbbell.sort_values(
    ["CRC", "Normal"], ascending=True
)
for iy, (_, row) in enumerate(dumbbell.iterrows()):
    ax.plot(
        [row["Normal"], row["CRC"]], [iy, iy],
        color="#B8B8B8", linewidth=1.05, zorder=1
    )
ax.scatter(
    dumbbell["Normal"], np.arange(len(dumbbell)),
    s=30, facecolors="white", edgecolors=NORMAL_BLUE,
    linewidths=1.15, label="Normal", zorder=3
)
ax.scatter(
    dumbbell["CRC"], np.arange(len(dumbbell)),
    s=30, color=CRC, edgecolors="white",
    linewidths=0.45, label="CRC", zorder=4
)
ax.set_yticks(
    range(len(dumbbell)), dumbbell.index.tolist(), fontsize=5.3
)
maximum = max(dumbbell[["Normal", "CRC"]].to_numpy().max(), 1e-8)
ax.set_xlim(-maximum * 0.035, maximum * 1.08)
ax.set_title("Biomarker-linked epithelial signaling")
ax.set_xlabel("Communication probability")
ax.ticklabel_format(
    axis="x", style="sci", scilimits=(-3, -3), useMathText=True
)
ax.tick_params(axis="x", labelsize=6.2)
ax.grid(axis="x", color="#E5E5E5", lw=0.7)
ax.legend(
    loc="lower right", frameon=False, fontsize=6.2,
    handletextpad=0.35
)
panel_tag(ax, "f")

# G: patient-level focal biomarker-signaling correlations.
ax = fig.add_subplot(gs[2, 0:5])
g["abs_rho"] = g["spearman_rho_CRC_minus_normal_change"].abs()
top_genes = (
    g.groupby("ligand_receptor_gene")["abs_rho"]
    .max().sort_values(ascending=False).head(18).index
)
show = g.loc[g["ligand_receptor_gene"].isin(top_genes)].copy()
mat = show.pivot(
    index="biomarker", columns="ligand_receptor_gene",
    values="spearman_rho_CRC_minus_normal_change"
).reindex(index=["GTF2IRD1", "KIAA1671"])
sns.heatmap(
    mat, ax=ax, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
    linewidths=0.45, linecolor="white",
    cbar_kws={"label": "Spearman ρ", "shrink": 0.65}
)
for iy, biomarker in enumerate(mat.index):
    for ix, gene in enumerate(mat.columns):
        q_row = show.loc[
            show["biomarker"].eq(biomarker)
            & show["ligand_receptor_gene"].eq(gene)
        ]
        if not q_row.empty and q_row.iloc[0]["spearman_BH_FDR"] < 0.05:
            ax.text(
                ix + 0.5, iy + 0.5, "*",
                ha="center", va="center",
                fontsize=7, fontweight="bold"
            )
ax.set(
    title="Patient-level biomarker–signaling correlations",
    xlabel="Ligand/receptor gene", ylabel=""
)
ax.tick_params(axis="x", rotation=55, labelsize=6.2)
ax.tick_params(axis="y", rotation=0, labelsize=7)
panel_tag(ax, "g")

# H: exploratory pseudotime and patient-aware biomarker trajectories.
ax = fig.add_subplot(gs[2, 5:8])
scatter = ax.scatter(
    plot_embed["UMAP1"], plot_embed["UMAP2"],
    c=plot_embed["graph_pseudotime"], s=1.8,
    cmap="viridis", alpha=0.72, rasterized=True
)
root = int(pt_stats["root_cluster"])
root_cells = embed.loc[
    embed["epithelial_leiden_cluster"].eq(root)
]
ax.scatter(
    root_cells["UMAP1"].median(),
    root_cells["UMAP2"].median(),
    marker="*", s=90, color="#F4A261",
    edgecolor="black", linewidth=0.5
)
ax.set(
    title="Exploratory epithelial state trajectory",
    xlabel="UMAP 1", ylabel="UMAP 2"
)
fig.colorbar(
    scatter, ax=ax, fraction=0.05, pad=0.14,
    label="Graph pseudotime", orientation="horizontal"
)
panel_tag(ax, "h")

for gene, span in [
    ("GTF2IRD1", (8, 10)),
    ("KIAA1671", (10, 12)),
]:
    ax = fig.add_subplot(gs[2, span[0]:span[1]])
    gene_units = trajectory_units.loc[
        trajectory_units["gene"].eq(gene)
    ].copy()
    gene_predictions = trajectory_predictions.loc[
        trajectory_predictions["gene"].eq(gene)
    ].copy()
    gene_statistics = trajectory_statistics.loc[
        trajectory_statistics["gene"].eq(gene)
    ].iloc[0]
    support_min = gene_statistics["observed_pseudotime_min"]
    support_max = gene_statistics["observed_pseudotime_max"]
    gene_predictions = gene_predictions.loc[
        gene_predictions["graph_pseudotime"].between(
            support_min, support_max
        )
    ]
    for tissue, label, color, marker in [
        ("Normal", "Normal", NORMAL_BLUE, "o"),
        ("CRC", "CRC", CRC, "s"),
    ]:
        observed = gene_units.loc[
            gene_units["tissue"].eq(tissue)
        ]
        fitted = gene_predictions.loc[
            gene_predictions["tissue"].eq(tissue)
        ].sort_values("graph_pseudotime")
        ax.scatter(
            observed["graph_pseudotime"],
            observed["expression"],
            s=5.5, color=color, alpha=0.10,
            edgecolors="none", rasterized=True, zorder=1,
        )
        x = fitted["graph_pseudotime"].to_numpy(dtype=float)
        y = fitted["fit"].to_numpy(dtype=float)
        lower = fitted["confidence_95_lower"].to_numpy(dtype=float)
        upper = fitted["confidence_95_upper"].to_numpy(dtype=float)
        ax.fill_between(
            x, lower, upper, color=color, alpha=0.18,
            linewidth=0, zorder=2,
        )
        ax.plot(
            x, y, color=color, linewidth=2.05,
            label=label, zorder=3,
        )
    p_value = gene_statistics["overall_curve_P"]
    p_text = (
        r"$P_{\mathrm{curve}} < 2.2 \times 10^{-16}$"
        if p_value < 2.2e-16
        else rf"$P_{{\mathrm{{curve}}}} = {p_value:.2g}$"
    )
    ax.set_title(
        gene, fontsize=9.2, fontstyle="italic", fontweight="bold"
    )
    ax.text(
        0.04, 0.96,
        f"Paired-patient GAMM\n{p_text}",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=6.2,
    )
    # Display only the patient-level support used by the GAMM. The underlying
    # graph pseudotime remains scaled to 0-100, but extending the panel beyond
    # the eligible patient-tissue-bin range creates unsupported blank space.
    ax.set_xlim(support_min, support_max)
    tick_start = int(np.ceil(support_min / 10.0) * 10)
    tick_end = int(np.floor(support_max / 10.0) * 10)
    ax.set_xticks(np.arange(tick_start, tick_end + 1, 10))
    ax.set_xlabel("Epithelial graph pseudotime")
    if gene == "GTF2IRD1":
        ax.set_ylabel(
            "Patient-level expression\n"
            r"$\log_{2}(\mathrm{CPM}+0.5)$"
        )
        ax.legend(
            loc="upper right", frameon=False, fontsize=6.2,
            handlelength=1.6,
        )
    else:
        ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=6.2)
    ax.grid(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

for axis in fig.axes:
    if hasattr(axis, "spines"):
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
fig.subplots_adjust(
    left=0.055, right=0.985, bottom=0.055, top=0.97
)
for extension in ["png", "pdf"]:
    fig.savefig(
        FIG
        / f"Figure_7_GSE178341_label_free_epithelial_state_remodeling."
        f"{extension}",
        dpi=320 if extension == "png" else None,
        bbox_inches="tight",
    )
plt.close(fig)
print("Built label-free Figure 7")
