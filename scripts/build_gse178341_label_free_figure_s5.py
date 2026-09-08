#!/usr/bin/env python3
"""Build supplementary diagnostics for the label-free epithelial analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
EPI = ROOT / "results/single_cell/GSE178341/label_free_epithelial"
F8 = ROOT / "results/single_cell/GSE178341/label_free_figure7"
OUT = ROOT / "manuscript/figures/supplementary"
NORMAL = "#2A9D8F"
NORMAL_BLUE = "#2F6BBD"
CRC = "#D1495B"
SEED = 20260723


OUT.mkdir(parents=True, exist_ok=True)
score = pd.read_csv(
    EPI / "epithelial_label_free_marker_scores.tsv", sep="\t"
).set_index("cluster")
confusion = pd.read_csv(
    EPI / "epithelial_label_free_vs_publisher_subtype_confusion.tsv",
    sep="\t", index_col=0
)
confusion = confusion.div(confusion.sum(axis=1), axis=0)
composition = pd.read_csv(
    F8 / "label_free_epithelial_subcluster_paired_statistics.tsv",
    sep="\t",
)
embed = pd.read_csv(
    F8 / "label_free_epithelial_embedding_pseudotime.tsv.gz", sep="\t"
)
embed["tissue"] = embed["SPECIMEN_TYPE"].map({"N": "Normal", "T": "CRC"})
patient_pt = pd.read_csv(
    F8 / "label_free_epithelial_patient_pseudotime.tsv", sep="\t"
)
pt_stats = pd.read_csv(
    F8 / "label_free_epithelial_patient_pseudotime_statistics.tsv",
    sep="\t",
).iloc[0]

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
sns.heatmap(
    score, ax=axes[0, 0], cmap="RdBu_r", center=0,
    linewidths=0.4, cbar_kws={"label": "Marker evidence"}
)
axes[0, 0].set_title("A  Canonical marker scores")
axes[0, 0].set_xlabel("Marker-defined state")
axes[0, 0].set_ylabel("De novo cluster")

sns.heatmap(
    confusion, ax=axes[0, 1], cmap="Blues", vmin=0, vmax=1,
    cbar_kws={"label": "Within-cluster proportion"}
)
axes[0, 1].set_title("B  Post-hoc publisher-subtype concordance")
axes[0, 1].set_xlabel("Publisher subtype (audit only)")
axes[0, 1].set_ylabel("De novo subcluster")
axes[0, 1].tick_params(axis="x", rotation=55, labelsize=6)

axes[0, 2].scatter(
    composition["median_normal"], composition["median_CRC"],
    c=composition["paired_wilcoxon_BH_FDR"].lt(0.05),
    cmap="coolwarm", s=60, edgecolor="black", linewidth=0.4
)
limit = max(
    composition["median_normal"].max(),
    composition["median_CRC"].max(),
)
axes[0, 2].plot([0, limit], [0, limit], "--", color="gray")
for _, row in composition.iterrows():
    axes[0, 2].text(
        row["median_normal"], row["median_CRC"],
        row["epithelial_subcluster_label"].split(" ", 1)[0],
        fontsize=7
    )
axes[0, 2].set(
    title="C  Patient-level subcluster redistribution",
    xlabel="Median normal fraction", ylabel="Median CRC fraction"
)

sns.histplot(
    data=embed, x="graph_pseudotime", hue="tissue",
    palette={"Normal": NORMAL, "CRC": CRC}, element="step",
    stat="density", common_norm=False, ax=axes[1, 0]
)
axes[1, 0].set_title("D  Graph-pseudotime distributions")

sns.boxplot(
    data=embed, x="epithelial_subcluster_label",
    y="graph_pseudotime", hue="tissue",
    palette={"Normal": NORMAL, "CRC": CRC},
    fliersize=0, ax=axes[1, 1]
)
axes[1, 1].tick_params(axis="x", rotation=55, labelsize=6)
axes[1, 1].set_title("E  Pseudotime by de novo subcluster")
axes[1, 1].legend(fontsize=7)

wide = patient_pt.pivot(
    index="PID", columns="SPECIMEN_TYPE", values="median_pseudotime"
).dropna(subset=["N", "T"])
rng = np.random.default_rng(SEED)
for _, row in wide.iterrows():
    axes[1, 2].plot(
        [0, 1], [row["N"], row["T"]],
        color="#A7A7A7", lw=0.75, alpha=0.7,
    )
axes[1, 2].scatter(
    rng.normal(0, 0.025, len(wide)), wide["N"],
    s=20, color=NORMAL_BLUE, alpha=0.85, label="Normal",
)
axes[1, 2].scatter(
    1 + rng.normal(0, 0.025, len(wide)), wide["T"],
    s=20, color=CRC, alpha=0.85, label="CRC",
)
axes[1, 2].set_xticks([0, 1], ["Normal", "CRC"])
axes[1, 2].set_ylabel("Patient median graph pseudotime")
axes[1, 2].set_title(
    "F  Paired patient pseudotime shift\n"
    f"P = {pt_stats['paired_wilcoxon_P']:.2g}; "
    f"n = {int(pt_stats['n_pairs'])} pairs"
)

fig.suptitle(
    "Figure S5  Label-free GSE178341 epithelial-state diagnostics",
    fontsize=15, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
for extension in ["png", "pdf"]:
    fig.savefig(
        OUT / f"Figure_S5_GSE178341_epithelial_diagnostics.{extension}",
        dpi=300 if extension == "png" else None,
        bbox_inches="tight",
    )
plt.close(fig)
print("Built updated Figure S5")
