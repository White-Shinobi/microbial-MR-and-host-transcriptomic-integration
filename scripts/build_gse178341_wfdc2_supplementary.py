#!/usr/bin/env python3
"""Build the patient-level WFDC2 pseudobulk supplementary figure for GSE178341."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PSEUDOBULK = RESULT / "pseudobulk"
SOURCE = ROOT / "manuscript/source_data"
FIGURE = ROOT / "manuscript/figures/supplementary"

CELL_ORDER = ["T/NK/ILC", "Epithelial", "Myeloid", "B", "Mast", "Stromal", "Plasma"]
CELL_LABELS = {
    "T/NK/ILC": "T/NK/ILC",
    "Epithelial": "Epithelial",
    "Myeloid": "Myeloid",
    "B": "B",
    "Mast": "Mast",
    "Stromal": "Stromal",
    "Plasma": "Plasma",
}
COLORS = {"Normal": "#B2182B", "CRC": "#2166AC"}
RNG = np.random.default_rng(20260722)


def star(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def violin(ax, values, position, color):
    values = np.asarray(values, dtype=float)
    if len(values) > 1 and not np.allclose(values, values[0]):
        parts = ax.violinplot(values, positions=[position], widths=0.34, showextrema=False, bw_method=0.30)
        body = parts["bodies"][0]
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_linewidth(0.45)
        body.set_alpha(0.88)
    q1, med, q3 = np.quantile(values, [0.25, 0.50, 0.75])
    ax.plot([position, position], [q1, q3], color="black", lw=2.2, solid_capstyle="butt", zorder=4)
    ax.scatter([position], [med], s=16, c="white", edgecolors="black", linewidths=0.55, zorder=5)


SOURCE.mkdir(parents=True, exist_ok=True)
FIGURE.mkdir(parents=True, exist_ok=True)

pb = pd.read_csv(
    PSEUDOBULK / "label_free_patient_tissue_cell_type_biomarker_values.tsv.gz",
    sep="\t",
)
edger = pd.read_csv(
    PSEUDOBULK / "label_free_paired_edger_biomarker_results_35_tests.tsv",
    sep="\t",
)

rows = []
for cell_type in CELL_ORDER:
    part = pb.loc[
        pb["gene"].eq("WFDC2")
        & pb["cell_type"].eq(cell_type)
        & pb["eligible_min_20_cells"],
        ["PID", "SPECIMEN_TYPE", "log2_cpm"],
    ]
    wide = part.pivot(index="PID", columns="SPECIMEN_TYPE", values="log2_cpm")
    if not {"T", "N"}.issubset(wide.columns):
        continue
    wide = wide.dropna(subset=["T", "N"])
    model = edger.loc[edger["gene"].eq("WFDC2") & edger["cell_type"].eq(cell_type)].iloc[0]
    for pid, vals in wide.iterrows():
        for tissue in ["N", "T"]:
            rows.append({
                "gene": "WFDC2",
                "cell_type_code": cell_type,
                "cell_type": CELL_LABELS[cell_type] + " cells",
                "PID": pid,
                "tissue": tissue,
                "group": "Normal" if tissue == "N" else "CRC",
                "log2_CPM_plus_0_5": vals[tissue],
                "n_paired_patients": int(model["n_paired_patients"]),
                "edgeR_log2_fold_change_CRC_vs_normal": model["edgeR_log2_fold_change_CRC_vs_normal"],
                "edgeR_QL_P": model["edgeR_QL_P"],
                "edgeR_BH_FDR_35_tests": model["edgeR_BH_FDR_35_tests"],
                "target_gene_passed_filterByExpr": bool(model["target_gene_passed_filterByExpr"]),
            })

plot_data = pd.DataFrame(rows)
source_path = SOURCE / "Figure_S6_GSE178341_WFDC2_paired_pseudobulk_source.tsv"
plot_data.to_csv(source_path, sep="\t", index=False)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
fig, ax = plt.subplots(figsize=(11.2, 4.8), facecolor="white")
for idx, cell_type in enumerate(CELL_ORDER):
    d = plot_data.loc[plot_data["cell_type_code"].eq(cell_type)]
    for group, offset in [("Normal", -0.19), ("CRC", 0.19)]:
        values = d.loc[d["group"].eq(group), "log2_CPM_plus_0_5"].to_numpy()
        violin(ax, values, idx + offset, COLORS[group])
        jitter = RNG.normal(0, 0.025, len(values))
        ax.scatter(idx + offset + jitter, values, s=11, c=COLORS[group], edgecolors="white", linewidths=0.25, alpha=0.82, zorder=3)
    first = d.iloc[0]
    ymax = d["log2_CPM_plus_0_5"].max()
    ax.text(idx, ymax + 0.30, star(first["edgeR_BH_FDR_35_tests"]), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.text(idx, ymax + 0.08, f"n={int(first['n_paired_patients'])}", ha="center", va="bottom", fontsize=6.5, color="#555555")

ax.set_xticks(range(len(CELL_ORDER)), [CELL_LABELS[c] for c in CELL_ORDER], rotation=25, ha="right")
ax.set_ylabel("Patient pseudobulk expression\nlog2(CPM + 0.5)")
ax.set_title("WFDC2 expression across major cell types in paired GSE178341 tissues", fontweight="bold")
ax.grid(axis="y", color="#E2E2E2", lw=0.55)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.scatter([], [], c=COLORS["Normal"], label="Normal")
ax.scatter([], [], c=COLORS["CRC"], label="CRC")
ax.legend(frameon=False, loc="upper left", ncol=2)
fig.text(0.995, 0.015, "Stars: edgeR quasi-likelihood BH-FDR across 35 gene × cell-type tests", ha="right", va="bottom", fontsize=7, color="#444444")
fig.tight_layout(rect=[0.02, 0.05, 0.99, 0.98])

for suffix in ["png", "pdf"]:
    fig.savefig(FIGURE / f"Figure_S6_GSE178341_WFDC2_pseudobulk.{suffix}", dpi=400 if suffix == "png" else None, bbox_inches="tight")
plt.close(fig)

print(source_path)
print(FIGURE / "Figure_S6_GSE178341_WFDC2_pseudobulk.pdf")
