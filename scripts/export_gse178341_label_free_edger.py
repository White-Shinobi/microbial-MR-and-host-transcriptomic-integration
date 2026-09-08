#!/usr/bin/env python3
"""Export paired label-free GSE178341 pseudobulk matrices for edgeR."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results/single_cell/GSE178341/label_free_atlas"
PB = ATLAS / "pseudobulk"
EXPORT = PB / "edger_input"
MIN_CELLS = 20
CELL_TYPES = {
    "Epithelial": "Epithelial",
    "T/NK/ILC": "T_NK_ILC",
    "Myeloid": "Myeloid",
    "B": "B",
    "Plasma": "Plasma",
    "Stromal": "Stromal",
    "Mast": "Mast",
}


EXPORT.mkdir(parents=True, exist_ok=True)
features = pd.read_csv(ATLAS / "features.tsv.gz", sep="\t")
sample_meta = pd.read_csv(
    PB / "label_free_pseudobulk_sample_metadata.tsv", sep="\t"
)
counts = sparse.load_npz(
    PB / "label_free_pseudobulk_raw_counts_genes_by_samples.npz"
).tocsc()
if counts.shape != (len(features), len(sample_meta)):
    raise RuntimeError("Pseudobulk matrix and metadata dimensions differ")

audit = []
for cell_type, slug in CELL_TYPES.items():
    meta = sample_meta.loc[
        sample_meta["cell_type"].eq(cell_type)
        & sample_meta["n_cells"].ge(MIN_CELLS)
    ].copy()
    paired = (
        meta.groupby("PID")["SPECIMEN_TYPE"]
        .agg(set)
        .loc[lambda x: x.apply(lambda value: {"N", "T"}.issubset(value))]
        .index
    )
    meta = meta.loc[meta["PID"].isin(paired)].copy()
    meta["tissue_order"] = meta["SPECIMEN_TYPE"].map({"N": 0, "T": 1})
    meta = meta.sort_values(["PID", "tissue_order"]).drop(
        columns="tissue_order"
    )
    if not (meta.groupby("PID").size() == 2).all():
        raise RuntimeError(f"Non-paired samples remained for {cell_type}")
    meta["sample_name"] = (
        meta["PID"].astype(str) + "__" + meta["SPECIMEN_TYPE"].astype(str)
    )
    columns = meta["pseudobulk_index"].to_numpy(dtype=int)
    dense = counts[:, columns].toarray()
    dense = np.rint(dense).astype(np.int32)
    frame = features[
        ["feature_index", "ensembl_id", "gene_symbol"]
    ].copy()
    for column, sample_name in enumerate(meta["sample_name"]):
        frame[sample_name] = dense[:, column]
    frame.to_csv(
        EXPORT / f"{slug}_counts.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    meta[
        ["sample_name", "PID", "SPECIMEN_TYPE", "cell_type", "n_cells"]
    ].to_csv(
        EXPORT / f"{slug}_samples.tsv", sep="\t", index=False
    )
    audit.append(
        {
            "cell_type": cell_type,
            "file_slug": slug,
            "paired_patients": len(paired),
            "samples": len(meta),
            "minimum_cells": int(meta["n_cells"].min()),
            "maximum_cells": int(meta["n_cells"].max()),
        }
    )

pd.DataFrame(audit).to_csv(
    EXPORT / "label_free_edger_input_audit.tsv", sep="\t", index=False
)
print(pd.DataFrame(audit).to_string(index=False))
