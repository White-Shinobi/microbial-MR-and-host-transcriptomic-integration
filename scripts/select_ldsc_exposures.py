#!/usr/bin/env python3
"""Select V2 MR exposures without applying an LDSC eligibility gate.

Every GWAS satisfying the prespecified genome-wide-significance candidate rule
is eligible.  LDSC estimates are descriptive and are used only to rank exact
cross-cohort duplicates (same biological trait and phenotype class).  Missing
or invalid LDSC estimates never exclude an otherwise eligible GWAS.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
LDSC_RESULTS = ROOT / "results/ldsc/candidate_ldsc_results.tsv"
INVENTORY = ROOT / "config/v2_candidate_exposure_inventory.tsv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit",
        type=Path,
        default=ROOT / "results/ldsc/v2_ldsc_selection_audit.tsv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "config/v2_ldsc_selected_exposures.tsv",
    )
    args = parser.parse_args()

    inventory = pd.read_csv(INVENTORY, sep="\t")
    ldsc = pd.read_csv(LDSC_RESULTS, sep="\t")
    ldsc_columns = [
        "accession",
        "status",
        "seconds",
        "ldsc_h2",
        "ldsc_h2_se",
        "ldsc_h2_z",
        "ldsc_intercept",
        "ldsc_intercept_se",
        "ldsc_mean_chisq",
        "ldsc_n_sumstats",
        "ldsc_n_merged",
    ]
    ldsc = ldsc[[column for column in ldsc_columns if column in ldsc.columns]]
    df = inventory.merge(ldsc, on="accession", how="left", validate="one_to_one")
    df["candidate_pass"] = df["candidate_rule_pass"].fillna(False).astype(bool)
    df["ldsc_valid_h2"] = pd.to_numeric(df["ldsc_h2"], errors="coerce").between(
        0, 1, inclusive="neither"
    )
    df["ldsc_z_pass"] = pd.to_numeric(df["ldsc_h2_z"], errors="coerce").gt(1.64)
    df["ldsc_gate_applied"] = False
    df["ldsc_gate_pass"] = pd.NA
    df["equivalence_key"] = df["species_equivalence_key"].fillna("")
    module = df["trait_family"].eq("module")
    df.loc[module, "equivalence_key"] = df.loc[
        module, "module_equivalence_key"
    ].fillna("")

    df["dedup_rank"] = pd.NA
    eligible = df[df["candidate_pass"]].copy()
    eligible["_ldsc_z_sort"] = pd.to_numeric(
        eligible["ldsc_h2_z"], errors="coerce"
    )
    eligible = eligible.sort_values(
        ["equivalence_key", "_ldsc_z_sort", "sample_size", "accession"],
        ascending=[True, False, False, True],
        na_position="last",
    )
    eligible["dedup_rank"] = eligible.groupby("equivalence_key").cumcount() + 1
    rank_map = eligible.set_index("accession")["dedup_rank"]
    df["dedup_rank"] = df["accession"].map(rank_map).astype("Int64")
    df["selected_for_mr"] = df["candidate_pass"] & df["dedup_rank"].eq(1)

    def reason(row: pd.Series) -> str:
        if not row["candidate_pass"]:
            return "Failed P<5e-8 candidate rule"
        if row["dedup_rank"] != 1:
            winner = eligible.loc[
                eligible["equivalence_key"].eq(row["equivalence_key"])
                & eligible["dedup_rank"].eq(1),
                "accession",
            ]
            return (
                "Equivalent trait retained from "
                + (winner.iloc[0] if len(winner) else "higher-h2-Z GWAS")
            )
        if pd.isna(row.get("ldsc_h2_z")):
            return "Selected for MR; LDSC Z unavailable and not used as a gate"
        return "Selected for MR; LDSC used descriptively only"

    df["selection_reason"] = df.apply(reason, axis=1)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.audit, sep="\t", index=False, na_rep="")

    selected = df[df["selected_for_mr"]].copy()
    manifest = pd.DataFrame(
        {
            "accession": selected["accession"],
            "family": selected["trait_family"],
            "feature_id": selected["trait_id"],
            "feature_name": selected["trait_name"],
            "model": selected["abundance_model"],
            "phenotype_class": selected["phenotype_class"],
            "analysis_role": "primary",
            "exposure_n": selected["sample_size"],
            "cohort": selected["cohort"],
            "ldsc_h2": selected["ldsc_h2"],
            "ldsc_h2_se": selected["ldsc_h2_se"],
            "ldsc_h2_z": selected["ldsc_h2_z"],
            "source_gwas_path": selected["local_path"],
        }
    ).sort_values(["family", "feature_name", "accession"])
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest, sep="\t", index=False, na_rep="")

    unified = ROOT / "data/raw/gwas/microbiome_selected"
    unified.mkdir(parents=True, exist_ok=True)
    for row in manifest.itertuples(index=False):
        destination = unified / row.accession
        source_path = Path(row.source_gwas_path)
        if not source_path.is_absolute():
            source_path = ROOT / source_path
        source = source_path.parent
        if destination.is_symlink():
            if destination.resolve() != source.resolve():
                destination.unlink()
                destination.symlink_to(source, target_is_directory=True)
        elif destination.exists():
            raise RuntimeError(f"Refusing to replace existing path: {destination}")
        else:
            destination.symlink_to(source, target_is_directory=True)

    print(
        df.groupby(["cohort", "trait_family"])
        .agg(
            candidates=("accession", "size"),
            candidate_rule_pass=("candidate_pass", "sum"),
            selected_for_mr=("selected_for_mr", "sum"),
        )
        .to_string()
    )
    print(f"Wrote {args.audit}")
    print(f"Wrote {args.manifest}")


if __name__ == "__main__":
    main()
