#!/usr/bin/env python3
"""Apply family-specific BH-FDR to primary IVW and audit MR interpretation."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MR_DIR = ROOT / "results/mr/crc/v2_primary"


def bh_adjust(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return out
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out.loc[order] = np.minimum(adjusted, 1)
    return out


def as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mr-dir", type=Path, default=DEFAULT_MR_DIR)
    args = parser.parse_args()
    mr_file = args.mr_dir / "strict_ref1_mr_results.tsv"
    audit_file = args.mr_dir / "strict_ref1_exposure_audit.tsv"
    if not mr_file.is_file():
        raise SystemExit(f"Missing {mr_file}")

    raw_snapshot = args.mr_dir / "strict_ref1_mr_results.before_v2_fdr.tsv"
    if not raw_snapshot.exists():
        shutil.copy2(mr_file, raw_snapshot)

    mr = pd.read_csv(mr_file, sep="\t")
    primary = as_bool(mr["primary_ivw"])
    mr["nominal_p_lt_0_05"] = primary & mr["pval"].lt(0.05)
    mr["v2_BH_FDR"] = np.nan
    for family, index in mr[primary].groupby("family").groups.items():
        mr.loc[index, "v2_BH_FDR"] = bh_adjust(mr.loc[index, "pval"])
    mr["species_model_class"] = np.where(
        mr["family"].eq("species"),
        np.where(
            mr["model"].astype(str).str.contains("Logistic", case=False, na=False),
            "presence_absence",
            "abundance",
        ),
        "not_applicable",
    )
    mr["species_model_stratified_BH_FDR"] = np.nan
    species_primary = mr[primary & mr["family"].eq("species")]
    for model_class, index in species_primary.groupby("species_model_class").groups.items():
        mr.loc[index, "species_model_stratified_BH_FDR"] = bh_adjust(
            mr.loc[index, "pval"]
        )
    mr["v2_fdr_signal"] = (
        primary & mr["v2_BH_FDR"].lt(0.05)
    )
    # Compatibility with the read-only V1 downstream scripts: in V2 this field
    # means the prespecified FDR gate, not the historical nominal-P gate.
    mr["ref1_candidate"] = mr["v2_fdr_signal"]
    mr.to_csv(mr_file, sep="\t", index=False, na_rep="")
    mr.to_csv(args.mr_dir / "v2_mr_results.tsv", sep="\t", index=False, na_rep="")

    ivw = mr[primary].copy()
    hetero = pd.read_csv(args.mr_dir / "strict_ref1_heterogeneity.tsv", sep="\t")
    egger = pd.read_csv(args.mr_dir / "strict_ref1_egger_pleiotropy.tsv", sep="\t")
    steiger = pd.read_csv(args.mr_dir / "strict_ref1_steiger.tsv", sep="\t")
    presso = pd.read_csv(args.mr_dir / "strict_ref1_mr_presso_summary.tsv", sep="\t")

    ivw_hetero = hetero[
        hetero["method"].astype(str).str.contains("Inverse variance weighted", regex=False)
    ][["accession", "Q", "Q_df", "Q_pval"]].rename(
        columns={"Q": "ivw_Q", "Q_df": "ivw_Q_df", "Q_pval": "ivw_Q_p"}
    )
    egger_small = egger[
        ["accession", "egger_intercept", "se", "pval"]
    ].rename(
        columns={
            "se": "egger_intercept_se",
            "pval": "egger_intercept_p",
        }
    )
    steiger_cols = [
        c
        for c in [
            "accession",
            "correct_causal_direction",
            "steiger_pval",
            "snp_r2.exposure",
            "snp_r2.outcome",
        ]
        if c in steiger.columns
    ]
    sensitivity = (
        ivw.merge(ivw_hetero, on="accession", how="left")
        .merge(egger_small, on="accession", how="left")
        .merge(steiger[steiger_cols], on="accession", how="left")
        .merge(
            presso[["accession", "global_p", "distortion_p", "status"]],
            on="accession",
            how="left",
        )
    )
    sensitivity["heterogeneity_model_consistent"] = np.where(
        sensitivity["ivw_Q_p"].lt(0.05),
        sensitivity["ivw_effect_model"].eq("random_effects"),
        sensitivity["ivw_effect_model"].eq("fixed_effects"),
    )
    sensitivity["egger_intercept_flag"] = sensitivity["egger_intercept_p"].lt(0.05)
    sensitivity["mr_presso_global_flag"] = pd.to_numeric(
        sensitivity["global_p"], errors="coerce"
    ).lt(0.05)
    if "correct_causal_direction" in sensitivity:
        sensitivity["steiger_direction_supported"] = (
            sensitivity["correct_causal_direction"].astype("boolean")
            & pd.to_numeric(sensitivity["steiger_pval"], errors="coerce").lt(0.05)
        )
    else:
        sensitivity["steiger_direction_supported"] = pd.NA
    sensitivity["interpretation"] = np.where(
        sensitivity["v2_fdr_signal"],
        "FDR-supported genetically proxied association; not proof of causality or mediation",
        "No family-wise FDR support",
    )
    sensitivity.to_csv(
        args.mr_dir / "v2_primary_ivw_sensitivity_audit.tsv",
        sep="\t",
        index=False,
        na_rep="",
    )

    audit = pd.read_csv(audit_file, sep="\t")
    summary = pd.DataFrame(
        [
            ("deduplicated_candidate_exposures", len(audit)),
            ("mr_eligible_pre_harmonization", int(audit["eligible_pre_harmonization"].sum())),
            ("mr_eligible_post_harmonization", int(audit["eligible_final"].sum())),
            ("primary_ivw_tests", len(ivw)),
            ("nominal_primary_ivw_p_lt_0_05", int(ivw["nominal_p_lt_0_05"].sum())),
            ("family_specific_bh_fdr_q_lt_0_05", int(ivw["v2_fdr_signal"].sum())),
            (
                "species_abundance_primary_ivw_tests",
                int(
                    (
                        ivw["family"].eq("species")
                        & ivw["species_model_class"].eq("abundance")
                    ).sum()
                ),
            ),
            (
                "species_presence_absence_primary_ivw_tests",
                int(
                    (
                        ivw["family"].eq("species")
                        & ivw["species_model_class"].eq("presence_absence")
                    ).sum()
                ),
            ),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(args.mr_dir / "v2_mr_run_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(
        ivw[
            [
                "family",
                "feature_name",
                "accession",
                "nsnp",
                "OR",
                "CI_95_low",
                "CI_95_high",
                "pval",
                "v2_BH_FDR",
                "v2_fdr_signal",
            ]
        ]
        .sort_values(["family", "pval"])
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
