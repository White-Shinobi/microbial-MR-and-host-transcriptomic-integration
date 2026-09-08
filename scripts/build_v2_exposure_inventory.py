#!/usr/bin/env python3
"""Build the V2 Swedish/HUNT exposure universe without using CRC results.

Species abundance and Swedish presence/absence traits are eligible.  The
phenotype class is retained explicitly so that abundance and presence/absence
are never deduplicated against one another.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_METADATA_DIR = Path(os.environ.get(
    "SOURCE_METADATA_DIR", ROOT / "config" / "source_metadata"
))
HUNT_XLSX = ROOT / "data/raw/metadata/hunt/HUNT_Supplementary_Tables_1_35.xlsx"


def read_sheet(name: str, header_row: int) -> pd.DataFrame:
    out = pd.read_excel(HUNT_XLSX, sheet_name=name, header=header_row)
    out = out.dropna(axis=1, how="all").dropna(axis=0, how="all")
    return out


def normalize_species(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    return re.sub(r"\s+", " ", text.strip()).casefold()


def accession(number: int) -> str:
    return f"GCST{number:08d}"


def gwas_url(acc: str) -> str:
    n = int(acc[4:])
    lo = (n // 1000) * 1000 + 1
    hi = lo + 999
    bucket = f"GCST{lo:08d}-GCST{hi:08d}"
    return (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        f"{bucket}/{acc}/harmonised/{acc}.h.tsv.gz"
    )


def local_path(cohort: str, acc: str, phenotype_class: str) -> Path:
    folder = "microbiome_swedish" if cohort == "Swedish" else "microbiome_hunt"
    shared_path = ROOT / "data/raw/gwas" / folder / acc / f"{acc}.h.tsv.gz"
    if shared_path.is_file():
        return shared_path
    if cohort == "Swedish" and phenotype_class == "presence_absence":
        return (
            ROOT
            / "data/raw/gwas/microbiome_swedish_v2_presence"
            / acc
            / f"{acc}.h.tsv.gz"
        )
    return shared_path


def build_swedish() -> pd.DataFrame:
    src = SOURCE_METADATA_DIR / "swedish_gws_eligible_mr_exposures.tsv"
    df = pd.read_csv(src, sep="\t")
    keep = df["family"].isin(["species", "metabolic_module"])
    df = df.loc[keep].copy()
    df["cohort"] = "Swedish"
    df["trait_family"] = df["family"].replace({"metabolic_module": "module"})
    df["trait_id"] = df["feature_id"]
    df["trait_name"] = df["feature_name"]
    df["abundance_model"] = df["model"]
    df["phenotype_class"] = np.where(
        df["trait_family"].eq("species"),
        df["phenotype"].replace(
            {"abundance": "abundance", "presence": "presence_absence"}
        ),
        "module",
    )
    df["sample_size"] = 16017
    df["candidate_gws_count"] = df["gws_variant_count"].astype(int)
    df["candidate_lead_rsid"] = df["lead_rsid"]
    df["candidate_lead_p"] = df["lead_p"].astype(float)
    df["candidate_source"] = df["eligibility_source"]
    df["published_ldsc_h2"] = np.nan
    df["published_ldsc_se"] = np.nan
    df["published_ldsc_z"] = np.nan
    df["published_ldsc_p"] = np.nan
    return df


def build_hunt() -> pd.DataFrame:
    species_all = read_sheet("Supplementary Table 1", 5)
    species_all = species_all.rename(columns={"MGS": "trait_id", "species": "trait_name"})
    species_all = species_all[["trait_id", "trait_name"]].reset_index(drop=True)
    species_all = species_all[
        species_all["trait_id"].astype(str).str.match(r"^hMGS\.\d+$", na=False)
    ].reset_index(drop=True)
    if len(species_all) != 546:
        raise ValueError(f"Expected 546 HUNT species, found {len(species_all)}")
    species_all["accession"] = [accession(90666543 + i) for i in range(len(species_all))]
    species_all["_key"] = species_all["trait_name"].map(normalize_species)

    species_hits = read_sheet("Supplementary Table 6", 3)
    species_hits = species_hits.rename(columns={"Outcome species": "trait_name"})
    species_hits = species_hits[
        species_hits["trait_name"].notna()
        & pd.to_numeric(species_hits["Chr"], errors="coerce").notna()
    ].copy()
    species_hits["_key"] = species_hits["trait_name"].map(normalize_species)
    hit_summary = (
        species_hits.groupby("_key", as_index=False)
        .agg(
            candidate_gws_count=("P-value", "size"),
            candidate_lead_rsid=("RSID", "first"),
            candidate_lead_p=("P-value", "min"),
        )
    )
    if len(hit_summary) != 98:
        raise ValueError(f"Expected 98 unique HUNT species hits, found {len(hit_summary)}")

    ldsc = read_sheet("Supplementary Table 5 ", 3)
    ldsc = ldsc.rename(
        columns={
            "Species": "trait_name",
            "Heritability (h2)": "published_ldsc_h2",
            "SE": "published_ldsc_se",
            "P-value": "published_ldsc_p",
        }
    )
    ldsc["_key"] = ldsc["trait_name"].map(normalize_species)
    ldsc["published_ldsc_z"] = (
        pd.to_numeric(ldsc["published_ldsc_h2"], errors="coerce")
        / pd.to_numeric(ldsc["published_ldsc_se"], errors="coerce")
    )

    sp = species_all.merge(hit_summary, on="_key", how="inner", validate="one_to_one")
    sp = sp.merge(
        ldsc[
            [
                "_key",
                "published_ldsc_h2",
                "published_ldsc_se",
                "published_ldsc_z",
                "published_ldsc_p",
            ]
        ],
        on="_key",
        how="left",
        validate="one_to_one",
    )
    sp["cohort"] = "HUNT"
    sp["trait_family"] = "species"
    sp["abundance_model"] = "continuous relative abundance"
    sp["phenotype_class"] = "abundance"
    sp["sample_size"] = 12652
    sp["candidate_source"] = "HUNT Supplementary Tables 1, 5 and 6"

    module_all = read_sheet("Supplementary Table 24", 3)
    module_all = module_all.rename(
        columns={"module": "trait_id", "name": "trait_name", "type": "module_type"}
    )
    module_all = module_all[["trait_id", "trait_name", "module_type"]].reset_index(drop=True)
    module_all = module_all[
        module_all["trait_id"].astype(str).str.match(r"^M\d{5}$", na=False)
    ].reset_index(drop=True)
    if len(module_all) != 461:
        raise ValueError(f"Expected 461 HUNT modules, found {len(module_all)}")
    module_all["accession"] = [accession(90667089 + i) for i in range(len(module_all))]

    module_hits = pd.DataFrame(
        [
            ("M00443", "rs4988235", 2.6e-18),
            ("M00233", "rs4988235", 3.1e-18),
            ("M00244", "rs182549", 4.0e-15),
            ("M00168", "rs182549", 2.7e-11),
            ("M00228", "rs635634", 1.3e-10),
            ("M00719", "rs35866622", 1.2e-10),
            ("M00518", "rs804264", 6.1e-11),
            ("M00747", "rs804273", 1.0e-10),
        ],
        columns=["trait_id", "candidate_lead_rsid", "candidate_lead_p"],
    )
    mo = module_all.merge(module_hits, on="trait_id", how="inner", validate="one_to_one")
    if len(mo) != 8:
        raise ValueError(f"Expected 8 HUNT module hits, found {len(mo)}")
    mo["candidate_gws_count"] = 1
    mo["cohort"] = "HUNT"
    mo["trait_family"] = "module"
    mo["abundance_model"] = "continuous relative abundance"
    mo["phenotype_class"] = "module"
    mo["sample_size"] = 12652
    mo["candidate_source"] = "HUNT main-text Table 2 and Supplementary Table 24"
    mo["published_ldsc_h2"] = np.nan
    mo["published_ldsc_se"] = np.nan
    mo["published_ldsc_z"] = np.nan
    mo["published_ldsc_p"] = np.nan
    return pd.concat([sp, mo], ignore_index=True, sort=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "config/v2_candidate_exposure_inventory.tsv",
    )
    args = parser.parse_args()

    df = pd.concat([build_swedish(), build_hunt()], ignore_index=True, sort=False)
    df["species_equivalence_key"] = np.where(
        df["trait_family"].eq("species"),
        df["trait_name"].map(normalize_species)
        + "::"
        + df["phenotype_class"].astype(str),
        "",
    )
    df["module_equivalence_key"] = np.where(
        df["trait_family"].eq("module"),
        df["cohort"].str.lower() + ":" + df["trait_id"].astype(str),
        "",
    )
    df["gwas_url"] = df["accession"].map(gwas_url)
    df["local_path"] = [
        str(local_path(cohort, acc, phenotype_class))
        for cohort, acc, phenotype_class in zip(
            df["cohort"], df["accession"], df["phenotype_class"]
        )
    ]
    df["raw_file_present"] = [Path(p).is_file() for p in df["local_path"]]
    df["candidate_rule_pass"] = (
        pd.to_numeric(df["candidate_gws_count"], errors="coerce").fillna(0).gt(0)
        & pd.to_numeric(df["candidate_lead_p"], errors="coerce").lt(5e-8)
    )

    columns = [
        "cohort",
        "accession",
        "trait_family",
        "trait_id",
        "trait_name",
        "abundance_model",
        "phenotype_class",
        "sample_size",
        "candidate_gws_count",
        "candidate_lead_rsid",
        "candidate_lead_p",
        "candidate_rule_pass",
        "published_ldsc_h2",
        "published_ldsc_se",
        "published_ldsc_z",
        "published_ldsc_p",
        "species_equivalence_key",
        "module_equivalence_key",
        "candidate_source",
        "raw_file_present",
        "local_path",
        "gwas_url",
    ]
    df = df[columns].sort_values(["trait_family", "trait_name", "cohort", "accession"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False, na_rep="")

    summary = (
        df.groupby(["cohort", "trait_family"], dropna=False)
        .agg(candidates=("accession", "size"), raw_present=("raw_file_present", "sum"))
        .reset_index()
    )
    print(summary.to_string(index=False))
    print(f"Total candidates: {len(df)}")
    print(f"Candidate-rule failures: {(~df['candidate_rule_pass']).sum()}")
    overlaps = (
        df[df["trait_family"].eq("species")]
        .groupby("species_equivalence_key")["cohort"]
        .nunique()
    )
    print(f"Exact Swedish/HUNT species overlaps: {(overlaps > 1).sum()}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
