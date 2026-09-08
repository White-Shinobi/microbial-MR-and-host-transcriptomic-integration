#!/usr/bin/env python3
"""Audit whether GWAS-SSF coordinates agree with coordinates encoded in variant_id."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "config/v2_candidate_exposure_inventory.tsv"
VARIANT_RE = re.compile(r"^(?:chr)?([^_]+)_(\d+)_([ACGT]+)_([ACGT]+)$", re.I)


def audit_one(path: Path, max_rows: int) -> dict[str, object]:
    compared = 0
    chromosome_match = 0
    position_match = 0
    liftover_rows = 0
    harmonised_code_rows = 0
    absolute_differences: list[int] = []
    examples: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        needed = {"chromosome", "base_pair_location", "variant_id"}
        if not needed.issubset(index):
            return {"rows_compared": 0, "audit_status": "required_columns_missing"}
        for line in handle:
            values = line.rstrip("\n").split("\t")
            match = VARIANT_RE.match(values[index["variant_id"]])
            if not match:
                continue
            try:
                reported_pos = int(values[index["base_pair_location"]])
                encoded_pos = int(match.group(2))
            except ValueError:
                continue
            reported_chr = values[index["chromosome"]].removeprefix("chr")
            encoded_chr = match.group(1).removeprefix("chr")
            conversion = (
                values[index["hm_coordinate_conversion"]]
                if "hm_coordinate_conversion" in index
                else ""
            )
            try:
                hm_code = int(values[index["hm_code"]]) if "hm_code" in index else 0
            except ValueError:
                hm_code = 0
            compared += 1
            liftover_rows += int(conversion == "lo")
            harmonised_code_rows += int(1 <= hm_code <= 13)
            chr_ok = reported_chr == encoded_chr
            pos_ok = reported_pos == encoded_pos
            chromosome_match += int(chr_ok)
            position_match += int(pos_ok)
            absolute_differences.append(abs(reported_pos - encoded_pos))
            if len(examples) < 3 and (not chr_ok or not pos_ok):
                examples.append(
                    f"{reported_chr}:{reported_pos}!={encoded_chr}:{encoded_pos}"
                )
            if compared >= max_rows:
                break
    if not compared:
        return {"rows_compared": 0, "audit_status": "no_parseable_variant_ids"}
    mismatch_rate = 1 - position_match / compared
    liftover_fraction = liftover_rows / compared
    harmonised_code_fraction = harmonised_code_rows / compared
    if mismatch_rate <= 0.01:
        status = "coordinates_agree"
    elif liftover_fraction >= 0.95 and harmonised_code_fraction >= 0.95:
        status = "expected_liftover_difference"
    else:
        status = "coordinate_mismatch_requires_review"
    return {
        "rows_compared": compared,
        "chromosome_match_fraction": chromosome_match / compared,
        "position_match_fraction": position_match / compared,
        "position_mismatch_fraction": mismatch_rate,
        "liftover_fraction": liftover_fraction,
        "harmonised_hm_code_fraction": harmonised_code_fraction,
        "median_absolute_position_difference": pd.Series(absolute_differences).median(),
        "mismatch_examples": ";".join(examples),
        "audit_status": status,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rows", type=int, default=10000)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/data_qc/gwas_coordinate_audit.tsv",
    )
    args = parser.parse_args()
    inventory = pd.read_csv(INVENTORY, sep="\t")
    rows = []
    for row in inventory.itertuples(index=False):
        path = Path(row.local_path)
        if not path.is_file():
            result = {"rows_compared": 0, "audit_status": "raw_missing"}
        else:
            result = audit_one(path, args.max_rows)
        rows.append(
            {
                "cohort": row.cohort,
                "accession": row.accession,
                "trait_family": row.trait_family,
                "trait_name": row.trait_name,
                "path": str(path),
                **result,
            }
        )
    output = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, sep="\t", index=False, na_rep="")
    print(
        output.groupby(["cohort", "audit_status"]).size().to_string()
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
