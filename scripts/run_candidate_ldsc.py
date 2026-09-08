#!/usr/bin/env python3
"""Munge all V2 candidate GWAS and estimate SNP heritability with LDSC."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import re
import sys
import subprocess
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path(os.environ.get("LDSC_PYTHON", sys.executable))
LDSC_HOME = Path(os.environ.get("LDSC_HOME", ROOT / "tools" / "ldsc"))
LDSC = LDSC_HOME / "ldsc.py"
MUNGE = LDSC_HOME / "munge_sumstats.py"
PYDEPS = os.environ.get("LDSC_PYDEPS", "")
LDREF = Path(os.environ.get(
    "LDSC_LD_REFERENCE", ROOT / "data" / "reference" / "ldsc" / "eur_w_ld_chr"
))
INVENTORY = ROOT / "config/v2_candidate_exposure_inventory.tsv"
RESULT_DIR = ROOT / "results/ldsc"
LOG_DIR = RESULT_DIR / "logs"


def parse_h2_log(path: Path) -> dict[str, object] | None:
    if not path.is_file():
        return None
    text = path.read_text(errors="replace")
    h2 = re.search(
        r"Total Observed scale h2:\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
        r"\s*\(([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\)",
        text,
    )
    intercept = re.search(
        r"Intercept:\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
        r"\s*\(([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\)",
        text,
    )
    mean_chi = re.search(
        r"Mean Chi\^2:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
        text,
    )
    n_snps = re.search(r"Read summary statistics for\s+(\d+)\s+SNPs", text)
    n_merged = re.search(
        r"After merging with regression SNP LD,\s+(\d+)\s+SNPs remain", text
    )
    if not h2:
        return None
    estimate = float(h2.group(1))
    se = float(h2.group(2))
    return {
        "ldsc_h2": estimate,
        "ldsc_h2_se": se,
        "ldsc_h2_z": estimate / se if se > 0 else float("nan"),
        "ldsc_intercept": float(intercept.group(1)) if intercept else float("nan"),
        "ldsc_intercept_se": float(intercept.group(2)) if intercept else float("nan"),
        "ldsc_mean_chisq": float(mean_chi.group(1)) if mean_chi else float("nan"),
        "ldsc_n_sumstats": int(n_snps.group(1)) if n_snps else None,
        "ldsc_n_merged": int(n_merged.group(1)) if n_merged else None,
    }


def run_command(command: list[str], stdout_path: Path) -> None:
    env = os.environ.copy()
    if PYDEPS:
        env["PYTHONPATH"] = PYDEPS
    with stdout_path.open("w") as stdout:
        subprocess.run(
            [str(PYTHON), *map(str, command)],
            cwd=ROOT,
            env=env,
            stdout=stdout,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_one(row: dict[str, object]) -> dict[str, object]:
    acc = str(row["accession"])
    raw = Path(str(row["local_path"]))
    if not raw.is_absolute():
        raw = ROOT / raw
    out_prefix = RESULT_DIR / acc
    h2_log = out_prefix.with_suffix(".log")
    existing = parse_h2_log(h2_log)
    if existing:
        return {
            "accession": acc,
            "status": "reused_completed_ldsc",
            "seconds": 0,
            **existing,
        }
    if not raw.is_file():
        return {"accession": acc, "status": "raw_missing", "seconds": 0}

    start = time.time()
    munged_dir = RESULT_DIR / "munged" / acc
    munged_dir.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    munged_prefix = munged_dir / acc
    munged = munged_prefix.with_suffix(".sumstats.gz")

    if not munged.is_file():
        run_command(
            [
                MUNGE,
                "--sumstats",
                raw,
                "--snp",
                "rsid",
                "--a1",
                "effect_allele",
                "--a2",
                "other_allele",
                "--p",
                "p_value",
                "--signed-sumstats",
                "beta,0",
                "--frq",
                "effect_allele_frequency",
                "--N",
                str(int(row["sample_size"])),
                "--merge-alleles",
                LDREF / "w_hm3.snplist",
                "--chunksize",
                "500000",
                "--out",
                munged_prefix,
            ],
            LOG_DIR / f"{acc}.munge.stdout.log",
        )
    run_command(
        [
            LDSC,
            "--h2",
            munged,
            "--ref-ld-chr",
            f"{LDREF}/",
            "--w-ld-chr",
            f"{LDREF}/",
            "--out",
            out_prefix,
        ],
        LOG_DIR / f"{acc}.h2.stdout.log",
    )
    parsed = parse_h2_log(h2_log)
    if not parsed:
        raise RuntimeError(f"{acc}: LDSC completed without a parseable h2 estimate")
    return {
        "accession": acc,
        "status": "completed_ldsc",
        "seconds": round(time.time() - start, 2),
        **parsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--out", type=Path, default=RESULT_DIR / "candidate_ldsc_results.tsv"
    )
    args = parser.parse_args()
    inventory = pd.read_csv(INVENTORY, sep="\t")
    rows = inventory.to_dict("records")
    results: list[dict[str, object]] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, row): row for row in rows}
        for future in concurrent.futures.as_completed(futures):
            acc = str(futures[future]["accession"])
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "accession": acc,
                    "status": "failed",
                    "seconds": 0,
                    "error": str(exc),
                }
            results.append(result)
            print(
                f"{result['status']}\t{acc}\t{result.get('seconds', 0)}s",
                flush=True,
            )

    results_df = pd.DataFrame(results)
    merged = inventory.merge(results_df, on="accession", how="left", validate="one_to_one")
    merged["ldsc_gate_pass"] = (
        pd.to_numeric(merged["ldsc_h2"], errors="coerce").gt(0)
        & pd.to_numeric(merged["ldsc_h2"], errors="coerce").lt(1)
        & pd.to_numeric(merged["ldsc_h2_z"], errors="coerce").gt(1.64)
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False, na_rep="")
    print(merged.groupby(["cohort", "trait_family", "status"]).size().to_string())
    print(f"Wrote {args.out}")
    if merged["status"].eq("failed").any():
        raise SystemExit("One or more LDSC jobs failed; see result table and logs.")


if __name__ == "__main__":
    main()
