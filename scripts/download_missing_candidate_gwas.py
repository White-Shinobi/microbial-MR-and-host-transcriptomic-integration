#!/usr/bin/env python3
"""Download only missing V2 candidate GWAS files; never overwrite shared V1 data."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import os
import re
import subprocess
import time
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "config/v2_candidate_exposure_inventory.tsv"


def md5sum(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def expected_md5(meta_path: Path) -> str | None:
    if not meta_path.is_file():
        return None
    match = re.search(
        r"^data_file_md5sum:\s*([0-9a-fA-F]{32})\s*$",
        meta_path.read_text(errors="replace"),
        flags=re.MULTILINE,
    )
    return match.group(1).lower() if match else None


def fetch_small(url: str, destination: Path, required: bool = True) -> None:
    if destination.exists():
        return
    tmp = destination.with_name(destination.name + ".part")
    try:
        with urllib.request.urlopen(url, timeout=60) as response, tmp.open("wb") as out:
            out.write(response.read())
    except Exception:
        tmp.unlink(missing_ok=True)
        if required:
            raise
        return
    os.replace(tmp, destination)


def download_one(row: dict[str, object]) -> dict[str, object]:
    acc = str(row["accession"])
    path = Path(str(row["local_path"]))
    if not path.is_absolute():
        path = ROOT / path
    url = str(row["gwas_url"])
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = path.parent / f"{acc}.h.tsv.gz-meta.yaml"
    checksum_file = path.parent / f"{acc}.harmonised.md5sum.txt"
    base = url.rsplit("/", 1)[0]
    fetch_small(f"{base}/{meta.name}", meta)
    fetch_small(f"{base}/{checksum_file.name}", checksum_file, required=False)
    want = expected_md5(meta)

    if path.is_file():
        return {
            "accession": acc,
            "status": "reused_existing",
            "path": str(path),
            "expected_md5": want,
            "observed_md5": "",
            "seconds": 0,
        }

    tmp = path.with_name(path.name + ".part")
    start = time.time()
    command = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "5",
        "--retry-delay",
        "2",
        "--silent",
        "--show-error",
        "--continue-at",
        "-",
        "--output",
        str(tmp),
        url,
    ]
    subprocess.run(command, check=True)
    got = md5sum(tmp)
    if want and got != want:
        raise RuntimeError(f"{acc}: MD5 mismatch, expected {want}, observed {got}")
    os.replace(tmp, path)
    return {
        "accession": acc,
        "status": "downloaded",
        "path": str(path),
        "expected_md5": want,
        "observed_md5": got,
        "seconds": round(time.time() - start, 2),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--log",
        type=Path,
        default=ROOT / "logs/candidate_gwas_download.tsv",
    )
    args = parser.parse_args()

    inventory = pd.read_csv(INVENTORY, sep="\t")
    rows = inventory.to_dict("records")
    results: list[dict[str, object]] = []
    failures: list[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download_one, row): row for row in rows}
        for future in concurrent.futures.as_completed(futures):
            acc = str(futures[future]["accession"])
            try:
                result = future.result()
                results.append(result)
                print(
                    f"{result['status']}\t{acc}\t{result['seconds']}s",
                    flush=True,
                )
            except Exception as exc:
                failures.append(f"{acc}: {exc}")
                print(f"FAILED\t{acc}\t{exc}", flush=True)

    args.log.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).sort_values("accession").to_csv(
        args.log, sep="\t", index=False
    )
    print(f"Completed: {len(results)}; failures: {len(failures)}")
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
