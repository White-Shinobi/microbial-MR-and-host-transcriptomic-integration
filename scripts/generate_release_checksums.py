#!/usr/bin/env python3
"""Generate SHA-256 checksums for the public repository contents."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "SHA256SUMS.txt"
SKIP_PARTS = {".git", ".venv", "venv", "__pycache__"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    files = [
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path != OUTPUT
        and not any(part in SKIP_PARTS for part in path.parts)
    ]
    lines = [f"{sha256(path)}  {path.relative_to(ROOT)}" for path in sorted(files)]
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.name} for {len(files)} files.")


if __name__ == "__main__":
    main()
