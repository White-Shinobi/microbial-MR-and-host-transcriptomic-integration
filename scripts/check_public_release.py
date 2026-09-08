#!/usr/bin/env python3
"""Fail if a public release contains likely local paths, emails or secrets."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKIP_PARTS = {".git", ".venv", "venv", "__pycache__"}
TEXT_SUFFIXES = {
    "", ".R", ".r", ".py", ".sh", ".awk", ".md", ".txt", ".tsv",
    ".csv", ".yml", ".yaml", ".json", ".env", ".cff", ".gitignore",
}
PATTERNS = {
    "macOS user path": re.compile("/" + "Users" + r"/[^/\s]+/"),
    "Linux home path": re.compile("/" + "home" + r"/[^/\s]+/"),
    "Windows user path": re.compile(r"[A-Za-z]:[\\/]Users[\\/][^\\/\s]+[\\/]"),
    "email address": re.compile(r"(?<![\w.-])[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "credential assignment": re.compile(
        r"(?i)(?:api[_-]?key|access[_-]?token|client[_-]?secret|password)\s*[:=]\s*[\"'][^\"']+[\"']"
    ),
}


def main() -> None:
    findings: list[str] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or any(part in SKIP_PARTS for part in path.parts):
            continue
        if path.suffix not in TEXT_SUFFIXES and path.name not in {"LICENSE", ".gitignore"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for label, pattern in PATTERNS.items():
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"{path.relative_to(ROOT)}:{line}: {label}")
    if findings:
        raise SystemExit("\n".join(findings))
    print("Privacy scan passed: no local user paths, emails, private keys, or inline credentials found.")


if __name__ == "__main__":
    main()
