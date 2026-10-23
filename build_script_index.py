#!/usr/bin/env python3
"""Build a deterministic inventory for the curated analysis-code package."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "SCRIPT_INDEX.tsv"
EXCLUDE = {"assemble_code_package.sh", "build_script_index.py"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


scripts = sorted(
    path
    for path in ROOT.rglob("*")
    if path.is_file()
    and path.name not in EXCLUDE
    and path.suffix.lower() in {".py", ".r", ".sh", ".txt"}
)

with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["category", "script", "language", "bytes", "sha256"])
    for path in scripts:
        rel = path.relative_to(ROOT)
        writer.writerow(
            [
                rel.parts[0],
                rel.as_posix(),
                path.suffix.lower().lstrip("."),
                path.stat().st_size,
                sha256(path),
            ]
        )

print(f"Indexed {len(scripts)} scripts in {OUTPUT}")

