#!/usr/bin/env python3
"""Repeated exact rarefaction sensitivity analysis for Reviewer 1.2.

The script rarefies the bacterial species-leaf Bracken counts without
replacement, independently for every sample, and recomputes two participant-
level quantities for all 60 ciprofloxacin-exposed participants:

1. mean A1--A4 Shannon diversity; and
2. actual-day-normalized B6--C6 AUC of Bray--Curtis displacement from the
   participant-specific A1--A4 compositional centroid (D5).

The first replicate uses seed 20261009 and is required to reproduce the
audited ``rarefied_primary`` table produced by ``core_reanalysis.py``.  Later
replicates use consecutive seeds.  Rarefaction is applied to the complete
1,040-sample matrix before the exposed samples are selected, matching the
Methods-faithful branch used elsewhere in the response package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "raw"
OUT_DIR = ROOT / "reviewer_response" / "R1.2"
DEPTH = 2_693_000
BASE_SEED = 20_261_009
BASELINE = ("A1", "A2", "A3", "A4")
RECOVERY = ("B6", "C1", "C2", "C3", "C4", "C5", "C6")

# Read-only objects inherited by forked workers.  The host is Linux; refusing
# non-fork execution prevents silent multi-gigabyte copies on spawn platforms.
COUNTS: np.ndarray | None = None
EXPOSED_COLUMN_INDEX: np.ndarray | None = None
SUBJECT_LAYOUT: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def species_leaf_mask(index: pd.Index) -> np.ndarray:
    series = index.to_series(index=index)
    terminal = series.str.rsplit("|", n=1).str[-1]
    return (series.str.startswith("k__Bacteria|") & terminal.str.startswith("s__")).to_numpy()


def normalized_auc(day: np.ndarray, value: np.ndarray) -> float:
    order = np.argsort(day, kind="mergesort")
    x = np.asarray(day, dtype=float)[order]
    y = np.asarray(value, dtype=float)[order]
    span = float(x[-1] - x[0])
    if len(x) < 2 or not np.isfinite(span) or span <= 0:
        return float("nan")
    return float(np.trapz(y, x=x) / span)


def shannon_from_equal_depth_counts(counts: np.ndarray, depth: int) -> np.ndarray:
    """Return Shannon entropy per column without allocating a float matrix."""
    out = np.empty(counts.shape[1], dtype=np.float64)
    log_depth = float(np.log(depth))
    for column in range(counts.shape[1]):
        values = counts[:, column]
        positive = values[values > 0].astype(np.float64, copy=False)
        out[column] = log_depth - float(np.dot(positive, np.log(positive)) / depth)
    return out


def run_replicate(replicate: int) -> pd.DataFrame:
    if COUNTS is None or EXPOSED_COLUMN_INDEX is None or SUBJECT_LAYOUT is None:
        raise RuntimeError("Worker data were not initialized")
    seed = BASE_SEED + replicate - 1
    rng = np.random.default_rng(seed)
    n_features = COUNTS.shape[0]
    n_exposed = len(EXPOSED_COLUMN_INDEX)
    rarefied_exposed = np.empty((n_features, n_exposed), dtype=np.int32)
    exposed_lookup = {int(source): target for target, source in enumerate(EXPOSED_COLUMN_INDEX)}

    # Iterate over all 1,040 samples so replicate 1 is byte-for-byte consistent
    # with the audited global-seed branch before selecting exposed participants.
    for column in range(COUNTS.shape[1]):
        draw = rng.multivariate_hypergeometric(COUNTS[:, column], DEPTH)
        target = exposed_lookup.get(column)
        if target is not None:
            rarefied_exposed[:, target] = draw

    shannon = shannon_from_equal_depth_counts(rarefied_exposed, DEPTH)
    rows: list[dict[str, object]] = []
    for subject, baseline_idx, recovery_idx, recovery_days in SUBJECT_LAYOUT:
        baseline_centroid = rarefied_exposed[:, baseline_idx].mean(axis=1, dtype=np.float64)
        bray = np.empty(len(recovery_idx), dtype=np.float64)
        for position, sample_idx in enumerate(recovery_idx):
            vector = rarefied_exposed[:, sample_idx]
            # Both the sample vector and the averaged baseline centroid sum to
            # DEPTH, so the Bray--Curtis denominator is exactly 2 * DEPTH.
            bray[position] = np.abs(vector - baseline_centroid).sum() / (2.0 * DEPTH)
        rows.append(
            {
                "replicate": replicate,
                "seed": seed,
                "subject": subject,
                "baseline_shannon": float(shannon[baseline_idx].mean()),
                "recovery_auc_bray_actual_day": normalized_auc(recovery_days, bray),
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--processes", type=int, default=min(16, os.cpu_count() or 1))
    parser.add_argument(
        "--output",
        type=Path,
        default=OUT_DIR / "cache" / "R1_2_repeated_rarefaction_subject_metrics.csv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=OUT_DIR / "logs" / "R1_2_repeated_rarefaction_manifest.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.replicates < 2:
        raise ValueError("At least two rarefaction replicates are required")
    if args.processes < 1:
        raise ValueError("--processes must be positive")
    if "fork" not in mp.get_all_start_methods():
        raise RuntimeError("This memory-safe parallel implementation requires the fork start method")

    started = time.time()
    meta_path = RAW_DIR / "meta.csv"
    counts_path = RAW_DIR / "combine_brk.mpa"
    meta = pd.read_csv(meta_path, dtype={"ids": str, "subject": str, "canon": str})
    mpa = pd.read_csv(counts_path, sep="\t", index_col=0)
    mpa.columns = mpa.columns.str.removesuffix(".brk.report")
    if set(meta["ids"]) != set(mpa.columns):
        raise ValueError("raw metadata and Bracken sample columns differ")
    sample_order = meta["ids"].tolist()
    species = mpa.loc[species_leaf_mask(mpa.index), sample_order].astype(np.int64)
    if species.shape != (11_773, 1_040):
        raise ValueError(f"Unexpected bacterial species-leaf matrix shape: {species.shape}")
    totals = species.sum(axis=0).to_numpy(dtype=np.int64)
    if int(totals.min()) < DEPTH:
        raise ValueError("At least one sample is below the prespecified rarefaction depth")

    exposed_meta = meta.loc[meta["abx"] == "Y"].copy()
    exposed_meta["day"] = pd.to_numeric(exposed_meta["day"], errors="raise")
    if exposed_meta["subject"].nunique() != 60 or len(exposed_meta) != 960:
        raise ValueError("Expected 60 exposed participants and 960 exposed samples")
    if set(exposed_meta.groupby("subject").size()) != {16}:
        raise ValueError("Each exposed participant must have exactly 16 samples")

    column_lookup = {sample: position for position, sample in enumerate(sample_order)}
    exposed_columns = np.array([column_lookup[sample] for sample in exposed_meta["ids"]], dtype=np.int64)
    exposed_position = {sample: position for position, sample in enumerate(exposed_meta["ids"])}
    layout: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []
    for subject, rows in exposed_meta.groupby("subject", sort=True):
        baseline_rows = rows.loc[rows["canon"].isin(BASELINE)].sort_values("canon")
        recovery_rows = rows.loc[rows["canon"].isin(RECOVERY)].sort_values("day")
        if len(baseline_rows) != 4 or len(recovery_rows) != 7:
            raise ValueError(f"Incomplete baseline/recovery grid for {subject}")
        layout.append(
            (
                str(subject),
                np.array([exposed_position[x] for x in baseline_rows["ids"]], dtype=np.int64),
                np.array([exposed_position[x] for x in recovery_rows["ids"]], dtype=np.int64),
                recovery_rows["day"].to_numpy(dtype=np.float64),
            )
        )

    global COUNTS, EXPOSED_COLUMN_INDEX, SUBJECT_LAYOUT
    COUNTS = species.to_numpy(dtype=np.int64, copy=False)
    EXPOSED_COLUMN_INDEX = exposed_columns
    SUBJECT_LAYOUT = layout
    del mpa, species

    processes = min(args.processes, args.replicates)
    context = mp.get_context("fork")
    if processes == 1:
        frames = [run_replicate(rep) for rep in range(1, args.replicates + 1)]
    else:
        with context.Pool(processes=processes) as pool:
            frames = pool.map(run_replicate, range(1, args.replicates + 1), chunksize=1)
    result = pd.concat(frames, ignore_index=True)
    expected_rows = args.replicates * 60
    if len(result) != expected_rows or result[["baseline_shannon", "recovery_auc_bray_actual_day"]].isna().any().any():
        raise RuntimeError("Repeated-rarefaction output failed completeness checks")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    result.to_csv(temporary, index=False)
    temporary.replace(args.output)

    manifest = {
        "analysis": "R1.2 repeated exact rarefaction sensitivity",
        "created_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "replicates": args.replicates,
        "processes": processes,
        "seed_first": BASE_SEED,
        "seed_last": BASE_SEED + args.replicates - 1,
        "rarefaction_depth": DEPTH,
        "algorithm": "numpy.random.Generator.multivariate_hypergeometric; without replacement",
        "scope": "all 1,040 samples rarefied before retaining 60 exposed participants",
        "features": int(COUNTS.shape[0]),
        "samples_total": int(COUNTS.shape[1]),
        "samples_exposed": int(len(EXPOSED_COLUMN_INDEX)),
        "participants_exposed": len(SUBJECT_LAYOUT),
        "elapsed_seconds": round(time.time() - started, 3),
        "input_sha256": {
            str(counts_path.relative_to(ROOT)): sha256(counts_path),
            str(meta_path.relative_to(ROOT)): sha256(meta_path),
        },
        "output": display_path(args.output),
        "output_sha256": sha256(args.output),
        "software": {
            "python": sys.version,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
        },
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest_tmp = args.manifest.with_suffix(args.manifest.suffix + ".tmp")
    manifest_tmp.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    manifest_tmp.replace(args.manifest)
    print(json.dumps({"rows": len(result), "output": str(args.output), "elapsed_seconds": manifest["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
