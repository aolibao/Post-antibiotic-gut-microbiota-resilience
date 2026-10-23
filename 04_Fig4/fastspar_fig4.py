#!/usr/bin/env python3
"""Recompute the twelve Figure 4 genus networks with FastSpar.

The manuscript provides Bracken counts but not the original FastSpar inputs,
feature filter, bootstrap seeds, or result matrices.  This script therefore
implements a deterministic, explicitly auditable reconstruction:

* bacterial terminal ``g__`` rows are read directly from
  ``raw/combine_brk.mpa``;
* one fixed set of the 275 most abundant genera across all CIP samples is used
  for every network (rather than selecting a different feature universe in
  each panel);
* HDR/LDR membership follows the endpoint ranking in the supplied
  ``species.csv`` table, which is the grouping used by the source-faithful
  Figure 2/4 reconstruction;
* Baseline and End pool their original sample columns (60 and 45 columns per
  group); B6, C1, C2, and C3 each contain 15 columns per group;
* FastSpar correlation, 500 bootstraps (configurable), exact bootstrap P
  values, and within-network Benjamini-Hochberg Q values are exported;
* both the caption rule (|rho| > 0.6 and raw P < 0.05) and the Methods rule
  (|rho| > 0.6 and BH Q < 0.05) are retained as separate edge lists.

This is a method-faithful rerun, not a claim that the unavailable original
filtering and seed choices have been recovered exactly.
"""

from __future__ import annotations

import argparse
import concurrent.futures
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


TIME_ORDER = [f"A{i}" for i in range(1, 5)] + [f"B{i}" for i in range(1, 7)] + [f"C{i}" for i in range(1, 7)]
BASELINE = ["A1", "A2", "A3", "A4"]
ENDPOINT = ["C4", "C5", "C6"]
NETWORK_STAGES: dict[str, list[str]] = {
    "Baseline": BASELINE,
    "B6": ["B6"],
    "C1": ["C1"],
    "C2": ["C2"],
    "C3": ["C3"],
    "End": ENDPOINT,
}
GROUP_ORDER = ["HDR", "LDR"]
DEFAULT_SEED = 20251010
DEFAULT_FEATURES = 275
DEFAULT_BOOTSTRAPS = 500
DEFAULT_RHO_THRESHOLD = 0.6
DEFAULT_ALPHA = 0.05


@dataclass(frozen=True)
class CompletedFastSparResults:
    """Validated final outputs consumed by the Figure 4 renderer."""

    root: Path
    manifest: dict[str, object]
    summary: pd.DataFrame
    caption_edges: pd.DataFrame
    methods_edges: pd.DataFrame
    panel_summaries: dict[str, dict[str, object]]


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjustment with NaN preservation."""
    raw = np.asarray(list(values), dtype=float)
    adjusted = np.full(raw.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(raw))
    if not len(finite):
        return adjusted
    p_values = raw[finite]
    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order]
    scale = len(ranked) / np.arange(1, len(ranked) + 1, dtype=float)
    corrected = np.minimum.accumulate((ranked * scale)[::-1])[::-1]
    restored = np.empty_like(corrected)
    restored[order] = np.minimum(corrected, 1.0)
    adjusted[finite] = restored
    return adjusted


def terminal_rank(path: str) -> str:
    return str(path).rsplit("|", 1)[-1]


def terminal_genus(path: str) -> str:
    terminal = terminal_rank(path)
    if not terminal.startswith("g__"):
        raise ValueError(f"not a terminal genus path: {path}")
    return terminal.removeprefix("g__").replace("_", " ").strip()


def load_metadata(root: Path) -> pd.DataFrame:
    path = root / "raw" / "meta.csv"
    meta = pd.read_csv(path, dtype={"ids": str, "subject": str, "canon": str, "abx": str})
    required = {"ids", "subject", "canon", "abx"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {sorted(missing)}")
    if meta["ids"].duplicated().any():
        raise ValueError(f"{path} contains duplicate sample IDs")
    exposed = meta.loc[meta["abx"].eq("Y")]
    if exposed["subject"].nunique() != 60 or len(exposed) != 960:
        raise ValueError(
            "expected 60 CIP subjects x 16 samples (960 rows), found "
            f"{exposed['subject'].nunique()} subjects and {len(exposed)} rows"
        )
    observed_times = set(exposed["canon"].dropna())
    if set(TIME_ORDER).difference(observed_times):
        raise ValueError(f"CIP metadata lacks canonical times: {sorted(set(TIME_ORDER).difference(observed_times))}")
    per_subject = exposed.groupby("subject")["canon"].nunique()
    if not per_subject.eq(16).all():
        bad = per_subject.loc[~per_subject.eq(16)].to_dict()
        raise ValueError(f"CIP subjects do not each have 16 canonical times: {bad}")
    return meta


def legacy_group_assignments(root: Path, meta: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the endpoint ranking used by the source-faithful figures."""
    path = root / "species.csv"
    table = pd.read_csv(path, index_col=0)
    exposed = meta.loc[meta["abx"].eq("Y")].copy()
    missing = set(exposed["ids"]).difference(table.columns)
    if missing:
        raise ValueError(f"{path} lacks {len(missing)} exposed samples; examples: {sorted(missing)[:5]}")
    values = table.loc[:, exposed["ids"]].astype(float)
    proportions = values.div(values.sum(axis=0).replace(0.0, np.nan), axis=1)
    array = proportions.to_numpy(dtype=float)
    log_array = np.zeros_like(array)
    np.log(array, out=log_array, where=array > 0)
    shannon = pd.Series(-(array * log_array).sum(axis=0), index=proportions.columns)
    exposed["grouping_shannon"] = exposed["ids"].map(shannon)
    baseline = (
        exposed.loc[exposed["canon"].isin(BASELINE)]
        .groupby("subject")["grouping_shannon"]
        .mean()
    )
    endpoint = (
        exposed.loc[exposed["canon"].isin(ENDPOINT)]
        .groupby("subject")["grouping_shannon"]
        .mean()
        .sort_values(kind="mergesort")
    )
    result = pd.DataFrame(
        {
            "subject": endpoint.index,
            "baseline_grouping_shannon": baseline.reindex(endpoint.index).to_numpy(),
            "endpoint_grouping_shannon": endpoint.to_numpy(),
            "endpoint_rank": np.arange(1, len(endpoint) + 1),
            "recovery_group": "Middle",
        }
    )
    result.loc[:14, "recovery_group"] = "LDR"
    result.loc[len(result) - 15 :, "recovery_group"] = "HDR"
    return result


def manifest_group_assignments(
    path: Path,
    meta: pd.DataFrame,
    definition_id: str,
) -> pd.DataFrame:
    """Load one fixed recovery definition from the reviewer-response manifest."""
    table = pd.read_csv(path, dtype={"subject": str, "definition_id": str})
    required = {"subject", "definition_id", "recovery_group"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path} lacks grouping columns: {sorted(missing)}")
    selected = table.loc[table["definition_id"].eq(definition_id)].copy()
    exposed_subjects = set(meta.loc[meta["abx"].eq("Y"), "subject"].astype(str))
    if len(selected) != 60 or selected["subject"].nunique() != 60:
        raise ValueError(
            f"{definition_id}: expected 60 unique assignment rows, found {len(selected)}"
        )
    if set(selected["subject"]) != exposed_subjects:
        raise ValueError(f"{definition_id}: grouping subjects differ from exposed metadata")
    counts = selected["recovery_group"].value_counts().to_dict()
    if counts != {"Middle": 30, "HDR": 15, "LDR": 15}:
        raise ValueError(f"{definition_id}: invalid recovery-group counts {counts}")
    rank_source = (
        "metric_rank_ascending"
        if "metric_rank_ascending" in selected.columns
        else "rank"
    )
    if rank_source not in selected.columns:
        raise ValueError(f"{path} lacks a stable rank column")
    result = pd.DataFrame(
        {
            "subject": selected["subject"].astype(str),
            "rank": pd.to_numeric(selected[rank_source], errors="raise").astype(int),
            "recovery_group": selected["recovery_group"].astype(str),
        }
    )
    for optional in (
        "metric_value",
        "metric_column",
        "definition_short",
        "definition_full",
        "high_recovery_direction",
        "dataset_branch",
    ):
        if optional in selected.columns:
            result[optional] = selected[optional].to_numpy()
    return result.sort_values(["rank", "subject"], kind="mergesort").reset_index(drop=True)


def load_direct_genus_counts(root: Path, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return direct bacterial g__ counts and a genus-to-source-path table."""
    path = root / "raw" / "combine_brk.mpa"
    sample_ids = meta["ids"].tolist()
    wanted = ["#Classification", *[f"{sample}.brk.report" for sample in sample_ids]]
    table = pd.read_csv(path, sep="\t", usecols=wanted)
    table = table.rename(
        columns={
            "#Classification": "taxonomy_path",
            **{
                column: column.removesuffix(".brk.report")
                for column in table.columns
                if column != "#Classification"
            },
        }
    )
    if set(sample_ids) != set(table.columns[1:]):
        raise ValueError("raw metadata IDs and Bracken sample columns differ")
    paths = table["taxonomy_path"].astype(str)
    selected = paths.str.startswith("k__Bacteria|") & paths.map(terminal_rank).str.startswith("g__")
    genus_rows = table.loc[selected, ["taxonomy_path", *sample_ids]].copy()
    genus_rows.insert(1, "genus", genus_rows["taxonomy_path"].map(terminal_genus))
    numeric = genus_rows.loc[:, sample_ids].apply(pd.to_numeric, errors="raise")
    if (numeric.to_numpy() < 0).any():
        raise ValueError("Bracken genus counts contain negative values")

    # Terminal genus names are unique in the supplied table.  Grouping is kept
    # as a defensive fallback and the full path mapping records any collapse.
    numeric.insert(0, "genus", genus_rows["genus"].to_numpy())
    counts = numeric.groupby("genus", sort=False)[sample_ids].sum()
    path_map = (
        genus_rows.groupby("genus", sort=False)["taxonomy_path"]
        .agg(lambda values: " || ".join(dict.fromkeys(map(str, values))))
        .rename("taxonomy_path")
        .to_frame()
    )
    counts.index.name = "genus"
    return counts.astype(np.int64), path_map


def select_global_features(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    path_map: pd.DataFrame,
    n_features: int,
) -> tuple[list[str], pd.DataFrame]:
    exposed_ids = meta.loc[meta["abx"].eq("Y"), "ids"].tolist()
    exposed = counts.loc[:, exposed_ids]
    totals = exposed.sum(axis=0)
    if (totals <= 0).any():
        raise ValueError("one or more CIP genus-count columns have zero total")
    relative = exposed.div(totals, axis=1)
    ranking = pd.DataFrame(
        {
            "genus": exposed.index,
            "global_mean_relative_abundance": relative.mean(axis=1).to_numpy(),
            "global_prevalence": (exposed > 0).mean(axis=1).to_numpy(),
            "global_total_count": exposed.sum(axis=1).to_numpy(dtype=np.int64),
        }
    )
    ranking = ranking.sort_values(
        ["global_mean_relative_abundance", "genus"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    if len(ranking) < n_features:
        raise ValueError(f"requested {n_features} genera but only {len(ranking)} are available")
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    ranking = ranking.merge(path_map.reset_index(), on="genus", how="left", validate="one_to_one")
    selected = ranking.head(n_features).copy()
    return selected["genus"].tolist(), selected


def panel_key(group: str, stage: str) -> str:
    return f"{group}__{stage}"


def build_panel_inputs(
    output: Path,
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    groups: pd.DataFrame,
    features: Sequence[str],
) -> pd.DataFrame:
    inputs = output / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    group_lookup = groups.set_index("subject")["recovery_group"]
    selected_meta = meta.loc[meta["abx"].eq("Y")].copy()
    selected_meta["recovery_group"] = selected_meta["subject"].map(group_lookup)
    rows: list[dict[str, object]] = []
    for group in GROUP_ORDER:
        for stage, canons in NETWORK_STAGES.items():
            panel_meta = selected_meta.loc[
                selected_meta["recovery_group"].eq(group)
                & selected_meta["canon"].isin(canons)
            ].copy()
            panel_meta["canon_order"] = panel_meta["canon"].map(TIME_ORDER.index)
            panel_meta = panel_meta.sort_values(["subject", "canon_order"], kind="mergesort")
            expected = 60 if stage == "Baseline" else 45 if stage == "End" else 15
            if len(panel_meta) != expected:
                raise ValueError(f"{group}/{stage}: expected {expected} samples, found {len(panel_meta)}")
            sample_ids = panel_meta["ids"].tolist()
            matrix = counts.loc[list(features), sample_ids].copy()
            zero_total = matrix.sum(axis=1).eq(0)
            if zero_total.any():
                dropped = matrix.index[zero_total].tolist()
                log(f"{group}/{stage}: dropping {len(dropped)} all-zero global features")
                matrix = matrix.loc[~zero_total]
            input_path = inputs / f"{panel_key(group, stage)}.tsv"
            matrix.index.name = "#OTU ID"
            matrix.to_csv(input_path, sep="\t", index=True)
            rows.append(
                {
                    "group": group,
                    "stage": stage,
                    "canons": ",".join(canons),
                    "n_subjects": int(panel_meta["subject"].nunique()),
                    "n_samples": int(len(sample_ids)),
                    "candidate_features": int(len(matrix)),
                    "dropped_all_zero_features": int(zero_total.sum()),
                    "input": str(input_path.resolve()),
                    "input_sha256": sha256_file(input_path),
                }
            )
    return pd.DataFrame(rows)


def resolve_binary(name: str, bin_dir: Path | None) -> Path:
    candidates: list[Path] = []
    if bin_dir is not None:
        candidates.append(bin_dir / name)
    found = shutil.which(name)
    if found:
        candidates.append(Path(found))
    candidates.append(Path("/home/qqt/.conda/envs/fastspar_env/bin") / name)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise FileNotFoundError(f"cannot find executable {name}; use --fastspar-bin-dir")


def run_command(command: Sequence[str], *, quiet: bool = False) -> subprocess.CompletedProcess[str]:
    normalized = [str(part) for part in command]
    # Each bootstrap fit is already parallelized explicitly by this driver.
    # Pin implicit numerical/OpenMP pools to one thread so a panel-level run
    # cannot multiply ``panel_jobs * jobs`` by the host BLAS thread count.
    command_env = os.environ.copy()
    for variable in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        command_env[variable] = "1"
    completed = subprocess.run(
        normalized,
        check=False,
        text=True,
        env=command_env,
        stdout=subprocess.DEVNULL if quiet else None,
        stderr=subprocess.PIPE if quiet else None,
    )
    if completed.returncode:
        detail = (completed.stderr or "").strip()
        if len(detail) > 4000:
            detail = detail[-4000:]
        suffix = f"\n{detail}" if detail else ""
        raise RuntimeError(
            f"command exited {completed.returncode}: {shlex.join(normalized)}{suffix}"
        )
    return completed


def natural_bootstrap_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.tsv$", path.name)
    if not match:
        raise ValueError(f"cannot parse bootstrap index from {path}")
    return int(match.group(1))


def read_square_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{path} is not square: {matrix.shape}")
    if set(matrix.index) != set(matrix.columns):
        raise ValueError(f"{path} row and column feature labels differ")
    return matrix.loc[matrix.index, matrix.index].astype(float)


def adjust_matrix(p_values: pd.DataFrame) -> pd.DataFrame:
    n_features = len(p_values)
    upper = np.triu_indices(n_features, 1)
    q_values = np.full((n_features, n_features), np.nan, dtype=float)
    adjusted = bh_adjust(p_values.to_numpy(dtype=float)[upper])
    q_values[upper] = adjusted
    q_values[(upper[1], upper[0])] = adjusted
    return pd.DataFrame(q_values, index=p_values.index, columns=p_values.columns)


def edge_table(
    correlation: pd.DataFrame,
    p_values: pd.DataFrame,
    q_values: pd.DataFrame,
    *,
    rho_threshold: float,
    alpha: float,
    significance: str,
) -> pd.DataFrame:
    names = correlation.index.to_numpy(dtype=str)
    rho = correlation.to_numpy(dtype=float)
    raw = p_values.loc[correlation.index, correlation.columns].to_numpy(dtype=float)
    q = q_values.loc[correlation.index, correlation.columns].to_numpy(dtype=float)
    upper = np.triu_indices(len(names), 1)
    metric = raw if significance == "raw_p" else q
    keep = (
        np.isfinite(rho[upper])
        & np.isfinite(metric[upper])
        & (np.abs(rho[upper]) > rho_threshold)
        & (metric[upper] < alpha)
    )
    first = upper[0][keep]
    second = upper[1][keep]
    result = pd.DataFrame(
        {
            "genus1": names[first],
            "genus2": names[second],
            "rho": rho[first, second],
            "abs_rho": np.abs(rho[first, second]),
            "p_raw": raw[first, second],
            "q_BH_within_network": q[first, second],
        }
    )
    result["sign"] = np.where(result["rho"] >= 0, "positive", "negative")
    return result.sort_values(
        ["abs_rho", "p_raw", "genus1", "genus2"],
        ascending=[False, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)


def node_count(edges: pd.DataFrame) -> int:
    if edges.empty:
        return 0
    return len(set(edges["genus1"]).union(edges["genus2"]))


def fit_bootstrap(
    fastspar: Path,
    bootstrap_input: Path,
    correlation_output: Path,
    *,
    iterations: int,
    exclusion_iterations: int,
    exclusion_threshold: float,
    threads: int,
    seed: int,
) -> tuple[int, float]:
    index = natural_bootstrap_index(bootstrap_input)
    started = time.perf_counter()
    command = [
        fastspar,
        "-c",
        bootstrap_input,
        "-r",
        correlation_output,
        "-a",
        os.devnull,
        "-i",
        str(iterations),
        "-x",
        str(exclusion_iterations),
        "-e",
        str(exclusion_threshold),
        "-t",
        str(threads),
        "-s",
        str(seed + index + 1),
        "-y",
    ]
    run_command(command, quiet=True)
    return index, time.perf_counter() - started


def safe_remove_tree(path: Path, required_parent: Path) -> None:
    resolved = path.resolve()
    parent = required_parent.resolve()
    if resolved == parent or parent not in resolved.parents:
        raise ValueError(f"refusing to remove unsafe path {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def run_panel(
    row: pd.Series,
    result_root: Path,
    work_root: Path,
    binaries: dict[str, Path],
    args: argparse.Namespace,
    bootstraps: int,
) -> dict[str, object]:
    group = str(row["group"])
    stage = str(row["stage"])
    key = panel_key(group, stage)
    input_path = Path(str(row["input"]))
    result_dir = result_root / key
    result_dir.mkdir(parents=True, exist_ok=True)
    panel_work = work_root / key

    correlation_path = result_dir / "correlation.tsv"
    covariance_path = result_dir / "covariance.tsv"
    pvalue_path = result_dir / "pvalues_raw.tsv"
    qvalue_path = result_dir / "qvalues_BH.tsv"
    raw_edges_path = result_dir / "edges_caption_raw_p.tsv"
    q_edges_path = result_dir / "edges_methods_BH_q.tsv"
    summary_path = result_dir / "summary.json"
    expected_outputs = [
        correlation_path,
        covariance_path,
        pvalue_path,
        qvalue_path,
        raw_edges_path,
        q_edges_path,
    ]

    if summary_path.exists() and all(path.exists() for path in expected_outputs) and not args.force:
        previous = json.loads(summary_path.read_text(encoding="utf-8"))
        if (
            int(previous.get("bootstraps", -1)) == bootstraps
            and int(previous.get("seed", -1)) == args.seed
            and int(previous.get("iterations", -1)) == args.iterations
            and int(previous.get("exclusion_iterations", -1)) == args.exclusion_iterations
            and float(previous.get("exclusion_threshold", math.nan)) == args.exclusion_threshold
            and float(previous.get("rho_threshold_strict_gt", math.nan)) == args.rho_threshold
            and float(previous.get("alpha_strict_lt", math.nan)) == args.alpha
            and previous.get("input_sha256") == row["input_sha256"]
        ):
            log(f"{key}: using completed outputs")
            return previous

    # A non-matching completed run or an interrupted run must not leak stale
    # bootstrap tables into this run.  Only this panel's scoped work directory
    # is removed; final outputs stay in place and are atomically overwritten by
    # the FastSpar programs.
    safe_remove_tree(panel_work, work_root)
    panel_work.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    log(f"{key}: fitting original FastSpar matrix")
    run_command(
        [
            binaries["fastspar"],
            "-c",
            input_path,
            "-r",
            correlation_path,
            "-a",
            covariance_path,
            "-i",
            str(args.iterations),
            "-x",
            str(args.exclusion_iterations),
            "-e",
            str(args.exclusion_threshold),
            "-t",
            str(args.threads_per_fit),
            "-s",
            str(args.seed),
            "-y",
        ]
    )

    bootstrap_dir = panel_work / "bootstrap_tables"
    correlation_dir = panel_work / "bootstrap_correlations"
    bootstrap_dir.mkdir(parents=True, exist_ok=True)
    correlation_dir.mkdir(parents=True, exist_ok=True)
    bootstrap_prefix = bootstrap_dir / "bootstrap"
    correlation_prefix = correlation_dir / "correlation_"
    log(f"{key}: generating {bootstraps} bootstrap tables")
    run_command(
        [
            binaries["fastspar_bootstrap"],
            "-c",
            input_path,
            "-n",
            str(bootstraps),
            "-p",
            bootstrap_prefix,
            "-t",
            str(args.bootstrap_threads),
            "-s",
            str(args.seed),
        ]
    )
    bootstrap_files = sorted(bootstrap_dir.glob("bootstrap_*.tsv"), key=natural_bootstrap_index)
    if len(bootstrap_files) != bootstraps:
        raise RuntimeError(f"{key}: expected {bootstraps} bootstrap tables, found {len(bootstrap_files)}")

    log(f"{key}: fitting {bootstraps} bootstrap matrices with {args.jobs} concurrent jobs")
    fit_times: list[float] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = []
        for bootstrap_path in bootstrap_files:
            index = natural_bootstrap_index(bootstrap_path)
            output_path = correlation_dir / f"correlation_{index}.tsv"
            futures.append(
                pool.submit(
                    fit_bootstrap,
                    binaries["fastspar"],
                    bootstrap_path,
                    output_path,
                    iterations=args.iterations,
                    exclusion_iterations=args.exclusion_iterations,
                    exclusion_threshold=args.exclusion_threshold,
                    threads=args.threads_per_fit,
                    seed=args.seed,
                )
            )
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            _, duration = future.result()
            fit_times.append(duration)
            completed += 1
            if completed % max(1, min(50, bootstraps // 10)) == 0 or completed == bootstraps:
                log(f"{key}: completed {completed}/{bootstraps} bootstrap fits")

    correlation_files = list(correlation_dir.glob("correlation_*.tsv"))
    if len(correlation_files) != bootstraps:
        raise RuntimeError(f"{key}: expected {bootstraps} bootstrap correlations, found {len(correlation_files)}")
    log(f"{key}: calculating exact bootstrap P values")
    run_command(
        [
            binaries["fastspar_pvalues"],
            "-c",
            input_path,
            "-r",
            correlation_path,
            "-p",
            correlation_prefix,
            "-n",
            str(bootstraps),
            "-o",
            pvalue_path,
            "-t",
            str(args.pvalue_threads),
        ]
    )

    correlation = read_square_matrix(correlation_path)
    p_values = read_square_matrix(pvalue_path).reindex(index=correlation.index, columns=correlation.columns)
    if p_values.isna().all(axis=None):
        raise RuntimeError(f"{key}: FastSpar P-value matrix is entirely missing")
    q_values = adjust_matrix(p_values)
    q_values.to_csv(qvalue_path, sep="\t", index=True, index_label="#OTU ID")
    raw_edges = edge_table(
        correlation,
        p_values,
        q_values,
        rho_threshold=args.rho_threshold,
        alpha=args.alpha,
        significance="raw_p",
    )
    q_edges = edge_table(
        correlation,
        p_values,
        q_values,
        rho_threshold=args.rho_threshold,
        alpha=args.alpha,
        significance="BH_q",
    )
    raw_edges.to_csv(raw_edges_path, sep="\t", index=False)
    q_edges.to_csv(q_edges_path, sep="\t", index=False)

    elapsed = time.perf_counter() - started
    summary: dict[str, object] = {
        "group": group,
        "stage": stage,
        "canons": row["canons"],
        "n_subjects": int(row["n_subjects"]),
        "n_samples": int(row["n_samples"]),
        "candidate_features": int(len(correlation)),
        "all_pairs": int(len(correlation) * (len(correlation) - 1) // 2),
        "raw_p_strong_edges": int(len(raw_edges)),
        "raw_p_strong_nodes": int(node_count(raw_edges)),
        "BH_q_strong_edges": int(len(q_edges)),
        "BH_q_strong_nodes": int(node_count(q_edges)),
        "rho_threshold_strict_gt": float(args.rho_threshold),
        "alpha_strict_lt": float(args.alpha),
        "bootstraps": int(bootstraps),
        "seed": int(args.seed),
        "bootstrap_fit_seed_rule": "seed + zero_based_bootstrap_index + 1",
        "iterations": int(args.iterations),
        "exclusion_iterations": int(args.exclusion_iterations),
        "exclusion_threshold": float(args.exclusion_threshold),
        "input": str(input_path.resolve()),
        "input_sha256": row["input_sha256"],
        "correlation_sha256": sha256_file(correlation_path),
        "pvalues_sha256": sha256_file(pvalue_path),
        "elapsed_seconds": float(elapsed),
        "bootstrap_fit_seconds_median": float(np.median(fit_times)) if fit_times else math.nan,
        "bootstrap_fit_seconds_max": float(np.max(fit_times)) if fit_times else math.nan,
        "method_boundary": (
            f"global top-{args.features} direct bacterial Bracken genus counts; pooled sample columns; "
            "original manuscript feature filter and seeds unavailable"
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    if not args.keep_intermediates:
        safe_remove_tree(panel_work, work_root)
    log(f"{key}: done in {elapsed / 60:.1f} min; raw={len(raw_edges)} edges, BH={len(q_edges)} edges")
    return summary


def selected_panel_rows(panel_inputs: pd.DataFrame, requested: Sequence[str] | None) -> pd.DataFrame:
    if not requested:
        return panel_inputs.copy()
    normalized = {item.replace(":", "__") for item in requested}
    available = {
        panel_key(str(row.group), str(row.stage))
        for row in panel_inputs.itertuples(index=False)
    }
    unknown = normalized.difference(available)
    if unknown:
        raise ValueError(f"unknown --panel values: {sorted(unknown)}; available: {sorted(available)}")
    keys = panel_inputs.apply(lambda row: panel_key(str(row["group"]), str(row["stage"])), axis=1)
    return panel_inputs.loc[keys.isin(normalized)].copy()


def combine_edges(result_root: Path, summaries: pd.DataFrame) -> None:
    for filename, combined_name in [
        ("edges_caption_raw_p.tsv", "edges_caption_raw_p_all_panels.tsv"),
        ("edges_methods_BH_q.tsv", "edges_methods_BH_q_all_panels.tsv"),
    ]:
        frames: list[pd.DataFrame] = []
        for row in summaries.itertuples(index=False):
            path = result_root / panel_key(str(row.group), str(row.stage)) / filename
            if not path.exists():
                continue
            frame = pd.read_csv(path, sep="\t")
            frame.insert(0, "stage", row.stage)
            frame.insert(0, "group", row.group)
            frames.append(frame)
        if frames:
            pd.concat(frames, ignore_index=True).to_csv(result_root / combined_name, sep="\t", index=False)


def load_completed_results(result_root: Path) -> CompletedFastSparResults:
    """Load and validate a complete 12-panel, non-benchmark FastSpar run.

    The renderer deliberately consumes the per-panel edge lists rather than
    recalculating thresholds from a correlation matrix.  This keeps the
    displayed caption rule and the exported Methods/BH sensitivity tied to the
    exact run recorded by :func:`run_panel`.
    """
    root = Path(result_root).resolve()
    manifest_path = root / "run_manifest.json"
    summary_path = root / "network_summary.tsv"
    missing_top_level = [
        str(path) for path in [manifest_path, summary_path] if not path.is_file()
    ]
    if missing_top_level:
        raise FileNotFoundError(
            "incomplete FastSpar result root; missing " + ", ".join(missing_top_level)
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("run_mode") != "full":
        raise ValueError(
            f"{root} is run_mode={manifest.get('run_mode')!r}; "
            "Figure 4 requires the completed full 12-panel results directory"
        )
    try:
        effective_bootstraps = int(manifest["effective_bootstraps"])
        run_seed = int(manifest["seed"])
        rho_threshold = float(manifest["rho_threshold_strict_gt"])
        alpha = float(manifest["alpha_strict_lt"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{manifest_path} lacks a valid final-run parameter set") from exc
    if effective_bootstraps <= 0:
        raise ValueError(f"{manifest_path} records a non-positive bootstrap count")
    required_manifest_text = [
        "taxonomic_source_sha256",
        "metadata_source_sha256",
        "grouping_source_sha256",
        "feature_rule",
        "sample_rule",
        "group_rule",
        "method_boundary",
        "fastspar_version",
    ]
    invalid_manifest_text = [
        field
        for field in required_manifest_text
        if not isinstance(manifest.get(field), str) or not str(manifest[field]).strip()
    ]
    if invalid_manifest_text:
        raise ValueError(
            f"{manifest_path} lacks required provenance fields: {invalid_manifest_text}"
        )

    summary = pd.read_csv(summary_path, sep="\t")
    required_summary = {
        "group",
        "stage",
        "candidate_features",
        "raw_p_strong_edges",
        "raw_p_strong_nodes",
        "BH_q_strong_edges",
        "BH_q_strong_nodes",
        "rho_threshold_strict_gt",
        "alpha_strict_lt",
        "bootstraps",
        "seed",
    }
    missing_columns = required_summary.difference(summary.columns)
    if missing_columns:
        raise ValueError(f"{summary_path} lacks columns: {sorted(missing_columns)}")
    if summary.duplicated(["group", "stage"]).any():
        raise ValueError(f"{summary_path} contains duplicate group/stage rows")

    expected_keys = {
        panel_key(group, stage)
        for group in GROUP_ORDER
        for stage in NETWORK_STAGES
    }
    observed_keys = {
        panel_key(str(row.group), str(row.stage))
        for row in summary.itertuples(index=False)
    }
    if observed_keys != expected_keys:
        raise ValueError(
            f"{summary_path} is not the complete 12-panel run; "
            f"missing={sorted(expected_keys - observed_keys)}, "
            f"unexpected={sorted(observed_keys - expected_keys)}"
        )
    if not summary["bootstraps"].astype(int).eq(effective_bootstraps).all():
        raise ValueError("panel bootstrap counts disagree with run_manifest.json")
    if not summary["seed"].astype(int).eq(run_seed).all():
        raise ValueError("panel seeds disagree with run_manifest.json")
    if not np.allclose(summary["rho_threshold_strict_gt"], rho_threshold):
        raise ValueError("panel rho thresholds disagree with run_manifest.json")
    if not np.allclose(summary["alpha_strict_lt"], alpha):
        raise ValueError("panel alpha thresholds disagree with run_manifest.json")

    edge_columns = {
        "genus1",
        "genus2",
        "rho",
        "abs_rho",
        "p_raw",
        "q_BH_within_network",
        "sign",
    }
    caption_frames: list[pd.DataFrame] = []
    methods_frames: list[pd.DataFrame] = []
    panel_summaries: dict[str, dict[str, object]] = {}
    output_names = [
        "correlation.tsv",
        "covariance.tsv",
        "pvalues_raw.tsv",
        "qvalues_BH.tsv",
        "edges_caption_raw_p.tsv",
        "edges_methods_BH_q.tsv",
        "summary.json",
    ]
    summary_lookup = summary.set_index(["group", "stage"], drop=False)
    for group in GROUP_ORDER:
        for stage in NETWORK_STAGES:
            key = panel_key(group, stage)
            panel_dir = root / key
            missing_outputs = [
                name for name in output_names if not (panel_dir / name).is_file()
            ]
            if missing_outputs:
                raise FileNotFoundError(
                    f"{key} is incomplete; missing final outputs: {missing_outputs}"
                )
            panel_summary_path = panel_dir / "summary.json"
            panel_summary = json.loads(panel_summary_path.read_text(encoding="utf-8"))
            if panel_summary.get("group") != group or panel_summary.get("stage") != stage:
                raise ValueError(f"{panel_summary_path} has inconsistent panel identity")
            if int(panel_summary.get("bootstraps", -1)) != effective_bootstraps:
                raise ValueError(f"{panel_summary_path} has an inconsistent bootstrap count")
            if int(panel_summary.get("seed", -1)) != run_seed:
                raise ValueError(f"{panel_summary_path} has an inconsistent seed")
            for filename, hash_field in [
                ("correlation.tsv", "correlation_sha256"),
                ("pvalues_raw.tsv", "pvalues_sha256"),
            ]:
                expected_hash = panel_summary.get(hash_field)
                if not isinstance(expected_hash, str) or sha256_file(panel_dir / filename) != expected_hash:
                    raise ValueError(f"{panel_dir / filename} does not match {panel_summary_path}")

            row = summary_lookup.loc[(group, stage)]
            caption = pd.read_csv(panel_dir / "edges_caption_raw_p.tsv", sep="\t")
            methods = pd.read_csv(panel_dir / "edges_methods_BH_q.tsv", sep="\t")
            for frame, label in [(caption, "caption/raw-P"), (methods, "Methods/BH-Q")]:
                missing_edge_columns = edge_columns.difference(frame.columns)
                if missing_edge_columns:
                    raise ValueError(
                        f"{key} {label} edge table lacks columns: {sorted(missing_edge_columns)}"
                    )
                if frame[["genus1", "genus2"]].isna().any(axis=None):
                    raise ValueError(f"{key} {label} edge table contains missing node names")
                for column in ["rho", "abs_rho", "p_raw", "q_BH_within_network"]:
                    frame[column] = pd.to_numeric(frame[column], errors="raise")
                unordered = frame.apply(
                    lambda item: "\0".join(sorted((str(item["genus1"]), str(item["genus2"])))),
                    axis=1,
                )
                if unordered.duplicated().any():
                    raise ValueError(f"{key} {label} edge table contains duplicate undirected edges")
                if len(frame) and not np.allclose(frame["abs_rho"], np.abs(frame["rho"])):
                    raise ValueError(f"{key} {label} edge table has inconsistent abs_rho values")
                if not (frame["abs_rho"] > rho_threshold).all():
                    raise ValueError(f"{key} {label} edge table violates the rho threshold")
            if not (caption["p_raw"] < alpha).all():
                raise ValueError(f"{key} caption edge table violates the raw-P threshold")
            if not (methods["q_BH_within_network"] < alpha).all():
                raise ValueError(f"{key} Methods edge table violates the BH-Q threshold")
            if len(caption) != int(row["raw_p_strong_edges"]):
                raise ValueError(f"{key} caption edge count disagrees with network_summary.tsv")
            if node_count(caption) != int(row["raw_p_strong_nodes"]):
                raise ValueError(f"{key} caption node count disagrees with network_summary.tsv")
            if len(methods) != int(row["BH_q_strong_edges"]):
                raise ValueError(f"{key} Methods edge count disagrees with network_summary.tsv")
            if node_count(methods) != int(row["BH_q_strong_nodes"]):
                raise ValueError(f"{key} Methods node count disagrees with network_summary.tsv")

            caption.insert(0, "stage", stage)
            caption.insert(0, "group", group)
            methods.insert(0, "stage", stage)
            methods.insert(0, "group", group)
            caption_frames.append(caption)
            methods_frames.append(methods)
            panel_summaries[key] = panel_summary

    return CompletedFastSparResults(
        root=root,
        manifest=manifest,
        summary=summary.copy(),
        caption_edges=pd.concat(caption_frames, ignore_index=True),
        methods_edges=pd.concat(methods_frames, ignore_index=True),
        panel_summaries=panel_summaries,
    )


def parser() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser(
        description="Prepare and run the 12 Figure 4 FastSpar genus networks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    cli.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parents[1])
    cli.add_argument(
        "--group-manifest",
        type=Path,
        help="CSV containing subject, definition_id and recovery_group for an alternative grouping.",
    )
    cli.add_argument(
        "--definition-id",
        help="Definition to select from --group-manifest (for example D1).",
    )
    cli.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "fastspar_fig4",
    )
    cli.add_argument("--features", type=int, default=DEFAULT_FEATURES)
    cli.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    cli.add_argument("--seed", type=int, default=DEFAULT_SEED)
    cli.add_argument("--iterations", type=int, default=50)
    cli.add_argument("--exclusion-iterations", type=int, default=10)
    cli.add_argument("--exclusion-threshold", type=float, default=0.1)
    cli.add_argument("--rho-threshold", type=float, default=DEFAULT_RHO_THRESHOLD)
    cli.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    cli.add_argument("--jobs", type=int, default=min(48, os.cpu_count() or 1))
    cli.add_argument(
        "--panel-jobs",
        type=int,
        default=1,
        help=(
            "number of network panels to run concurrently; the approximate "
            "FastSpar fit concurrency is panel-jobs x jobs"
        ),
    )
    cli.add_argument("--threads-per-fit", type=int, default=1)
    cli.add_argument("--bootstrap-threads", type=int, default=min(16, os.cpu_count() or 1))
    cli.add_argument("--pvalue-threads", type=int, default=min(32, os.cpu_count() or 1))
    cli.add_argument("--fastspar-bin-dir", type=Path)
    cli.add_argument(
        "--panel",
        action="append",
        help="run only GROUP:STAGE (repeatable); inputs for all 12 panels are still prepared",
    )
    cli.add_argument("--prepare-only", action="store_true", help="write all inputs/manifests without running FastSpar")
    cli.add_argument(
        "--benchmark",
        action="store_true",
        help="run one panel with --benchmark-bootstraps and estimate a full 12-panel run",
    )
    cli.add_argument("--benchmark-panel", default="HDR:Baseline")
    cli.add_argument("--benchmark-bootstraps", type=int, default=20)
    cli.add_argument("--keep-intermediates", action="store_true")
    cli.add_argument("--force", action="store_true", help="rerun panels even when matching completed outputs exist")
    return cli


def validate_args(args: argparse.Namespace) -> None:
    positive = {
        "features": args.features,
        "bootstraps": args.bootstraps,
        "iterations": args.iterations,
        "exclusion_iterations": args.exclusion_iterations,
        "jobs": args.jobs,
        "panel_jobs": args.panel_jobs,
        "threads_per_fit": args.threads_per_fit,
        "bootstrap_threads": args.bootstrap_threads,
        "pvalue_threads": args.pvalue_threads,
        "benchmark_bootstraps": args.benchmark_bootstraps,
    }
    bad = {name: value for name, value in positive.items() if value <= 0}
    if bad:
        raise ValueError(f"arguments must be positive: {bad}")
    if not 0 < args.alpha < 1:
        raise ValueError("--alpha must lie between 0 and 1")
    if not 0 <= args.rho_threshold < 1:
        raise ValueError("--rho-threshold must lie in [0, 1)")
    if args.benchmark and args.panel:
        raise ValueError("use --benchmark-panel, not --panel, with --benchmark")
    if bool(args.group_manifest) != bool(args.definition_id):
        raise ValueError("--group-manifest and --definition-id must be supplied together")


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    validate_args(args)
    root = args.data_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    log(f"Loading metadata and direct Bracken genus counts from {root}")
    meta = load_metadata(root)
    if args.group_manifest is None:
        groups = legacy_group_assignments(root, meta)
        grouping_source = (root / "species.csv").resolve()
        grouping_rule = "bottom/top 15 C4-C6 mean Shannon ranks from supplied species.csv"
        grouping_definition_id = "D0"
    else:
        grouping_source = args.group_manifest.resolve()
        groups = manifest_group_assignments(
            grouping_source,
            meta,
            str(args.definition_id),
        )
        grouping_definition_id = str(args.definition_id)
        first = groups.iloc[0]
        grouping_rule = (
            f"{grouping_definition_id}: {first.get('definition_full', first.get('definition_short', 'alternative recovery definition'))}; "
            "15 HDR, 15 LDR, 30 Middle from the fixed reviewer-response manifest"
        )
    counts, path_map = load_direct_genus_counts(root, meta)
    features, ranking = select_global_features(counts, meta, path_map, args.features)
    ranking.to_csv(output / f"global_top{args.features}_genera.tsv", sep="\t", index=False)
    groups.to_csv(output / "group_assignments.tsv", sep="\t", index=False)
    assignment_table = output / "group_assignments.tsv"
    panel_inputs = build_panel_inputs(output, counts, meta, groups, features)
    panel_inputs.to_csv(output / "panel_inputs.tsv", sep="\t", index=False)

    source_manifest = {
        "created_by": Path(__file__).name,
        "data_root": str(root),
        "taxonomic_source": str((root / "raw" / "combine_brk.mpa").resolve()),
        "taxonomic_source_sha256": sha256_file(root / "raw" / "combine_brk.mpa"),
        "metadata_source": str((root / "raw" / "meta.csv").resolve()),
        "metadata_source_sha256": sha256_file(root / "raw" / "meta.csv"),
        "grouping_source": str(grouping_source),
        "grouping_source_sha256": sha256_file(grouping_source),
        "grouping_definition_id": grouping_definition_id,
        "group_assignment_table": str(assignment_table.resolve()),
        "group_assignment_table_sha256": sha256_file(assignment_table),
        "feature_rule": f"top {args.features} bacterial terminal direct g__ rows by mean relative abundance across 960 CIP samples",
        "sample_rule": "Baseline=A1-A4 pooled; B6/C1/C2/C3 single time; End=C4-C6 pooled",
        "group_rule": grouping_rule,
        "method_boundary": "original feature filter, raw-vs-rarefied choice, pooling rule, and seeds were not supplied",
    }
    (output / "input_manifest.json").write_text(
        json.dumps(source_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if args.prepare_only:
        log(f"Prepared all 12 inputs in {output / 'inputs'}")
        return 0

    binaries = {
        name: resolve_binary(name, args.fastspar_bin_dir)
        for name in ["fastspar", "fastspar_bootstrap", "fastspar_pvalues"]
    }
    version = subprocess.run(
        [str(binaries["fastspar"]), "--version"],
        check=True,
        text=True,
        capture_output=True,
    )
    source_manifest["fastspar_binaries"] = {name: str(path) for name, path in binaries.items()}
    source_manifest["fastspar_version"] = (version.stdout or version.stderr).strip()

    if args.benchmark:
        selected = selected_panel_rows(panel_inputs, [args.benchmark_panel])
        bootstraps = args.benchmark_bootstraps
        result_root = output / "benchmark" / f"b{bootstraps}"
    else:
        selected = selected_panel_rows(panel_inputs, args.panel)
        bootstraps = args.bootstraps
        result_root = output / "results"
    work_root = output / "intermediates" / (
        f"benchmark_b{bootstraps}" if args.benchmark else f"full_b{bootstraps}"
    )
    result_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    source_manifest.update(
        {
            "run_mode": "benchmark" if args.benchmark else "full",
            "effective_bootstraps": bootstraps,
            "target_full_bootstraps": args.bootstraps,
            "seed": args.seed,
            "iterations": args.iterations,
            "exclusion_iterations": args.exclusion_iterations,
            "exclusion_threshold": args.exclusion_threshold,
            "rho_threshold_strict_gt": args.rho_threshold,
            "alpha_strict_lt": args.alpha,
            "jobs": args.jobs,
            "panel_jobs": args.panel_jobs,
            "threads_per_fit": args.threads_per_fit,
        }
    )
    (result_root / "run_manifest.json").write_text(
        json.dumps(source_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    summaries: list[dict[str, object]] = []
    run_started = time.perf_counter()
    rows = [row for _, row in selected.iterrows()]
    if args.panel_jobs == 1 or len(rows) == 1:
        for row in rows:
            summaries.append(run_panel(row, result_root, work_root, binaries, args, bootstraps))
    else:
        log(
            f"Running {len(rows)} panels with {args.panel_jobs} concurrent panel jobs "
            f"and up to {args.jobs} bootstrap fits per panel"
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.panel_jobs) as pool:
            future_to_key = {
                pool.submit(run_panel, row, result_root, work_root, binaries, args, bootstraps):
                panel_key(str(row["group"]), str(row["stage"]))
                for row in rows
            }
            completed_by_key: dict[str, dict[str, object]] = {}
            for future in concurrent.futures.as_completed(future_to_key):
                key = future_to_key[future]
                completed_by_key[key] = future.result()
                log(f"{key}: panel job completed")
        summaries = [
            completed_by_key[panel_key(str(row["group"]), str(row["stage"]))]
            for row in rows
        ]
    total_elapsed = time.perf_counter() - run_started
    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(result_root / "network_summary.tsv", sep="\t", index=False)
    combine_edges(result_root, summary_frame)

    if args.benchmark:
        projected = total_elapsed * (args.bootstraps / bootstraps) * 12
        benchmark = {
            "benchmark_panel": args.benchmark_panel,
            "benchmark_bootstraps": bootstraps,
            "benchmark_elapsed_seconds": total_elapsed,
            "target_bootstraps_per_panel": args.bootstraps,
            "target_panels": 12,
            "simple_linear_projection_seconds": projected,
            "simple_linear_projection_hours": projected / 3600,
            "projection_warning": "includes fixed per-panel overhead and assumes unchanged I/O/concurrency scaling",
        }
        (result_root / "benchmark_projection.json").write_text(
            json.dumps(benchmark, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        log(f"Benchmark complete; simple full-run projection: {projected / 3600:.2f} h")
    else:
        log(f"Completed {len(summary_frame)} panels in {total_elapsed / 3600:.2f} h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
