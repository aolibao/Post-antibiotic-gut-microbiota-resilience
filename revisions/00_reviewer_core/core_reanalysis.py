#!/usr/bin/env python3
"""Reviewer-focused reanalysis of longitudinal ciprofloxacin data.

This analysis is intentionally separate from ``reproduction/``.  The primary
branch follows the manuscript Methods statement by rarefying the full bacterial
species-leaf Bracken table to 2,693,000 reads (seed 20261009).  An unrarefied
branch is run as a sensitivity analysis because reverse calibration of the
published figures indicates that those figures used unrarefied leaf counts.

The sampling unit for cross-sectional inference is the participant (n=60
ciprofloxacin-exposed participants).  Repeated observations are modelled with
GEE or random-intercept mixed models.  Analyses based on the original HDR/LDR
labels are explicitly marked selection-coupled/descriptive because those labels
were reconstructed from endpoint Shannon diversity.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-reviewer-core")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from PIL import Image
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, spearmanr, t as student_t, ttest_ind
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Gaussian


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = ROOT / "reviewer_response" / "02_analysis"
TABLE_DIR = ROOT / "reviewer_response" / "03_tables"
FIGURE_DIR = ROOT / "reviewer_response" / "04_figures"
RAW_DIR = ROOT / "raw"

SEED = 20261009
RAREFACTION_DEPTH = 2_693_000
TIME_ORDER = [f"A{i}" for i in range(1, 5)] + [f"B{i}" for i in range(1, 7)] + [f"C{i}" for i in range(1, 7)]
BASELINE = ["A1", "A2", "A3", "A4"]
EXPOSURE = ["B1", "B2", "B3", "B4", "B5"]
RECOVERY = ["B6", "C1", "C2", "C3", "C4", "C5", "C6"]
ENDPOINT = ["C4", "C5", "C6"]
PRIMARY_RECOVERY_METRIC = "endpoint_mean_bray_to_baseline_centroid"
PERMANOVA_PERMUTATIONS = 9_999
CORRELATION_PERMUTATIONS = 999

BLUE = "#356A9A"
RED = "#B5524E"
TEAL = "#3E8C89"
GOLD = "#C7973A"
GREY = "#777777"
LIGHT_GREY = "#D9D9D9"


@dataclass
class BranchResult:
    name: str
    counts: pd.DataFrame
    genus_counts: pd.DataFrame
    species_to_genus: pd.Series
    sample_metrics: pd.DataFrame
    subject_metrics: pd.DataFrame
    stage_relative: dict[str, pd.DataFrame]
    stage_binary: dict[str, pd.DataFrame]
    enterotype_silhouette: float


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(pvalues), dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return out
    p = values[valid]
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[valid] = restored
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def terminal_rank(paths: pd.Index, prefix: str) -> pd.Index:
    values = (
        paths.to_series(index=paths)
        .str.rsplit("|", n=1)
        .str[-1]
        .str.removeprefix(prefix)
        .str.replace("_", " ", regex=False)
        .replace("", "Unclassified")
    )
    return pd.Index(values.to_numpy(), name=prefix[0])


def species_leaf_mask(paths: pd.Index) -> np.ndarray:
    series = paths.to_series(index=paths)
    terminal = series.str.rsplit("|", n=1).str[-1]
    return (series.str.startswith("k__Bacteria|") & terminal.str.startswith("s__")).to_numpy()


def genus_from_path(path: str) -> str:
    for part in str(path).split("|"):
        if part.startswith("g__"):
            value = part.removeprefix("g__").replace("_", " ").strip()
            return value or "Unclassified"
    return "Unclassified"


def alpha_from_counts(counts: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    values = counts.to_numpy(dtype=float)
    totals = values.sum(axis=0)
    proportions = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
    logs = np.zeros_like(proportions)
    np.log(proportions, out=logs, where=proportions > 0)
    shannon = -(proportions * logs).sum(axis=0)
    richness = (values > 0).sum(axis=0).astype(float)
    simpson = 1.0 - np.square(proportions).sum(axis=0)
    evenness = np.divide(shannon, np.log(richness), out=np.full_like(shannon, np.nan), where=richness > 1)
    idx = counts.columns
    return (
        pd.Series(shannon, index=idx),
        pd.Series(richness, index=idx),
        pd.Series(simpson, index=idx),
        pd.Series(evenness, index=idx),
    )


def bray(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.abs(a + b).sum()
    return float(np.abs(a - b).sum() / denom) if denom > 0 else 0.0


def binary_jaccard(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    union = np.logical_or(a, b).sum()
    return float(1.0 - np.logical_and(a, b).sum() / union) if union else 0.0


def normalized_auc(day: Sequence[float], value: Sequence[float]) -> float:
    x = np.asarray(day, dtype=float)
    y = np.asarray(value, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        return np.nan
    order = np.argsort(x, kind="mergesort")
    x, y = x[order], y[order]
    span = x[-1] - x[0]
    return float(np.trapz(y, x=x) / span) if span > 0 else np.nan


def zscore(series: pd.Series) -> pd.Series:
    sd = series.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (series - series.mean()) / sd


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ np.square(distance) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    positive = eigenvalues > max(1e-12, np.max(np.abs(eigenvalues)) * 1e-12)
    coords = eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])
    explained = eigenvalues[positive] / eigenvalues[positive].sum()
    negative_fraction = (
        float(np.abs(eigenvalues[eigenvalues < 0]).sum() / eigenvalues[positive].sum())
        if positive.any()
        else np.nan
    )
    return coords, explained, negative_fraction


def _permuted_labels(labels: np.ndarray, rng: np.random.Generator, strata: np.ndarray | None) -> np.ndarray:
    if strata is None:
        return rng.permutation(labels)
    out = labels.copy()
    for stratum in np.unique(strata):
        idx = np.flatnonzero(strata == stratum)
        out[idx] = rng.permutation(out[idx])
    return out


def permanova(
    distance: np.ndarray,
    labels: Sequence[str],
    permutations: int,
    seed: int,
    strata: Sequence[str] | None = None,
) -> tuple[float, float, float]:
    labels = np.asarray(labels)
    strata_arr = None if strata is None else np.asarray(strata)
    n = len(labels)
    tri = np.triu_indices(n, 1)
    ss_total = float(np.square(distance[tri]).sum() / n)

    def statistic(current: np.ndarray) -> tuple[float, float]:
        unique = np.unique(current)
        within = 0.0
        for level in unique:
            idx = np.flatnonzero(current == level)
            if len(idx) > 1:
                block = distance[np.ix_(idx, idx)]
                within += float(np.square(block[np.triu_indices(len(idx), 1)]).sum() / len(idx))
        among = ss_total - within
        df_between = len(unique) - 1
        df_within = n - len(unique)
        f_value = (among / df_between) / (within / df_within)
        return f_value, among / ss_total if ss_total > 0 else np.nan

    observed, r2 = statistic(labels)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        permuted = _permuted_labels(labels, rng, strata_arr)
        exceed += statistic(permuted)[0] >= observed - 1e-14
    return float(observed), float((exceed + 1) / (permutations + 1)), float(r2)


def permdisp(
    distance: np.ndarray,
    labels: Sequence[str],
    permutations: int,
    seed: int,
    strata: Sequence[str] | None = None,
) -> tuple[float, float, float, float]:
    labels = np.asarray(labels)
    strata_arr = None if strata is None else np.asarray(strata)
    coords, _, negative_fraction = pcoa(distance)

    def statistic(current: np.ndarray) -> tuple[float, float]:
        levels = np.unique(current)
        distances = np.empty(len(current), dtype=float)
        for level in levels:
            idx = np.flatnonzero(current == level)
            centroid = coords[idx].mean(axis=0)
            distances[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)
        grand = distances.mean()
        ss_between = sum(np.sum(current == level) * (distances[current == level].mean() - grand) ** 2 for level in levels)
        ss_within = sum(np.square(distances[current == level] - distances[current == level].mean()).sum() for level in levels)
        numerator = ss_between / (len(levels) - 1)
        denominator = ss_within / (len(current) - len(levels))
        if denominator <= np.finfo(float).eps:
            # The exact statistic is infinite when every within-group
            # distance-to-centroid is identical but group means differ.  Use
            # an epsilon-denominator limit so the permutation ordering is
            # retained while tables/tests remain finite and machine-readable.
            f_value = numerator / np.finfo(float).eps if numerator > np.finfo(float).eps else 0.0
        else:
            f_value = numerator / denominator
        means = [distances[current == level].mean() for level in levels]
        return float(f_value), float(max(means) - min(means))

    observed, mean_range = statistic(labels)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        permuted = _permuted_labels(labels, rng, strata_arr)
        exceed += statistic(permuted)[0] >= observed - 1e-14
    return observed, float((exceed + 1) / (permutations + 1)), mean_range, negative_fraction


def build_branch(
    name: str,
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    original_groups: pd.DataFrame,
    species_to_genus: pd.Series,
) -> BranchResult:
    exposed_meta = meta.loc[meta["abx"] == "Y"].copy()
    exposed_meta["day"] = pd.to_numeric(exposed_meta["day"], errors="coerce")
    exposed_meta["canon"] = pd.Categorical(exposed_meta["canon"], TIME_ORDER, ordered=True)
    sample_ids = exposed_meta["ids"].tolist()
    branch_counts = counts.loc[:, sample_ids]
    values = branch_counts.to_numpy(dtype=float)
    totals = values.sum(axis=0)
    proportions = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
    prop = pd.DataFrame(proportions, index=branch_counts.index, columns=branch_counts.columns)

    shannon, richness, simpson, evenness = alpha_from_counts(branch_counts)
    metrics = exposed_meta.copy()
    metrics["dataset_branch"] = name
    metrics["read_depth"] = metrics["ids"].map(pd.Series(totals, index=sample_ids))
    metrics["shannon"] = metrics["ids"].map(shannon)
    metrics["richness"] = metrics["ids"].map(richness)
    metrics["simpson"] = metrics["ids"].map(simpson)
    metrics["evenness"] = metrics["ids"].map(evenness)
    metrics = metrics.merge(original_groups[["recovery_group"]], left_on="subject", right_index=True, how="left")
    metrics = metrics.rename(columns={"recovery_group": "original_endpoint_group"})

    genus_counts = branch_counts.groupby(species_to_genus, sort=True).sum()
    genus_prop = genus_counts.div(genus_counts.sum(axis=0), axis=1)

    sample_index = {sample: i for i, sample in enumerate(sample_ids)}
    subject_rows: list[dict[str, object]] = []
    stage_relative: dict[str, dict[str, np.ndarray]] = {"baseline": {}, "antibiotic_end": {}, "endpoint": {}}
    stage_binary: dict[str, dict[str, np.ndarray]] = {"baseline": {}, "antibiotic_end": {}, "endpoint": {}}
    baseline_genus_by_subject: dict[str, pd.Series] = {}

    for subject, rows in metrics.groupby("subject", sort=True, observed=True):
        rows = rows.sort_values("day")
        base_rows = rows.loc[rows["canon"].isin(BASELINE)]
        recovery_rows = rows.loc[rows["canon"].isin(RECOVERY)].sort_values("day")
        endpoint_rows = rows.loc[rows["canon"].isin(ENDPOINT)]
        antibiotic_end_rows = rows.loc[rows["canon"].isin(["B5", "B6"])]
        base_ids = base_rows["ids"].tolist()
        endpoint_ids = endpoint_rows["ids"].tolist()
        recovery_ids = recovery_rows["ids"].tolist()
        base_prop = prop[base_ids].mean(axis=1).to_numpy()
        end_prop = prop[endpoint_ids].mean(axis=1).to_numpy()
        antibiotic_end_prop = prop[antibiotic_end_rows["ids"]].mean(axis=1).to_numpy()
        base_binary = (branch_counts[base_ids].to_numpy() > 0).mean(axis=1) >= 0.5
        end_binary = (branch_counts[endpoint_ids].to_numpy() > 0).mean(axis=1) >= 0.5
        antibiotic_end_binary = (branch_counts[antibiotic_end_rows["ids"]].to_numpy() > 0).mean(axis=1) >= 0.5
        stage_relative["baseline"][subject] = base_prop
        stage_relative["antibiotic_end"][subject] = antibiotic_end_prop
        stage_relative["endpoint"][subject] = end_prop
        stage_binary["baseline"][subject] = base_binary
        stage_binary["antibiotic_end"][subject] = antibiotic_end_binary
        stage_binary["endpoint"][subject] = end_binary

        baseline_shannon = float(base_rows["shannon"].mean())
        baseline_sd = float(base_rows["shannon"].std(ddof=1))
        endpoint_shannon = float(endpoint_rows["shannon"].mean())
        base_depth = float(base_rows["read_depth"].mean())
        endpoint_abs_dev = float(np.mean(np.abs(endpoint_rows["shannon"].to_numpy() - baseline_shannon)))

        sample_bray: dict[str, float] = {}
        sample_jaccard: dict[str, float] = {}
        for sample in rows["ids"]:
            vector = prop[sample].to_numpy()
            presence = branch_counts[sample].to_numpy() > 0
            sample_bray[sample] = bray(vector, base_prop)
            sample_jaccard[sample] = binary_jaccard(presence, base_binary)
        subject_mask = metrics["ids"].isin(rows["ids"])
        metrics.loc[subject_mask, "bray_to_personal_baseline_centroid"] = metrics.loc[
            subject_mask, "ids"
        ].map(sample_bray)
        metrics.loc[subject_mask, "jaccard_to_personal_baseline_consensus"] = metrics.loc[
            subject_mask, "ids"
        ].map(sample_jaccard)
        signed_deviation = metrics.loc[subject_mask, "shannon"] - baseline_shannon
        metrics.loc[subject_mask, "signed_shannon_deviation_from_baseline"] = signed_deviation
        metrics.loc[subject_mask, "absolute_shannon_deviation_from_baseline"] = signed_deviation.abs()

        recovery_bray = np.array([sample_bray[x] for x in recovery_ids])
        recovery_jaccard = np.array([sample_jaccard[x] for x in recovery_ids])
        recovery_abs_shannon = np.abs(recovery_rows["shannon"].to_numpy() - baseline_shannon)
        recovery_days = recovery_rows["day"].to_numpy(dtype=float)
        endpoint_bray_values = np.array([sample_bray[x] for x in endpoint_ids])
        endpoint_jaccard_values = np.array([sample_jaccard[x] for x in endpoint_ids])

        within_band = recovery_abs_shannon <= baseline_sd
        first_idx = np.flatnonzero(within_band)
        first_return_day = float(recovery_days[first_idx[0]]) if len(first_idx) else np.nan
        persistent_idx = [i for i in range(len(within_band)) if within_band[i:].all()]
        persistent_return_day = float(recovery_days[persistent_idx[0]]) if persistent_idx else np.nan
        last_day = float(np.nanmax(recovery_days))

        baseline_genus = genus_prop[base_ids].mean(axis=1)
        baseline_genus_by_subject[subject] = baseline_genus
        dominant_genus = str(baseline_genus.idxmax())
        bp_fraction = float(baseline_genus.get("Bacteroides", 0.0) + baseline_genus.get("Phocaeicola", 0.0))
        competitors = baseline_genus.drop(labels=[x for x in ["Bacteroides", "Phocaeicola"] if x in baseline_genus.index])
        bp_dominant = bool(bp_fraction > competitors.max())

        subject_rows.append(
            {
                "subject": subject,
                "dataset_branch": name,
                "baseline_shannon": baseline_shannon,
                "baseline_shannon_sd": baseline_sd,
                "baseline_richness": float(base_rows["richness"].mean()),
                "baseline_mean_read_depth": base_depth,
                "endpoint_shannon": endpoint_shannon,
                "endpoint_signed_shannon_change": endpoint_shannon - baseline_shannon,
                "endpoint_abs_shannon_deviation": endpoint_abs_dev,
                "endpoint_mean_bray_to_baseline_centroid": float(endpoint_bray_values.mean()),
                "endpoint_mean_jaccard_to_baseline_consensus": float(endpoint_jaccard_values.mean()),
                "baseline_endpoint_centroid_bray": bray(base_prop, end_prop),
                "baseline_endpoint_consensus_jaccard": binary_jaccard(base_binary, end_binary),
                "recovery_auc_abs_shannon_actual_day": normalized_auc(recovery_days, recovery_abs_shannon),
                "recovery_auc_bray_actual_day": normalized_auc(recovery_days, recovery_bray),
                "recovery_auc_jaccard_actual_day": normalized_auc(recovery_days, recovery_jaccard),
                "full_auc_abs_shannon_actual_day": normalized_auc(rows.loc[rows["canon"].isin(EXPOSURE + RECOVERY), "day"], np.abs(rows.loc[rows["canon"].isin(EXPOSURE + RECOVERY), "shannon"] - baseline_shannon)),
                "first_shannon_return_actual_day": first_return_day,
                "first_shannon_return_by_C6": bool(len(first_idx)),
                "persistent_shannon_return_actual_day": persistent_return_day,
                "persistent_shannon_return_by_C6": bool(len(persistent_idx)),
                "persistent_return_day_or_censor": persistent_return_day if np.isfinite(persistent_return_day) else last_day,
                "return_band_definition": "personal A1-A4 mean plus/minus 1 baseline SD",
                "dominant_genus": dominant_genus,
                "baseline_bacteroides_fraction": float(baseline_genus.get("Bacteroides", 0.0)),
                "baseline_phocaeicola_fraction": float(baseline_genus.get("Phocaeicola", 0.0)),
                "baseline_bacteroides_plus_phocaeicola_fraction": bp_fraction,
                "bacteroides_plus_phocaeicola_dominant": bp_dominant,
            }
        )

    subject = pd.DataFrame(subject_rows).set_index("subject")
    subject = subject.join(original_groups[["recovery_group", "rank", "grouping_endpoint_shannon"]])
    subject = subject.rename(columns={"recovery_group": "original_endpoint_group", "rank": "original_endpoint_rank"})

    # Exploratory three-cluster CLR enterotype, included as an adjustment/sensitivity
    # rather than treated as a biological ground truth.
    baseline_genus_frame = pd.DataFrame(baseline_genus_by_subject).T
    top_genus = baseline_genus_frame.mean(axis=0).nlargest(min(30, baseline_genus_frame.shape[1])).index
    pseudo = 0.5 / max(1.0, float(subject["baseline_mean_read_depth"].median()))
    log_values = np.log(baseline_genus_frame[top_genus].to_numpy(dtype=float) + pseudo)
    clr_values = log_values - log_values.mean(axis=1, keepdims=True)
    km = KMeans(n_clusters=3, n_init=100, random_state=SEED).fit(clr_values)
    silhouette = float(silhouette_score(clr_values, km.labels_))
    cluster_names: dict[int, str] = {}
    for cluster in sorted(set(km.labels_)):
        members = baseline_genus_frame.loc[np.asarray(km.labels_) == cluster]
        marker = str(members.mean(axis=0).idxmax())
        cluster_names[cluster] = f"E{cluster + 1}:{marker}"
    subject["exploratory_enterotype_k3"] = [cluster_names[x] for x in km.labels_]

    recovery_metrics = [
        "endpoint_abs_shannon_deviation",
        "endpoint_mean_bray_to_baseline_centroid",
        "endpoint_mean_jaccard_to_baseline_consensus",
        "recovery_auc_bray_actual_day",
    ]
    ranks = pd.DataFrame({column: subject[column].rank(method="average", ascending=True) for column in recovery_metrics})
    subject["recovery_composite_mean_rank"] = ranks.mean(axis=1)
    subject["recovery_composite_percentile_lower_is_better"] = (
        subject["recovery_composite_mean_rank"].rank(method="average") - 1
    ) / (len(subject) - 1)
    ordered = subject["recovery_composite_mean_rank"].sort_values(kind="mergesort")
    subject["new_recovery_extreme_n15"] = "Middle"
    subject.loc[ordered.head(15).index, "new_recovery_extreme_n15"] = "HighRecovery"
    subject.loc[ordered.tail(15).index, "new_recovery_extreme_n15"] = "LowRecovery"

    metrics = metrics.drop(columns=["original_endpoint_group"]).merge(
        subject[["original_endpoint_group", "new_recovery_extreme_n15"]],
        left_on="subject",
        right_index=True,
        how="left",
    )

    relative_frames = {
        stage: pd.DataFrame.from_dict(vectors, orient="index", columns=branch_counts.index)
        for stage, vectors in stage_relative.items()
    }
    binary_frames = {
        stage: pd.DataFrame.from_dict(vectors, orient="index", columns=branch_counts.index).astype(bool)
        for stage, vectors in stage_binary.items()
    }
    return BranchResult(
        name=name,
        counts=branch_counts,
        genus_counts=genus_counts,
        species_to_genus=species_to_genus,
        sample_metrics=metrics,
        subject_metrics=subject,
        stage_relative=relative_frames,
        stage_binary=binary_frames,
        enterotype_silhouette=silhouette,
    )


CONTINUOUS_OUTCOMES = [
    "endpoint_abs_shannon_deviation",
    "endpoint_mean_bray_to_baseline_centroid",
    "endpoint_mean_jaccard_to_baseline_consensus",
    "recovery_auc_abs_shannon_actual_day",
    "recovery_auc_bray_actual_day",
    "recovery_auc_jaccard_actual_day",
    "baseline_endpoint_centroid_bray",
    "baseline_endpoint_consensus_jaccard",
]


def fit_hc3_models(branch: BranchResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = branch.subject_metrics.copy()
    frame["z_baseline_shannon"] = zscore(frame["baseline_shannon"])
    frame["z_baseline_bp"] = zscore(frame["baseline_bacteroides_plus_phocaeicola_fraction"])
    frame["z_log_baseline_depth"] = zscore(np.log(frame["baseline_mean_read_depth"]))
    rows: list[dict[str, object]] = []
    formulas = {
        "unadjusted": "z_baseline_shannon",
        "adjusted_BP_continuous": "z_baseline_shannon + z_baseline_bp",
        "adjusted_BP_dominance": "z_baseline_shannon + C(bacteroides_plus_phocaeicola_dominant)",
        "adjusted_exploratory_enterotype": "z_baseline_shannon + C(exploratory_enterotype_k3)",
        "adjusted_BP_and_enterotype": "z_baseline_shannon + z_baseline_bp + C(exploratory_enterotype_k3)",
    }
    if branch.name == "unrarefied_sensitivity":
        formulas["adjusted_BP_and_read_depth"] = "z_baseline_shannon + z_baseline_bp + z_log_baseline_depth"
    for outcome in CONTINUOUS_OUTCOMES:
        frame[f"z_{outcome}"] = zscore(frame[outcome])
        for model_name, rhs in formulas.items():
            analysis = frame[[f"z_{outcome}"] + [x for x in frame.columns if x in rhs]].copy()
            # Formula evaluation handles categorical columns not captured above; use
            # the full frame after outcome-specific finite filtering.
            analysis = frame.loc[np.isfinite(frame[f"z_{outcome}"])].copy()
            fitted = smf.ols(f"z_{outcome} ~ {rhs}", data=analysis).fit(cov_type="HC3")
            ci = fitted.conf_int(alpha=0.05)
            for term in fitted.params.index:
                rows.append(
                    {
                        "dataset_branch": branch.name,
                        "outcome": outcome,
                        "model": model_name,
                        "term": term,
                        "n_subjects": int(fitted.nobs),
                        "estimate_standardized_outcome": float(fitted.params[term]),
                        "HC3_SE": float(fitted.bse[term]),
                        "statistic": float(fitted.tvalues[term]),
                        "df_residual": float(fitted.df_resid),
                        "p_value": float(fitted.pvalues[term]),
                        "CI95_low": float(ci.loc[term, 0]),
                        "CI95_high": float(ci.loc[term, 1]),
                        "adjusted_R2": float(fitted.rsquared_adj),
                        "outcome_mean_original_scale": float(frame[outcome].mean()),
                        "outcome_SD_original_scale": float(frame[outcome].std(ddof=1)),
                        "inference": "OLS with HC3 heteroscedasticity-robust covariance; participant is the unit",
                    }
                )
    table = pd.DataFrame(rows)
    table["q_BH_within_branch_model"] = table.groupby(["dataset_branch", "model"])["p_value"].transform(bh_adjust)

    # Binary return-to-band endpoint: HC3 logistic models report log odds and OR.
    binary_rows: list[dict[str, object]] = []
    for outcome in ["first_shannon_return_by_C6", "persistent_shannon_return_by_C6"]:
        for model_name, rhs in formulas.items():
            if model_name not in {"unadjusted", "adjusted_BP_continuous", "adjusted_exploratory_enterotype", "adjusted_BP_and_read_depth"}:
                continue
            try:
                fitted = smf.glm(f"{outcome} ~ {rhs}", data=frame, family=sm.families.Binomial()).fit(cov_type="HC3")
            except Exception as exc:  # complete separation is retained as an auditable failure
                binary_rows.append(
                    {
                        "dataset_branch": branch.name,
                        "outcome": outcome,
                        "model": model_name,
                        "term": "MODEL_FAILURE",
                        "n_subjects": len(frame),
                        "error": str(exc),
                    }
                )
                continue
            ci = fitted.conf_int()
            for term in fitted.params.index:
                binary_rows.append(
                    {
                        "dataset_branch": branch.name,
                        "outcome": outcome,
                        "model": model_name,
                        "term": term,
                        "n_subjects": int(fitted.nobs),
                        "estimate_log_odds": float(fitted.params[term]),
                        "odds_ratio": float(np.exp(fitted.params[term])),
                        "HC3_SE": float(fitted.bse[term]),
                        "p_value": float(fitted.pvalues[term]),
                        "CI95_log_odds_low": float(ci.loc[term, 0]),
                        "CI95_log_odds_high": float(ci.loc[term, 1]),
                        "CI95_OR_low": float(np.exp(ci.loc[term, 0])),
                        "CI95_OR_high": float(np.exp(ci.loc[term, 1])),
                        "events": int(frame[outcome].sum()),
                        "inference": "logistic GLM with HC3 covariance; return band is personal A1-A4 mean +/-1 SD",
                    }
                )
    binary = pd.DataFrame(binary_rows)
    if "p_value" in binary:
        binary["q_BH_within_branch"] = bh_adjust(binary["p_value"])
    return table, binary


def _joint_wald(fitted, term_predicate) -> tuple[float, int, float]:
    names = list(fitted.params.index)
    selected = [i for i, name in enumerate(names) if term_predicate(name)]
    if not selected:
        return np.nan, 0, np.nan
    restriction = np.zeros((len(selected), len(names)))
    for row, col in enumerate(selected):
        restriction[row, col] = 1.0
    test = fitted.wald_test(restriction, scalar=True)
    return float(np.asarray(test.statistic).squeeze()), len(selected), float(np.asarray(test.pvalue).squeeze())


def fit_repeated_models(branch: BranchResult) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = branch.sample_metrics.copy()
    data["canon"] = pd.Categorical(data["canon"], TIME_ORDER, ordered=True)
    outcome_map = {
        "shannon": "Shannon diversity",
        "signed_shannon_deviation_from_baseline": "signed Shannon deviation from personal baseline",
        "bray_to_personal_baseline_centroid": "Bray-Curtis to personal baseline centroid",
        "jaccard_to_personal_baseline_consensus": "binary Jaccard to personal baseline consensus",
    }
    model_specs = [
        ("all60_time", None, None, "confirmatory repeated-measures time effect in all exposed participants"),
        ("original_HDR_LDR", "original_endpoint_group", ["HDR", "LDR"], "selection-coupled/descriptive: group defined by endpoint Shannon"),
        ("new_recovery_extremes", "new_recovery_extreme_n15", ["HighRecovery", "LowRecovery"], "descriptive: group derived from recovery outcomes"),
    ]
    coefficient_rows: list[dict[str, object]] = []
    omnibus_rows: list[dict[str, object]] = []
    mixed_rows: list[dict[str, object]] = []
    time_term = "C(canon, Treatment(reference='A4'))"
    for outcome, outcome_label in outcome_map.items():
        for model_name, group_column, levels, interpretation in model_specs:
            subset = data.copy()
            if group_column is None:
                formula = f"{outcome} ~ {time_term}"
            else:
                subset = subset.loc[subset[group_column].isin(levels)].copy()
                subset[group_column] = pd.Categorical(subset[group_column], levels)
                formula = f"{outcome} ~ {time_term} * C({group_column})"
            gee = smf.gee(
                formula,
                groups="subject",
                data=subset,
                cov_struct=Exchangeable(),
                family=Gaussian(),
            ).fit(maxiter=200)
            ci = gee.conf_int()
            for term in gee.params.index:
                coefficient_rows.append(
                    {
                        "dataset_branch": branch.name,
                        "outcome": outcome,
                        "outcome_label": outcome_label,
                        "model": model_name,
                        "term": term,
                        "n_subjects": subset["subject"].nunique(),
                        "n_observations": len(subset),
                        "estimate": float(gee.params[term]),
                        "robust_SE": float(gee.bse[term]),
                        "z": float(gee.tvalues[term]),
                        "p_value": float(gee.pvalues[term]),
                        "CI95_low": float(ci.loc[term, 0]),
                        "CI95_high": float(ci.loc[term, 1]),
                        "working_correlation": "exchangeable",
                        "interpretation": interpretation,
                    }
                )
            if group_column is None:
                stat, df, pvalue = _joint_wald(gee, lambda term: term.startswith("C(canon"))
                tested = "all categorical time coefficients jointly zero"
            else:
                stat, df, pvalue = _joint_wald(gee, lambda term: ":" in term)
                tested = "all group-by-time interaction coefficients jointly zero"
            omnibus_rows.append(
                {
                    "dataset_branch": branch.name,
                    "outcome": outcome,
                    "model": model_name,
                    "n_subjects": subset["subject"].nunique(),
                    "n_observations": len(subset),
                    "wald_chi2": stat,
                    "df": df,
                    "p_value": pvalue,
                    "tested": tested,
                    "interpretation": interpretation,
                }
            )

            # Random-intercept MixedLM sensitivity.  Retain convergence/failure
            # explicitly; GEE above is the primary repeated-measures inference.
            if outcome == "shannon" or (model_name == "all60_time" and outcome == "bray_to_personal_baseline_centroid"):
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        mixed_model = smf.mixedlm(formula, subset, groups=subset["subject"])
                        optimizer = "lbfgs"
                        try:
                            mixed = mixed_model.fit(method=optimizer, reml=False, maxiter=500, disp=False)
                        except (np.linalg.LinAlgError, ValueError):
                            # A near-zero random-intercept variance can make the
                            # LBFGS Hessian singular. Powell is a prespecified
                            # fallback; boundary warnings remain in the audit.
                            optimizer = "powell_fallback_after_lbfgs"
                            mixed = mixed_model.fit(method="powell", reml=False, maxiter=1000, disp=False)
                    ci_m = mixed.conf_int()
                    for term in mixed.fe_params.index:
                        mixed_rows.append(
                            {
                                "dataset_branch": branch.name,
                                "outcome": outcome,
                                "model": model_name,
                                "term": term,
                                "n_subjects": subset["subject"].nunique(),
                                "n_observations": len(subset),
                                "estimate": float(mixed.fe_params[term]),
                                "SE": float(mixed.bse_fe[term]),
                                "p_value": float(mixed.pvalues[term]),
                                "CI95_low": float(ci_m.loc[term, 0]),
                                "CI95_high": float(ci_m.loc[term, 1]),
                                "converged": bool(mixed.converged),
                                "warning_count": len(caught),
                                "warnings": " | ".join(str(x.message) for x in caught),
                                "optimizer": optimizer,
                                "interpretation": "random-intercept MixedLM sensitivity",
                            }
                        )
                except Exception as exc:
                    mixed_rows.append(
                        {
                            "dataset_branch": branch.name,
                            "outcome": outcome,
                            "model": model_name,
                            "term": "MODEL_FAILURE",
                            "n_subjects": subset["subject"].nunique(),
                            "n_observations": len(subset),
                            "converged": False,
                            "warnings": str(exc),
                            "optimizer": "lbfgs_then_powell_if_supported",
                            "interpretation": "random-intercept MixedLM sensitivity",
                        }
                    )
    coefficients = pd.DataFrame(coefficient_rows)
    coefficients["q_BH_within_branch_model"] = coefficients.groupby(["dataset_branch", "model"])["p_value"].transform(bh_adjust)
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["q_BH_within_branch"] = omnibus.groupby("dataset_branch")["p_value"].transform(bh_adjust)
    return coefficients, omnibus, pd.DataFrame(mixed_rows)


def hedges_g(high: np.ndarray, low: np.ndarray) -> float:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    n1, n2 = len(high), len(low)
    pooled = math.sqrt(((n1 - 1) * np.var(high, ddof=1) + (n2 - 1) * np.var(low, ddof=1)) / (n1 + n2 - 2))
    if pooled == 0:
        return np.nan
    d = (np.mean(high) - np.mean(low)) / pooled
    correction = 1 - 3 / (4 * (n1 + n2) - 9)
    return float(correction * d)


def cutoff_and_overlap(branch: BranchResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    subject = branch.subject_metrics.copy()
    schemes = [("extreme_n10", 10), ("extreme_n15", 15), ("extreme_n20", 20), ("quartile_extremes", 15), ("tertile_extremes", 20)]
    outcomes = ["baseline_shannon", "endpoint_shannon"] + CONTINUOUS_OUTCOMES
    rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    selection_definitions = [
        ("original_endpoint_shannon", "endpoint_shannon", False, "selection-coupled/descriptive"),
        ("new_primary_recovery_bray", PRIMARY_RECOVERY_METRIC, True, "outcome-derived recovery groups; descriptive for selected metric"),
        ("new_composite_recovery", "recovery_composite_mean_rank", True, "outcome-derived recovery groups; descriptive"),
    ]
    original_low15 = set(subject.loc[subject["original_endpoint_group"] == "LDR"].index)
    original_high15 = set(subject.loc[subject["original_endpoint_group"] == "HDR"].index)
    for scheme, n in schemes:
        for selection_name, selector, lower_is_high_recovery, caveat in selection_definitions:
            ordered = subject[selector].sort_values(kind="mergesort")
            if lower_is_high_recovery:
                high_idx, low_idx = ordered.head(n).index, ordered.tail(n).index
                contrast = "high recovery minus low recovery"
            else:
                high_idx, low_idx = ordered.tail(n).index, ordered.head(n).index
                contrast = "high endpoint Shannon minus low endpoint Shannon"
            for outcome in outcomes:
                high = subject.loc[high_idx, outcome].dropna().to_numpy(dtype=float)
                low = subject.loc[low_idx, outcome].dropna().to_numpy(dtype=float)
                test = ttest_ind(high, low, equal_var=False)
                difference = float(np.mean(high) - np.mean(low))
                se = math.sqrt(np.var(high, ddof=1) / len(high) + np.var(low, ddof=1) / len(low))
                df_num = (np.var(high, ddof=1) / len(high) + np.var(low, ddof=1) / len(low)) ** 2
                df_den = (np.var(high, ddof=1) / len(high)) ** 2 / (len(high) - 1) + (np.var(low, ddof=1) / len(low)) ** 2 / (len(low) - 1)
                df = df_num / df_den if df_den > 0 else np.nan
                critical = student_t.ppf(0.975, df) if np.isfinite(df) else np.nan
                rows.append(
                    {
                        "dataset_branch": branch.name,
                        "scheme": scheme,
                        "n_per_extreme": n,
                        "selection_definition": selection_name,
                        "outcome": outcome,
                        "contrast": contrast,
                        "mean_high_label": float(np.mean(high)),
                        "mean_low_label": float(np.mean(low)),
                        "difference": difference,
                        "CI95_low": difference - critical * se,
                        "CI95_high": difference + critical * se,
                        "welch_t": float(test.statistic),
                        "welch_df": df,
                        "p_value": float(test.pvalue),
                        "hedges_g": hedges_g(high, low),
                        "interpretation": caveat,
                    }
                )
            overlap_rows.append(
                {
                    "dataset_branch": branch.name,
                    "scheme": scheme,
                    "n_per_extreme": n,
                    "new_group_definition": selection_name,
                    "high_recovery_or_endpoint_overlap_original_HDR": len(set(high_idx) & original_high15),
                    "low_recovery_or_endpoint_overlap_original_LDR": len(set(low_idx) & original_low15),
                    "original_HDR_size": len(original_high15),
                    "original_LDR_size": len(original_low15),
                    "interpretation": caveat,
                }
            )
    table = pd.DataFrame(rows)
    table["q_BH_within_branch_selection_outcome"] = table.groupby(
        ["dataset_branch", "selection_definition", "outcome"]
    )["p_value"].transform(bh_adjust)

    recovery_metrics = [
        "endpoint_abs_shannon_deviation",
        "endpoint_mean_bray_to_baseline_centroid",
        "endpoint_mean_jaccard_to_baseline_consensus",
        "recovery_auc_abs_shannon_actual_day",
        "recovery_auc_bray_actual_day",
        "recovery_auc_jaccard_actual_day",
        "recovery_composite_mean_rank",
    ]
    for metric in recovery_metrics:
        ordered = subject[metric].sort_values(kind="mergesort")
        for n in [10, 15, 20]:
            high_recovery, low_recovery = set(ordered.head(n).index), set(ordered.tail(n).index)
            overlap_rows.append(
                {
                    "dataset_branch": branch.name,
                    "scheme": f"metric_extreme_n{n}",
                    "n_per_extreme": n,
                    "new_group_definition": metric,
                    "high_recovery_or_endpoint_overlap_original_HDR": len(high_recovery & original_high15),
                    "low_recovery_or_endpoint_overlap_original_LDR": len(low_recovery & original_low15),
                    "original_HDR_size": len(original_high15),
                    "original_LDR_size": len(original_low15),
                    "interpretation": "overlap only; original labels are endpoint-Shannon selection-coupled",
                }
            )
    return table, pd.DataFrame(overlap_rows)


def pairwise_distance(frame: pd.DataFrame, metric: str) -> np.ndarray:
    if metric == "braycurtis":
        return squareform(pdist(frame.to_numpy(dtype=float), metric="braycurtis"))
    if metric == "jaccard":
        return squareform(pdist(frame.to_numpy(dtype=bool), metric="jaccard"))
    raise ValueError(metric)


def distance_tests(branch: BranchResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    subject = branch.subject_metrics
    rows: list[dict[str, object]] = []
    coordinate_rows: list[pd.DataFrame] = []
    original_ids = subject.index[subject["original_endpoint_group"].isin(["HDR", "LDR"])]
    new_ids = subject.index[subject["new_recovery_extreme_n15"].isin(["HighRecovery", "LowRecovery"])]

    for distance_name, matrices in [("Bray-Curtis", branch.stage_relative), ("binary Jaccard", branch.stage_binary)]:
        metric = "braycurtis" if distance_name == "Bray-Curtis" else "jaccard"
        for stage, frame in matrices.items():
            for comparison, ids, group_col, caveat in [
                ("original_HDR_vs_LDR", original_ids, "original_endpoint_group", "selection-coupled/descriptive"),
                ("new_recovery_extremes", new_ids, "new_recovery_extreme_n15", "outcome-derived groups; descriptive"),
            ]:
                selected = frame.loc[ids]
                labels = subject.loc[ids, group_col].to_numpy(dtype=str)
                distance = pairwise_distance(selected, metric)
                f_value, p_value, r2 = permanova(distance, labels, PERMANOVA_PERMUTATIONS, SEED + len(rows))
                disp_f, disp_p, disp_range, neg = permdisp(distance, labels, PERMANOVA_PERMUTATIONS, SEED + 1000 + len(rows))
                rows.append(
                    {
                        "dataset_branch": branch.name,
                        "distance": distance_name,
                        "stage": stage,
                        "comparison": comparison,
                        "n_subjects": len(ids),
                        "PERMANOVA_pseudoF": f_value,
                        "PERMANOVA_R2": r2,
                        "PERMANOVA_p": p_value,
                        "PERMDISP_F": disp_f,
                        "PERMDISP_p": disp_p,
                        "PERMDISP_group_mean_distance_range": disp_range,
                        "PCoA_negative_eigenvalue_fraction": neg,
                        "permutations": PERMANOVA_PERMUTATIONS,
                        "permutation_unit": "participant label",
                        "interpretation": caveat,
                    }
                )

        # Paired baseline-to-endpoint location and dispersion test in all 60 subjects.
        baseline = matrices["baseline"]
        endpoint = matrices["endpoint"]
        combined = pd.concat([baseline, endpoint], axis=0)
        labels = np.array(["baseline"] * len(baseline) + ["endpoint"] * len(endpoint))
        strata = np.array(list(baseline.index) + list(endpoint.index))
        distance = pairwise_distance(combined, metric)
        f_value, p_value, r2 = permanova(distance, labels, PERMANOVA_PERMUTATIONS, SEED + 2000 + len(rows), strata)
        disp_f, disp_p, disp_range, neg = permdisp(distance, labels, PERMANOVA_PERMUTATIONS, SEED + 3000 + len(rows), strata)
        rows.append(
            {
                "dataset_branch": branch.name,
                "distance": distance_name,
                "stage": "paired_baseline_vs_endpoint",
                "comparison": "all60_within_subject_phase",
                "n_subjects": len(baseline),
                "PERMANOVA_pseudoF": f_value,
                "PERMANOVA_R2": r2,
                "PERMANOVA_p": p_value,
                "PERMDISP_F": disp_f,
                "PERMDISP_p": disp_p,
                "PERMDISP_group_mean_distance_range": disp_range,
                "PCoA_negative_eigenvalue_fraction": neg,
                "permutations": PERMANOVA_PERMUTATIONS,
                "permutation_unit": "phase labels swapped within participant",
                "interpretation": "paired, participant-restricted permutation",
            }
        )
        coords, explained, _ = pcoa(distance)
        coord = pd.DataFrame(
            {
                "dataset_branch": branch.name,
                "distance": distance_name,
                "subject": strata,
                "phase": labels,
                "PCo1": coords[:, 0],
                "PCo2": coords[:, 1],
                "PCo1_explained_positive_fraction": explained[0],
                "PCo2_explained_positive_fraction": explained[1],
            }
        )
        coordinate_rows.append(coord)
    table = pd.DataFrame(rows)
    table["PERMANOVA_q_BH_within_branch"] = table.groupby("dataset_branch")["PERMANOVA_p"].transform(bh_adjust)
    table["PERMDISP_q_BH_within_branch"] = table.groupby("dataset_branch")["PERMDISP_p"].transform(bh_adjust)
    return table, pd.concat(coordinate_rows, ignore_index=True)


def _rank_correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan
    xr = rankdata(x)
    yr = rankdata(y)
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    denom = np.linalg.norm(xr) * np.linalg.norm(yr)
    return float(np.dot(xr, yr) / denom) if denom > 0 else np.nan


def compositional_correlations(branch: BranchResult) -> pd.DataFrame:
    metrics = branch.sample_metrics.set_index("ids")
    subjects = sorted(metrics["subject"].unique())
    genus = branch.genus_counts
    rel = genus.div(genus.sum(axis=0), axis=1)
    prevalence = (genus > 0).mean(axis=1)
    mean_rel = rel.mean(axis=1)
    focal_genera = prevalence.index[(prevalence >= 0.10) & (mean_rel >= 0.0001) & (prevalence.index != "Unclassified")]
    all_genus = genus.index

    count_values = branch.counts.to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        clogc_values = np.where(count_values > 0, count_values * np.log(count_values), 0.0)
    total_clogc = pd.Series(clogc_values.sum(axis=0), index=branch.counts.columns)
    genus_clogc = pd.DataFrame(clogc_values, index=branch.counts.index, columns=branch.counts.columns).groupby(
        branch.species_to_genus, sort=True
    ).sum()
    full_shannon, _, _, _ = alpha_from_counts(branch.counts)

    rows: list[dict[str, object]] = []
    for time_index, canon in enumerate(TIME_ORDER):
        ordered = metrics.loc[metrics["canon"].astype(str) == canon].sort_values("subject")
        if list(ordered["subject"]) != subjects:
            raise ValueError(f"{canon}: expected one sample for each exposed participant")
        ids = ordered.index.tolist()
        genus_matrix = genus.loc[all_genus, ids].T.to_numpy(dtype=float)
        log_genus = np.log(genus_matrix + 0.5)
        clr = log_genus - log_genus.mean(axis=1, keepdims=True)
        genus_position = {name: i for i, name in enumerate(all_genus)}
        totals = branch.counts[ids].sum(axis=0).to_numpy(dtype=float)
        total_c = total_clogc[ids].to_numpy(dtype=float)
        full_h = full_shannon[ids].to_numpy(dtype=float)
        permutation_rng = np.random.default_rng(SEED + 5000 + time_index)
        permutations = np.vstack([permutation_rng.permutation(len(ids)) for _ in range(CORRELATION_PERMUTATIONS)])

        for focal in focal_genera:
            focal_counts = genus.loc[focal, ids].to_numpy(dtype=float)
            remaining = totals - focal_counts
            focal_clogc = genus_clogc.loc[focal, ids].to_numpy(dtype=float)
            excluded_h = np.full(len(ids), np.nan)
            valid = remaining > 0
            excluded_h[valid] = np.log(remaining[valid]) - (total_c[valid] - focal_clogc[valid]) / remaining[valid]
            x = clr[:, genus_position[focal]]
            rho_corrected = _rank_correlation(x, excluded_h)
            rho_conventional = _rank_correlation(rel.loc[focal, ids].to_numpy(dtype=float), full_h)
            xr = rankdata(x) - np.mean(rankdata(x))
            yr = rankdata(excluded_h) - np.mean(rankdata(excluded_h))
            denom = np.linalg.norm(xr) * np.linalg.norm(yr)
            permuted_rho = (yr[permutations] @ xr) / denom if denom > 0 else np.full(CORRELATION_PERMUTATIONS, np.nan)
            permutation_p = float((1 + np.sum(np.abs(permuted_rho) >= abs(rho_corrected) - 1e-14)) / (CORRELATION_PERMUTATIONS + 1))
            analytic = spearmanr(x, excluded_h)
            rows.append(
                {
                    "dataset_branch": branch.name,
                    "canon": canon,
                    "genus": focal,
                    "n_subjects": len(ids),
                    "prevalence_all_exposed_samples": float(prevalence[focal]),
                    "mean_relative_abundance_all_exposed_samples": float(mean_rel[focal]),
                    "rho_relative_abundance_vs_unmodified_shannon": rho_conventional,
                    "rho_CLR_vs_focal_genus_excluded_shannon": rho_corrected,
                    "p_spearman_asymptotic": float(analytic.pvalue),
                    "p_subject_label_permutation": permutation_p,
                    "permutations": CORRELATION_PERMUTATIONS,
                    "CLR_pseudocount_reads": 0.5,
                    "focal_exclusion": "all species assigned to the focal genus removed before Shannon recomputation",
                    "sampling_unit": "participant; one sample per canonical time",
                }
            )
    table = pd.DataFrame(rows)
    table["q_BH_within_time"] = table.groupby(["dataset_branch", "canon"])["p_subject_label_permutation"].transform(bh_adjust)
    table["q_BH_global_16_times"] = table.groupby("dataset_branch")["p_subject_label_permutation"].transform(bh_adjust)
    return table


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.0,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def export_figure(fig: mpl.figure.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    temp = stem.with_name(stem.name + "_temp.png")
    fig.savefig(temp, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    with Image.open(temp) as opened:
        rgb = opened.convert("RGB")
        rgb.save(stem.with_suffix(".tif"), dpi=(600, 600), compression="tiff_lzw")
    temp.unlink()


def add_panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.14, 1.07, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def render_core_figure(
    primary: BranchResult,
    hc3: pd.DataFrame,
    cutoff: pd.DataFrame,
    overlap: pd.DataFrame,
    gee_omnibus: pd.DataFrame,
) -> None:
    subject = primary.subject_metrics
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.3), layout="constrained")

    ax = axes[0, 0]
    colors = np.where(subject["bacteroides_plus_phocaeicola_dominant"], RED, BLUE)
    ax.scatter(subject["baseline_shannon"], subject[PRIMARY_RECOVERY_METRIC], c=colors, s=20, alpha=0.82, edgecolor="white", linewidth=0.25)
    fit = hc3.loc[
        (hc3["dataset_branch"] == primary.name)
        & (hc3["outcome"] == PRIMARY_RECOVERY_METRIC)
        & (hc3["model"] == "adjusted_BP_continuous")
        & (hc3["term"] == "z_baseline_shannon")
    ].iloc[0]
    raw_fit = np.polyfit(subject["baseline_shannon"], subject[PRIMARY_RECOVERY_METRIC], 1)
    xx = np.linspace(subject["baseline_shannon"].min(), subject["baseline_shannon"].max(), 100)
    ax.plot(xx, np.polyval(raw_fit, xx), color="black", linewidth=0.9)
    ax.set(xlabel="Baseline Shannon (A1–A4 mean)", ylabel="Endpoint Bray–Curtis\nto personal baseline")
    ax.set_title(f"All exposed participants (n=60)\nAdjusted β={fit['estimate_standardized_outcome']:.2f}, p={fit['p_value']:.3g}")
    ax.scatter([], [], color=RED, label="Bacteroides + Phocaeicola dominant")
    ax.scatter([], [], color=BLUE, label="Other dominant genus")
    ax.legend(fontsize=5.6, loc="best")
    add_panel_label(ax, "A")

    ax = axes[0, 1]
    rho, pvalue = spearmanr(subject["original_endpoint_rank"], subject["recovery_composite_mean_rank"])
    ax.scatter(subject["original_endpoint_rank"], subject["recovery_composite_mean_rank"], c=colors, s=18, alpha=0.8, edgecolor="none")
    ax.set(xlabel="Original endpoint-Shannon rank", ylabel="New recovery composite rank\n(lower = closer return)")
    ax.set_title(f"Rank agreement: ρ={rho:.2f}, p={pvalue:.3g}")
    add_panel_label(ax, "B")

    ax = axes[0, 2]
    display_outcomes = [
        ("endpoint_abs_shannon_deviation", "|Δ Shannon|"),
        (PRIMARY_RECOVERY_METRIC, "Bray endpoint"),
        ("endpoint_mean_jaccard_to_baseline_consensus", "Jaccard endpoint"),
        ("recovery_auc_bray_actual_day", "Bray AUC"),
    ]
    for offset, (branch, color) in enumerate([("rarefied_primary", TEAL), ("unrarefied_sensitivity", GOLD)]):
        part = hc3.loc[(hc3["dataset_branch"] == branch) & (hc3["model"] == "adjusted_BP_continuous") & (hc3["term"] == "z_baseline_shannon")].set_index("outcome")
        y = np.arange(len(display_outcomes)) + (offset - 0.5) * 0.16
        estimates = np.array([part.loc[key, "estimate_standardized_outcome"] for key, _ in display_outcomes])
        low = np.array([part.loc[key, "CI95_low"] for key, _ in display_outcomes])
        high = np.array([part.loc[key, "CI95_high"] for key, _ in display_outcomes])
        ax.errorbar(estimates, y, xerr=[estimates - low, high - estimates], fmt="o", color=color, ms=4, lw=0.9, capsize=2, label=branch.replace("_", " "))
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(np.arange(len(display_outcomes)), [label for _, label in display_outcomes])
    ax.invert_yaxis()
    ax.set_xlabel("Standardized β for baseline Shannon\n(adjusted for Bacteroides + Phocaeicola)")
    ax.set_title("Continuous-model effects (95% CI)")
    ax.legend(fontsize=5.5)
    add_panel_label(ax, "C")

    ax = axes[1, 0]
    part = cutoff.loc[
        (cutoff["dataset_branch"] == primary.name)
        & (cutoff["selection_definition"] == "original_endpoint_shannon")
        & (cutoff["outcome"].isin([PRIMARY_RECOVERY_METRIC, "recovery_auc_bray_actual_day"]))
        & (cutoff["scheme"].isin(["extreme_n10", "extreme_n15", "extreme_n20", "quartile_extremes", "tertile_extremes"]))
    ].copy()
    x_order = ["extreme_n10", "extreme_n15", "extreme_n20", "quartile_extremes", "tertile_extremes"]
    labels = ["n=10", "n=15", "n=20", "quartile", "tertile"]
    for outcome, color, marker, label in [
        (PRIMARY_RECOVERY_METRIC, BLUE, "o", "Endpoint Bray"),
        ("recovery_auc_bray_actual_day", RED, "s", "Recovery Bray AUC"),
    ]:
        values = part.loc[part["outcome"] == outcome].set_index("scheme").reindex(x_order)
        ax.plot(range(len(x_order)), values["difference"], color=color, marker=marker, ms=4, label=label)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(len(x_order)), labels, rotation=35, ha="right")
    ax.set_ylabel("Original HDR − LDR")
    ax.set_title("Extreme-cutoff sensitivity\n(selection-coupled)")
    ax.legend(fontsize=5.5)
    add_panel_label(ax, "D")

    ax = axes[1, 1]
    selected = primary.sample_metrics.loc[primary.sample_metrics["original_endpoint_group"].isin(["HDR", "LDR"])]
    for group, color in [("HDR", RED), ("LDR", BLUE)]:
        part = selected.loc[selected["original_endpoint_group"] == group]
        summary = part.groupby("canon", observed=True)["signed_shannon_deviation_from_baseline"].agg(["mean", "sem"]).reindex(TIME_ORDER)
        x = np.arange(len(TIME_ORDER))
        ax.plot(x, summary["mean"], color=color, marker="o", ms=2.5, lw=1, label=group)
        ax.fill_between(x, summary["mean"] - 1.96 * summary["sem"], summary["mean"] + 1.96 * summary["sem"], color=color, alpha=0.15, lw=0)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(np.arange(len(TIME_ORDER)), TIME_ORDER, rotation=45, ha="right")
    ax.set_ylabel("Shannon − personal baseline")
    omnibus = gee_omnibus.loc[
        (gee_omnibus["dataset_branch"] == primary.name)
        & (gee_omnibus["outcome"] == "signed_shannon_deviation_from_baseline")
        & (gee_omnibus["model"] == "original_HDR_LDR")
    ].iloc[0]
    ax.set_title(f"GEE group × time: p={omnibus['p_value']:.3g}\nselection-coupled; mean ± 95% CI")
    ax.legend(fontsize=5.5)
    add_panel_label(ax, "E")

    ax = axes[1, 2]
    metrics = [
        "endpoint_abs_shannon_deviation",
        PRIMARY_RECOVERY_METRIC,
        "endpoint_mean_jaccard_to_baseline_consensus",
        "recovery_auc_abs_shannon_actual_day",
        "recovery_auc_bray_actual_day",
        "recovery_composite_mean_rank",
    ]
    labels_y = ["|Δ Shannon|", "Bray endpoint", "Jaccard endpoint", "Shannon AUC", "Bray AUC", "Composite"]
    ov = overlap.loc[
        (overlap["dataset_branch"] == primary.name)
        & (overlap["scheme"] == "metric_extreme_n15")
        & (overlap["new_group_definition"].isin(metrics))
    ].set_index("new_group_definition").reindex(metrics)
    image_data = ov[["high_recovery_or_endpoint_overlap_original_HDR", "low_recovery_or_endpoint_overlap_original_LDR"]].to_numpy(dtype=float)
    im = ax.imshow(image_data, cmap="Blues", vmin=0, vmax=15, aspect="auto")
    for i in range(image_data.shape[0]):
        for j in range(image_data.shape[1]):
            ax.text(j, i, f"{int(image_data[i, j])}/15", ha="center", va="center", fontsize=6, color="white" if image_data[i, j] > 8 else "black")
    ax.set_xticks([0, 1], ["High recovery\n∩ original HDR", "Low recovery\n∩ original LDR"], rotation=20, ha="right")
    ax.set_yticks(range(len(labels_y)), labels_y)
    ax.set_title("Original vs recovery-defined groups")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03, label="Overlap count")
    add_panel_label(ax, "F")

    fig.suptitle("Recovery reanalysis in all 60 ciprofloxacin-exposed participants", fontsize=9.2, fontweight="bold")
    export_figure(fig, FIGURE_DIR / "Reviewer_core_recovery_reanalysis")


def render_distance_figure(primary: BranchResult, tests: pd.DataFrame, coordinates: pd.DataFrame) -> None:
    coord = coordinates.loc[(coordinates["dataset_branch"] == primary.name) & (coordinates["distance"] == "Bray-Curtis")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.15, 1]}, layout="constrained")
    ax = axes[0]
    for subject, part in coord.groupby("subject"):
        part = part.set_index("phase").loc[["baseline", "endpoint"]]
        ax.plot(part["PCo1"], part["PCo2"], color=LIGHT_GREY, lw=0.45, zorder=1)
    for phase, color, marker in [("baseline", BLUE, "o"), ("endpoint", RED, "^")]:
        part = coord.loc[coord["phase"] == phase]
        ax.scatter(part["PCo1"], part["PCo2"], color=color, marker=marker, s=16, alpha=0.75, edgecolor="white", linewidth=0.2, label=phase, zorder=2)
    explained1 = coord["PCo1_explained_positive_fraction"].iloc[0] * 100
    explained2 = coord["PCo2_explained_positive_fraction"].iloc[0] * 100
    pair = tests.loc[(tests["dataset_branch"] == primary.name) & (tests["distance"] == "Bray-Curtis") & (tests["comparison"] == "all60_within_subject_phase")].iloc[0]
    ax.set(xlabel=f"PCo1 ({explained1:.1f}% positive inertia)", ylabel=f"PCo2 ({explained2:.1f}% positive inertia)")
    ax.set_title(f"Personal centroids, n=60\npaired PERMANOVA p={pair['PERMANOVA_p']:.4g}; PERMDISP p={pair['PERMDISP_p']:.4g}")
    ax.legend(fontsize=6)
    add_panel_label(ax, "A")

    ax = axes[1]
    part = tests.loc[(tests["dataset_branch"] == primary.name) & (tests["comparison"] == "original_HDR_vs_LDR")].copy()
    part["label"] = part["distance"].str.replace("binary ", "") + "\n" + part["stage"].str.replace("_", " ")
    y = np.arange(len(part))
    colors = np.where(part["PERMDISP_p"] < 0.05, RED, BLUE)
    ax.barh(y, part["PERMANOVA_R2"], color=colors, alpha=0.82)
    for yi, (_, row) in zip(y, part.iterrows()):
        ax.text(row["PERMANOVA_R2"] + 0.004, yi, f"P={row['PERMANOVA_p']:.3g}\nD={row['PERMDISP_p']:.3g}", va="center", fontsize=5.4)
    ax.set_yticks(y, part["label"])
    ax.invert_yaxis()
    ax.set_xlabel("PERMANOVA R²")
    ax.set_title("Original HDR/LDR (selection-coupled)\nlabels: location P; dispersion D")
    ax.set_xlim(0, max(0.16, float(part["PERMANOVA_R2"].max()) * 1.45))
    add_panel_label(ax, "B")
    export_figure(fig, FIGURE_DIR / "Reviewer_distance_PERMANOVA_PERMDISP")


def render_compositional_figure(correlations: pd.DataFrame) -> None:
    primary = correlations.loc[correlations["dataset_branch"] == "rarefied_primary"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), layout="constrained")
    ax = axes[0]
    ax.scatter(
        primary["rho_relative_abundance_vs_unmodified_shannon"],
        primary["rho_CLR_vs_focal_genus_excluded_shannon"],
        s=8,
        color=GREY,
        alpha=0.35,
        edgecolor="none",
    )
    focus = primary.loc[primary["genus"].isin(["Bacteroides", "Phocaeicola", "Faecalibacterium"])]
    focus_colors = {"Bacteroides": RED, "Phocaeicola": GOLD, "Faecalibacterium": TEAL}
    for genus, part in focus.groupby("genus"):
        ax.scatter(part["rho_relative_abundance_vs_unmodified_shannon"], part["rho_CLR_vs_focal_genus_excluded_shannon"], s=16, color=focus_colors[genus], label=genus, edgecolor="white", linewidth=0.25)
    lim = 1.02
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.6, ls="--")
    ax.axhline(0, color=LIGHT_GREY, lw=0.6)
    ax.axvline(0, color=LIGHT_GREY, lw=0.6)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="Conventional ρ:\nrelative abundance vs Shannon", ylabel="Corrected ρ:\nCLR vs focal-excluded Shannon")
    ax.set_title("Compositionality-aware reanalysis")
    ax.legend(fontsize=5.5)
    add_panel_label(ax, "A")

    ax = axes[1]
    summary = primary.groupby("canon", observed=True).agg(
        tested=("genus", "size"),
        nominal=("p_subject_label_permutation", lambda x: int((x < 0.05).sum())),
        within_time_FDR=("q_BH_within_time", lambda x: int((x < 0.05).sum())),
        global_FDR=("q_BH_global_16_times", lambda x: int((x < 0.05).sum())),
    ).reindex(TIME_ORDER)
    x = np.arange(len(TIME_ORDER))
    ax.plot(x, summary["nominal"], color=GREY, marker="o", ms=2.5, label="Permutation P<0.05")
    ax.plot(x, summary["within_time_FDR"], color=GOLD, marker="s", ms=2.5, label="Within-time BH q<0.05")
    ax.plot(x, summary["global_FDR"], color=TEAL, marker="^", ms=2.5, label="Global BH q<0.05")
    ax.set_xticks(x, TIME_ORDER, rotation=45, ha="right")
    ax.set_ylabel("Number of genus associations")
    ax.set_title(f"Subject-label permutation null\n{CORRELATION_PERMUTATIONS} permutations; n=60/time")
    ax.legend(fontsize=5.5)
    add_panel_label(ax, "B")
    export_figure(fig, FIGURE_DIR / "Reviewer_compositionality_CLR_focal_exclusion")


def write_key_results(
    primary: BranchResult,
    hc3: pd.DataFrame,
    gee_omnibus: pd.DataFrame,
    distance: pd.DataFrame,
    correlations: pd.DataFrame,
    overlap: pd.DataFrame,
) -> None:
    model = hc3.loc[
        (hc3["dataset_branch"] == primary.name)
        & (hc3["model"] == "adjusted_BP_continuous")
        & (hc3["term"] == "z_baseline_shannon")
    ].set_index("outcome")
    subject = primary.subject_metrics
    original_overlap = overlap.loc[
        (overlap["dataset_branch"] == primary.name)
        & (overlap["scheme"] == "metric_extreme_n15")
        & (overlap["new_group_definition"] == PRIMARY_RECOVERY_METRIC)
    ].iloc[0]
    group_gee = gee_omnibus.loc[
        (gee_omnibus["dataset_branch"] == primary.name)
        & (gee_omnibus["outcome"] == "signed_shannon_deviation_from_baseline")
        & (gee_omnibus["model"] == "original_HDR_LDR")
    ].iloc[0]
    paired_distance = distance.loc[
        (distance["dataset_branch"] == primary.name)
        & (distance["distance"] == "Bray-Curtis")
        & (distance["comparison"] == "all60_within_subject_phase")
    ].iloc[0]
    corrected_significant = correlations.loc[
        (correlations["dataset_branch"] == primary.name) & (correlations["q_BH_global_16_times"] < 0.05)
    ]
    lines = [
        "# Core reviewer reanalysis: response-ready numerical summary",
        "",
        "These statements are deliberately conservative and should be checked against the manuscript wording before insertion.",
        "",
        "## R1 major 1–3 / R2.3: recovery definition and full-cohort models",
        "",
        f"- The primary analysis retained all {len(subject)} ciprofloxacin-exposed participants and defined recovery as C4–C6 Bray–Curtis distance to each participant's own A1–A4 compositional centroid (lower is closer recovery).",
    ]
    for outcome, label in [
        (PRIMARY_RECOVERY_METRIC, "endpoint Bray–Curtis recovery distance"),
        ("endpoint_abs_shannon_deviation", "absolute endpoint Shannon deviation"),
        ("recovery_auc_bray_actual_day", "actual-day recovery Bray–Curtis AUC"),
    ]:
        row = model.loc[outcome]
        lines.append(
            f"- After adjustment for baseline Bacteroides + Phocaeicola abundance, a 1-SD higher baseline Shannon was associated with {label}: standardized beta={row['estimate_standardized_outcome']:.3f}, 95% CI [{row['CI95_low']:.3f}, {row['CI95_high']:.3f}], HC3 p={row['p_value']:.4g}, n={int(row['n_subjects'])}."
        )
    lines += [
        f"- For n=15 recovery-defined extremes, the high-recovery set overlapped the original HDR set in {int(original_overlap['high_recovery_or_endpoint_overlap_original_HDR'])}/15 participants and the low-recovery set overlapped original LDR in {int(original_overlap['low_recovery_or_endpoint_overlap_original_LDR'])}/15. This is overlap, not validation of the original endpoint-based labels.",
        f"- Persistent Shannon return by C6 occurred in {int(subject['persistent_shannon_return_by_C6'].sum())}/{len(subject)} participants under the prespecified personal baseline mean +/-1 SD rule; non-returners were censored descriptively at their last observed study day.",
        "",
        "## R1 major 6: repeated measures and distance inference",
        "",
        f"- The original endpoint-defined HDR/LDR Shannon trajectories were tested with participant-clustered exchangeable GEE: joint group-by-time Wald chi-square={group_gee['wald_chi2']:.3f}, df={int(group_gee['df'])}, p={group_gee['p_value']:.4g}. Because HDR/LDR was selected using endpoint Shannon, this result remains selection-coupled/descriptive.",
        f"- In all 60 participants, the paired baseline-versus-endpoint Bray–Curtis test used phase-label swaps within participant: PERMANOVA pseudo-F={paired_distance['PERMANOVA_pseudoF']:.3f}, R2={paired_distance['PERMANOVA_R2']:.3f}, p={paired_distance['PERMANOVA_p']:.4g}; accompanying PERMDISP F={paired_distance['PERMDISP_F']:.3f}, p={paired_distance['PERMDISP_p']:.4g}.",
        "",
        "## R1 major 7: compositionality",
        "",
        f"- CLR genus abundance was correlated with Shannon recomputed after removing all species assigned to the focal genus; inference used {CORRELATION_PERMUTATIONS} participant-label permutations per genus/time and global BH correction across 16 times. {len(corrected_significant)} genus-time associations met global q<0.05 in the rarefied primary branch.",
        "",
        "## Required limitations",
        "",
        "- The primary rarefied branch follows the Methods statement (2,693,000 reads; seed 20261009); the unrarefied branch is retained because the published figure coordinates are more consistent with unrarefied leaf counts. This Methods/figure conflict must be disclosed.",
        "- The exploratory k=3 CLR enterotype is a sensitivity covariate, not a validated biological enterotype. Bacteroides + Phocaeicola abundance/dominance is the more transparent adjustment.",
        "- Age was unavailable for 15/60 and sex for 2/60 exposed participants, so neither was silently included in the primary complete-case model.",
        "- The five unexposed controls are not used for inferential claims; n=5 permits descriptive display only.",
    ]
    (ANALYSIS_DIR / "core_key_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def branch_sensitivity_summary(subject_table: pd.DataFrame, hc3: pd.DataFrame) -> pd.DataFrame:
    """Summarize agreement between Methods-faithful and source-figure branches."""
    rows: list[dict[str, object]] = []
    primary = subject_table.loc[subject_table["dataset_branch"] == "rarefied_primary"].set_index("subject")
    sensitivity = subject_table.loc[subject_table["dataset_branch"] == "unrarefied_sensitivity"].set_index("subject")
    for outcome in CONTINUOUS_OUTCOMES + ["baseline_shannon", "endpoint_shannon", "recovery_composite_mean_rank"]:
        aligned = primary[[outcome]].join(sensitivity[[outcome]], lsuffix="_primary", rsuffix="_sensitivity")
        rho = spearmanr(aligned.iloc[:, 0], aligned.iloc[:, 1])
        rows.append(
            {
                "comparison_type": "subject_metric",
                "outcome": outcome,
                "model": "not_applicable",
                "term": "not_applicable",
                "primary_value": float(aligned.iloc[:, 0].mean()),
                "sensitivity_value": float(aligned.iloc[:, 1].mean()),
                "sensitivity_minus_primary": float((aligned.iloc[:, 1] - aligned.iloc[:, 0]).mean()),
                "subject_level_spearman_rho": float(rho.statistic),
                "subject_level_spearman_p": float(rho.pvalue),
                "n_subjects": len(aligned),
            }
        )
    selected = hc3.loc[
        (hc3["model"] == "adjusted_BP_continuous") & (hc3["term"] == "z_baseline_shannon")
    ]
    for outcome in CONTINUOUS_OUTCOMES:
        part = selected.loc[selected["outcome"] == outcome].set_index("dataset_branch")
        primary_value = float(part.loc["rarefied_primary", "estimate_standardized_outcome"])
        sensitivity_value = float(part.loc["unrarefied_sensitivity", "estimate_standardized_outcome"])
        rows.append(
            {
                "comparison_type": "HC3_standardized_beta",
                "outcome": outcome,
                "model": "adjusted_BP_continuous",
                "term": "z_baseline_shannon",
                "primary_value": primary_value,
                "sensitivity_value": sensitivity_value,
                "sensitivity_minus_primary": sensitivity_value - primary_value,
                "subject_level_spearman_rho": np.nan,
                "subject_level_spearman_p": np.nan,
                "n_subjects": int(part.loc["rarefied_primary", "n_subjects"]),
            }
        )
    return pd.DataFrame(rows)


def write_readme_and_manifest(
    meta: pd.DataFrame,
    species_counts: pd.DataFrame,
    branches: list[BranchResult],
    outputs: list[Path],
) -> None:
    exposed = meta.loc[meta["abx"] == "Y"]
    manifest = {
        "analysis": "reviewer-focused longitudinal recovery reanalysis",
        "created": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "seed": SEED,
        "primary_branch": "rarefied_primary",
        "rarefaction_depth": RAREFACTION_DEPTH,
        "rarefaction_scope": "full 1040-sample bacterial species-leaf matrix before selecting exposed participants",
        "sensitivity_branch": "unrarefied_sensitivity",
        "methods_figure_conflict": "Methods states rarefaction; source-raster calibration supports unrarefied figure inputs",
        "primary_recovery_outcome": PRIMARY_RECOVERY_METRIC,
        "primary_recovery_direction": "lower is closer return to personal A1-A4 baseline composition",
        "participants_exposed": int(exposed["subject"].nunique()),
        "participants_unexposed_descriptive_only": int(meta.loc[meta["abx"] == "N", "subject"].nunique()),
        "samples_exposed": int(len(exposed)),
        "samples_per_exposed_participant": sorted(exposed.groupby("subject").size().unique().tolist()),
        "bacterial_species_leaf_features": int(species_counts.shape[0]),
        "raw_input": str(RAW_DIR / "combine_brk.mpa"),
        "raw_input_sha256": sha256(RAW_DIR / "combine_brk.mpa"),
        "metadata_sha256": sha256(RAW_DIR / "meta.csv"),
        "analysis_script_sha256": sha256(ANALYSIS_DIR / "core_reanalysis.py"),
        "unit_test_sha256": sha256(ANALYSIS_DIR / "test_core_reanalysis.py"),
        "permanova_permutations": PERMANOVA_PERMUTATIONS,
        "compositional_null_permutations": CORRELATION_PERMUTATIONS,
        "branch_enterotype_silhouette": {branch.name: branch.enterotype_silhouette for branch in branches},
        "host_covariate_availability_exposed": {
            "age_subjects": int(exposed.groupby("subject")["age"].first().notna().sum()),
            "sex_subjects": int(exposed.groupby("subject")["sex"].first().notna().sum()),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": sm.__version__,
            "matplotlib": mpl.__version__,
            "scikit_learn": importlib.metadata.version("scikit-learn"),
            "Pillow": importlib.metadata.version("Pillow"),
        },
        "outputs": [str(path.relative_to(ROOT)) for path in outputs],
    }
    (TABLE_DIR / "core_reanalysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    readme = f"""# Core reviewer reanalysis

This folder contains a reproducible Python analysis responding to Reviewer 1 major comments 1–3, 6, and 7 and Reviewer 2 comment 3. It does not modify the original reproduction pipeline.

## Figure contract

Core conclusion: recovery must be evaluated as return toward each participant's own baseline, not as absolute endpoint diversity. The primary outcome is the participant-level mean C4–C6 Bray–Curtis distance to that participant's A1–A4 species-composition centroid. The evidence layout is a quantitative grid: full-cohort continuous models, cutoff/overlap sensitivity, repeated-measures models, valid distance permutations with PERMDISP, and compositionality-aware correlations.

## Run

```bash
MPLCONFIGDIR=/tmp/mpl-reviewer-core python reviewer_response/02_analysis/core_reanalysis.py
python -m pytest -q reviewer_response/02_analysis/test_core_reanalysis.py
```

The run is deterministic (`seed={SEED}`). The primary input is the bacterial species-leaf portion of `raw/combine_brk.mpa`; metadata are read from `raw/meta.csv`.

## Prespecified analysis decisions

- Primary branch: rarefy the full 1,040-sample species-leaf table without replacement to {RAREFACTION_DEPTH:,} reads, then retain all 60 exposed participants (960 samples).
- Sensitivity branch: unrarefied species-leaf counts. This is necessary because the manuscript Methods specify rarefaction, while source-figure coordinate calibration supports unrarefied inputs. Results must identify the branch.
- Personal baseline: A1–A4. Endpoint: C4–C6. Recovery interval for actual-day AUC: B6–C6.
- Primary recovery outcome: mean endpoint Bray–Curtis distance to the personal baseline composition centroid (lower is better recovery).
- Other recovery outcomes: absolute Shannon displacement; binary Jaccard to the personal baseline consensus (presence in at least 2/4 baseline samples); normalized actual-day AUCs; and time to first/persistent return within personal baseline mean ±1 baseline SD.
- Full-cohort models use all 60 participants and HC3 covariance. Baseline Bacteroides + Phocaeicola abundance/dominance and an explicitly exploratory k=3 CLR enterotype are adjustment sensitivities.
- Repeated measures use exchangeable GEE clustered by participant; random-intercept MixedLM is a sensitivity analysis. Original HDR/LDR tests are labelled selection-coupled/descriptive.
- PERMANOVA uses one centroid per participant/stage and participant-label permutations; paired baseline/end tests swap phase within participant. PERMDISP is always reported alongside PERMANOVA.
- Taxon–diversity correlations use genus CLR abundance and Shannon recomputed after excluding all species in the focal genus. P values use participant-label permutations and BH correction within time and globally across 16 times.

## Interpretive boundaries

- Original HDR/LDR labels were reconstructed from endpoint Shannon and cannot independently validate endpoint-Shannon or related trajectory contrasts.
- Recovery-defined extreme groups are useful illustrations but remain outcome-derived; the continuous n=60 models are primary.
- The k=3 enterotype is an exploratory clustering covariate, not a validated enterotype taxonomy.
- Age is available for only 45/60 exposed participants and sex for 58/60. Primary models avoid silent complete-case deletion.
- The five unexposed participants are descriptive only.
- A nonsignificant result is not evidence of equivalence, and a PERMANOVA association with significant PERMDISP cannot be interpreted solely as a centroid shift.

## Outputs

Tables are written to `reviewer_response/03_tables/`; figures (editable SVG, PDF, and LZW-compressed 600-dpi TIFF) to `reviewer_response/04_figures/`. `core_key_results.md` provides response-ready numerical sentences but retains all caveats.

## Method resource

Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065
"""
    (ANALYSIS_DIR / "README_core_reanalysis.md").write_text(readme, encoding="utf-8")


def validate_outputs(outputs: list[Path]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in outputs:
        exists = path.exists() and path.stat().st_size > 0
        record: dict[str, object] = {"path": str(path.relative_to(ROOT)), "exists_nonempty": exists}
        if exists and path.suffix.lower() in {".tif", ".tiff"}:
            with Image.open(path) as image:
                record.update(
                    {
                        "format": image.format,
                        "mode": image.mode,
                        "width_px": image.width,
                        "height_px": image.height,
                        "dpi_x": image.info.get("dpi", (np.nan, np.nan))[0],
                        "dpi_y": image.info.get("dpi", (np.nan, np.nan))[1],
                        "compression": image.info.get("compression", "unknown"),
                    }
                )
        rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-correlations", action="store_true", help="Development-only speed option; not for final delivery")
    args = parser.parse_args()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(ROOT / "reproduction"))
    import load_raw_bracken as raw_loader  # reuse the audited grouping/rarefaction implementation

    meta = pd.read_csv(RAW_DIR / "meta.csv", dtype={"ids": str, "subject": str, "canon": str})
    mpa = pd.read_csv(RAW_DIR / "combine_brk.mpa", sep="\t", index_col=0)
    mpa.columns = mpa.columns.str.removesuffix(".brk.report")
    if set(meta["ids"]) != set(mpa.columns):
        raise ValueError("raw metadata/sample columns differ")
    sample_order = meta["ids"].tolist()
    species_counts = mpa.loc[species_leaf_mask(mpa.index), sample_order].astype(np.int64)
    if species_counts.shape[0] != 11_773:
        warnings.warn(f"Expected 11,773 bacterial species leaves, observed {species_counts.shape[0]}")
    if int(species_counts.sum(axis=0).min()) < RAREFACTION_DEPTH:
        raise ValueError("At least one sample is below the rarefaction depth")
    species_to_genus = pd.Series([genus_from_path(path) for path in species_counts.index], index=species_counts.index, name="genus")
    original_groups = raw_loader._legacy_grouping(ROOT, meta)

    rarefied_counts = raw_loader._rarefy(species_counts, RAREFACTION_DEPTH, SEED)
    branches = [
        build_branch("rarefied_primary", rarefied_counts, meta, original_groups, species_to_genus),
        build_branch("unrarefied_sensitivity", species_counts, meta, original_groups, species_to_genus),
    ]
    del mpa, rarefied_counts

    subject_table = pd.concat([branch.subject_metrics.reset_index() for branch in branches], ignore_index=True)
    sample_table = pd.concat([branch.sample_metrics for branch in branches], ignore_index=True)
    hc3_tables, binary_tables, gee_coef_tables, gee_omnibus_tables, mixed_tables = [], [], [], [], []
    cutoff_tables, overlap_tables, distance_tables, coordinate_tables = [], [], [], []
    correlation_tables = []
    for branch in branches:
        hc3, binary = fit_hc3_models(branch)
        gee_coef, gee_omnibus, mixed = fit_repeated_models(branch)
        cutoff, overlap = cutoff_and_overlap(branch)
        distance, coords = distance_tests(branch)
        hc3_tables.append(hc3)
        binary_tables.append(binary)
        gee_coef_tables.append(gee_coef)
        gee_omnibus_tables.append(gee_omnibus)
        mixed_tables.append(mixed)
        cutoff_tables.append(cutoff)
        overlap_tables.append(overlap)
        distance_tables.append(distance)
        coordinate_tables.append(coords)
        if not args.skip_correlations:
            correlation_tables.append(compositional_correlations(branch))

    hc3 = pd.concat(hc3_tables, ignore_index=True)
    binary = pd.concat(binary_tables, ignore_index=True)
    gee_coef = pd.concat(gee_coef_tables, ignore_index=True)
    gee_omnibus = pd.concat(gee_omnibus_tables, ignore_index=True)
    mixed = pd.concat(mixed_tables, ignore_index=True)
    cutoff = pd.concat(cutoff_tables, ignore_index=True)
    overlap = pd.concat(overlap_tables, ignore_index=True)
    distance = pd.concat(distance_tables, ignore_index=True)
    coordinates = pd.concat(coordinate_tables, ignore_index=True)
    correlations = pd.concat(correlation_tables, ignore_index=True) if correlation_tables else pd.DataFrame()

    tables = {
        "core_subject_recovery_metrics.csv": subject_table,
        "core_sample_longitudinal_metrics.csv": sample_table,
        "core_continuous_models_HC3.csv": hc3,
        "core_return_to_baseline_logistic_HC3.csv": binary,
        "core_repeated_GEE_coefficients.csv": gee_coef,
        "core_repeated_GEE_omnibus.csv": gee_omnibus,
        "core_repeated_MixedLM_sensitivity.csv": mixed,
        "core_extreme_cutoff_sensitivity.csv": cutoff,
        "core_original_vs_recovery_group_overlap.csv": overlap,
        "core_distance_PERMANOVA_PERMDISP.csv": distance,
        "core_distance_PCoA_coordinates.csv": coordinates,
        "core_rarefied_unrarefied_sensitivity_summary.csv": branch_sensitivity_summary(subject_table, hc3),
    }
    if not correlations.empty:
        tables["core_compositional_CLR_focal_excluded_correlations.csv"] = correlations
    for filename, frame in tables.items():
        frame.to_csv(TABLE_DIR / filename, index=False)

    set_style()
    render_core_figure(branches[0], hc3, cutoff, overlap, gee_omnibus)
    render_distance_figure(branches[0], distance, coordinates)
    if not correlations.empty:
        render_compositional_figure(correlations)
        write_key_results(branches[0], hc3, gee_omnibus, distance, correlations, overlap)

    output_paths = [TABLE_DIR / filename for filename in tables]
    for stem in [
        "Reviewer_core_recovery_reanalysis",
        "Reviewer_distance_PERMANOVA_PERMDISP",
        "Reviewer_compositionality_CLR_focal_exclusion",
    ]:
        for suffix in [".pdf", ".svg", ".tif"]:
            path = FIGURE_DIR / f"{stem}{suffix}"
            if path.exists() or stem != "Reviewer_compositionality_CLR_focal_exclusion" or not args.skip_correlations:
                output_paths.append(path)
    output_paths.extend(
        [
            ANALYSIS_DIR / "core_reanalysis.py",
            ANALYSIS_DIR / "test_core_reanalysis.py",
            ANALYSIS_DIR / "README_core_reanalysis.md",
            ANALYSIS_DIR / "core_key_results.md",
            ROOT / "reviewer_response" / "05_response" / "core_reanalysis_response_blocks.md",
        ]
    )
    write_readme_and_manifest(meta, species_counts, branches, output_paths)
    output_paths.append(TABLE_DIR / "core_reanalysis_manifest.json")
    qa = validate_outputs(output_paths)
    qa.to_csv(TABLE_DIR / "core_output_QA.tsv", sep="\t", index=False)
    if not qa["exists_nonempty"].all():
        missing = qa.loc[~qa["exists_nonempty"], "path"].tolist()
        raise RuntimeError(f"Missing/empty outputs: {missing}")
    print(json.dumps({"status": "complete", "outputs": len(output_paths), "primary_n": 60}, indent=2))


if __name__ == "__main__":
    main()
