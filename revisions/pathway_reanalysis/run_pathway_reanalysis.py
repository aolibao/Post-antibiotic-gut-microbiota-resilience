#!/usr/bin/env python3
"""Reanalyse three inferred human microbial pathways at the subject level.

The pathway inputs contain baseline samples only (A1--A4).  Consequently,
all inferential rows in this workflow are first aggregated to one row per
subject.  Endpoint-defined HDR/LDR contrasts are retained only as explicitly
selection-coupled descriptions.  The primary analysis relates baseline
inferred functional potential to a continuous, within-subject compositional
return metric in all 60 ciprofloxacin-exposed participants.

This script deliberately uses the phrase ``inferred functional potential``.
CPM values inferred from metagenomic gene content do not measure NAD
concentrations, pathway activity, metabolic flux, or a causal mechanism.
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
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
from matplotlib.lines import Line2D
from PIL import Image
from scipy.stats import norm, rankdata
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning


SEED = 20261009
BOOTSTRAPS = 5000
BASELINE = ("A1", "A2", "A3", "A4")
ENDPOINT = ("C4", "C5", "C6")
PATHWAYS: Mapping[str, str] = {
    "PWY-7761": "NAD salvage pathway II\n(PNC IV cycle)",
    "PYRIDNUCSYN-PWY": "NAD de novo biosynthesis I\n(from aspartate)",
    "PWY-5265": "Peptidoglycan biosynthesis II\n(staphylococci)",
}
PATHWAY_SHORT: Mapping[str, str] = {
    "PWY-7761": "NAD salvage II",
    "PYRIDNUCSYN-PWY": "NAD de novo I",
    "PWY-5265": "Peptidoglycan II",
}
PATHWAY_COLORS: Mapping[str, str] = {
    "PWY-7761": "#4477AA",
    "PYRIDNUCSYN-PWY": "#AA6F3C",
    "PWY-5265": "#6B8E6E",
}
NAD_GENE_KOS: Mapping[str, tuple[str, ...]] = {
    "nadD": ("K00969",),
    "nudC": ("K03426",),
    # KEGG has two bacterial pncC orthology entries for the same enzymatic
    # activity; both are summed to audit the manuscript's gene-level claim.
    "pncC": ("K03742", "K03743"),
}
GROUP_COLORS = {"LDR": "#3B6EA8", "HDR": "#C7614A", "Middle": "#A4A4A4"}


@dataclass(frozen=True)
class Paths:
    project_root: Path
    tables: Path
    figures: Path

    @property
    def raw(self) -> Path:
        return self.project_root / "raw"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    """Benjamini--Hochberg adjustment with NaN preservation."""
    p = np.asarray(list(values), dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return adjusted
    pv = p[valid]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.clip(q, 0.0, 1.0)
    adjusted[valid] = restored
    return adjusted


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    sd = values.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"Cannot standardize degenerate variable {series.name!r}")
    return (values - values.mean()) / sd


def rank_normal(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    out = pd.Series(np.nan, index=values.index, dtype=float)
    mask = values.notna()
    n = int(mask.sum())
    if n:
        ranks = rankdata(values.loc[mask], method="average")
        out.loc[mask] = norm.ppf((ranks - 0.5) / n)
    return out


def _pathway_row(index: pd.Index, pathway: str) -> str:
    matches = [str(x) for x in index if str(x) == pathway or str(x).startswith(pathway + ":")]
    if len(matches) != 1:
        raise ValueError(f"Expected one unstratified row for {pathway}; found {len(matches)}")
    return matches[0]


def load_pathways(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return metadata, target unstratified CPM and all target stratified rows."""
    meta = pd.read_csv(paths.raw / "meta.csv", dtype=str)
    unstrat = pd.read_csv(paths.raw / "pathway_cpm_unstratified.tsv", sep="\t", index_col=0)
    strat = pd.read_csv(paths.raw / "pathway_cpm_stratified.tsv", sep="\t", index_col=0)
    for frame in (unstrat, strat):
        frame.columns = frame.columns.astype(str).str.removesuffix("_Abundance-CPM")
    selected_rows = {_id: _pathway_row(unstrat.index, _id) for _id in PATHWAYS}
    target_unstrat = pd.DataFrame(
        {pathway: unstrat.loc[row].astype(float) for pathway, row in selected_rows.items()}
    )
    keep = strat.index.to_series().astype(str).map(
        lambda value: any(value.startswith(pathway + ":") or value.startswith(pathway + "|") for pathway in PATHWAYS)
    )
    target_strat = strat.loc[keep.to_numpy()].astype(float)
    if not np.isfinite(target_unstrat.to_numpy()).all() or (target_unstrat.to_numpy() < 0).any():
        raise ValueError("Target unstratified pathway CPM contains invalid values")
    if not np.isfinite(target_strat.to_numpy()).all() or (target_strat.to_numpy() < 0).any():
        raise ValueError("Target stratified pathway CPM contains invalid values")
    if set(target_unstrat.index) - set(meta["ids"]):
        raise ValueError("Pathway samples are absent from metadata")
    baseline_meta = meta.set_index("ids").loc[target_unstrat.index].reset_index()
    if "ids" not in baseline_meta.columns and "index" in baseline_meta.columns:
        baseline_meta = baseline_meta.rename(columns={"index": "ids"})
    if set(baseline_meta["canon"]) != set(BASELINE):
        raise ValueError("Pathway table is not restricted to A1--A4")
    counts = baseline_meta.groupby("subject").size()
    if baseline_meta["subject"].nunique() != 60 or not (counts == 4).all():
        raise ValueError("Expected 60 exposed subjects with four baseline pathway samples each")
    if not (baseline_meta["abx"] == "Y").all():
        raise ValueError("Unexpected non-exposed participant in pathway table")
    return baseline_meta, target_unstrat, target_strat


def _rank_mask(index: pd.Index, prefix: str) -> np.ndarray:
    strings = index.to_series(index=index).astype(str)
    terminal = strings.str.rsplit("|", n=1).str[-1]
    return (strings.str.startswith("k__Bacteria|") & terminal.str.startswith(prefix)).to_numpy()


def _terminal_names(index: pd.Index, prefix: str) -> pd.Index:
    names = (
        index.to_series(index=index)
        .astype(str)
        .str.rsplit("|", n=1)
        .str[-1]
        .str.removeprefix(prefix)
        .str.replace("_", " ", regex=False)
        .str.strip()
        .replace("", "Unclassified")
    )
    return pd.Index(names.to_numpy())


def _legacy_endpoint_groups(project_root: Path, meta: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct manuscript HDR/LDR labels exactly as the existing loader does."""
    legacy = pd.read_csv(project_root / "species.csv", index_col=0)
    exposed_ids = meta.loc[meta["abx"] == "Y", "ids"].tolist()
    legacy = legacy.loc[:, exposed_ids]
    prop = legacy.div(legacy.sum(axis=0).replace(0, np.nan), axis=1)
    shannon = -(prop * np.log(prop.where(prop > 0))).sum(axis=0)
    frame = meta.loc[meta["abx"] == "Y", ["ids", "subject", "canon"]].copy()
    frame["legacy_shannon"] = frame["ids"].map(shannon)
    endpoint = (
        frame.loc[frame["canon"].isin(ENDPOINT)]
        .groupby("subject")["legacy_shannon"]
        .mean()
        .sort_values(kind="mergesort")
    )
    baseline = frame.loc[frame["canon"].isin(BASELINE)].groupby("subject")["legacy_shannon"].mean()
    out = pd.DataFrame({"legacy_baseline_shannon": baseline, "legacy_endpoint_shannon": endpoint})
    out["hdr_ldr_group"] = "Middle"
    out.loc[endpoint.index[:15], "hdr_ldr_group"] = "LDR"
    out.loc[endpoint.index[-15:], "hdr_ldr_group"] = "HDR"
    out["endpoint_rank"] = endpoint.rank(method="first").astype(int)
    out.index.name = "subject"
    return out


def load_taxonomy_and_recovery(paths: Paths, meta: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    """Build baseline confounders and independent within-subject recovery metrics."""
    mpa = pd.read_csv(paths.raw / "combine_brk.mpa", sep="\t", index_col=0)
    mpa.columns = mpa.columns.astype(str).str.removesuffix(".brk.report")
    sample_order = meta["ids"].tolist()
    if set(sample_order) != set(mpa.columns):
        raise ValueError("Metadata and Bracken sample IDs differ")
    species = mpa.loc[_rank_mask(mpa.index, "s__"), sample_order].astype(float)
    genus = mpa.loc[_rank_mask(mpa.index, "g__"), sample_order].astype(float)
    if (species.to_numpy() < 0).any() or (genus.to_numpy() < 0).any():
        raise ValueError("Negative Bracken counts")

    species_depth = species.sum(axis=0)
    species_prop = species.div(species_depth.replace(0, np.nan), axis=1)
    richness = (species > 0).sum(axis=0).astype(float)
    genus_prop = genus.div(genus.sum(axis=0).replace(0, np.nan), axis=1)
    genus_named = genus_prop.groupby(_terminal_names(genus.index, "g__"), sort=False).sum()
    for genus_name in ("Bacteroides", "Phocaeicola"):
        if genus_name not in genus_named.index:
            genus_named.loc[genus_name] = 0.0

    baseline = meta.loc[(meta["abx"] == "Y") & meta["canon"].isin(BASELINE), ["ids", "subject"]].copy()
    baseline["bacteroides_fraction"] = baseline["ids"].map(genus_named.loc["Bacteroides"])
    baseline["phocaeicola_fraction"] = baseline["ids"].map(genus_named.loc["Phocaeicola"])
    baseline["richness"] = baseline["ids"].map(richness)
    baseline["bacterial_reads"] = baseline["ids"].map(species_depth)
    covariates = baseline.groupby("subject", sort=True).agg(
        bacteroides_fraction=("bacteroides_fraction", "mean"),
        phocaeicola_fraction=("phocaeicola_fraction", "mean"),
        richness=("richness", "mean"),
        bacterial_reads=("bacterial_reads", "mean"),
    )
    covariates["bact_phoc_fraction"] = covariates["bacteroides_fraction"] + covariates["phocaeicola_fraction"]

    # Endpoint similarity to each participant's A1--A4 baseline centroid.
    recovery_rows: list[dict[str, float | str]] = []
    sample_to_col = {sample: i for i, sample in enumerate(species_prop.columns)}
    presence = species.to_numpy() > 0
    for subject, rows in meta.loc[meta["abx"] == "Y"].groupby("subject", sort=True):
        baseline_ids = rows.loc[rows["canon"].isin(BASELINE), "ids"].tolist()
        endpoint_ids = rows.loc[rows["canon"].isin(ENDPOINT), "ids"].tolist()
        if len(baseline_ids) != 4 or len(endpoint_ids) != 3:
            raise ValueError(f"Unexpected longitudinal design for {subject}")
        centroid = species_prop[baseline_ids].mean(axis=1).to_numpy()
        baseline_binary = [presence[:, sample_to_col[sample]] for sample in baseline_ids]
        bray: list[float] = []
        jaccard: list[float] = []
        for sample in endpoint_ids:
            current = species_prop[sample].to_numpy()
            denom = np.abs(current + centroid).sum()
            bray.append(0.0 if denom == 0 else float(np.abs(current - centroid).sum() / denom))
            current_binary = presence[:, sample_to_col[sample]]
            distances = []
            for base in baseline_binary:
                union = np.logical_or(base, current_binary).sum()
                distances.append(0.0 if union == 0 else 1.0 - np.logical_and(base, current_binary).sum() / union)
            jaccard.append(float(np.mean(distances)))
        recovery_rows.append(
            {
                "subject": subject,
                "endpoint_bray_distance": float(np.mean(bray)),
                "endpoint_jaccard_distance": float(np.mean(jaccard)),
                "bray_recovery": 1.0 - float(np.mean(bray)),
                "jaccard_recovery": 1.0 - float(np.mean(jaccard)),
            }
        )
    recovery = pd.DataFrame(recovery_rows).set_index("subject")
    recovery = recovery.rename(
        columns={
            "endpoint_bray_distance": "endpoint_bray_distance_unrarefied_direct",
            "endpoint_jaccard_distance": "endpoint_jaccard_distance_unrarefied_direct",
            "bray_recovery": "bray_recovery_unrarefied_direct",
            "jaccard_recovery": "jaccard_recovery_unrarefied_direct",
        }
    )
    out = covariates.join(recovery, how="inner").join(_legacy_endpoint_groups(paths.project_root, meta), how="inner")

    # Align the pathway audit to the formal reviewer reanalysis: the Methods-
    # faithful rarefied branch defines the primary recovery outcome, while the
    # unrarefied calculation above remains an explicit sensitivity/audit field.
    core_path = paths.project_root / "reviewer_response" / "03_tables" / "core_subject_recovery_metrics.csv"
    if not core_path.exists():
        raise FileNotFoundError(
            "Run reviewer_response/02_analysis/core_reanalysis.py first; "
            "core_subject_recovery_metrics.csv is required for the rarefied primary recovery outcome"
        )
    core_all = pd.read_csv(core_path, dtype={"subject": str, "dataset_branch": str})
    core = core_all.loc[core_all["dataset_branch"] == "rarefied_primary"].set_index("subject")
    core_columns = [
        "endpoint_mean_bray_to_baseline_centroid",
        "endpoint_mean_jaccard_to_baseline_consensus",
        "recovery_auc_bray_actual_day",
        "recovery_auc_jaccard_actual_day",
        "baseline_richness",
        "baseline_bacteroides_fraction",
        "baseline_phocaeicola_fraction",
        "baseline_bacteroides_plus_phocaeicola_fraction",
        "original_endpoint_group",
        "original_endpoint_rank",
    ]
    if core.shape[0] != 60 or any(column not in core.columns for column in core_columns):
        raise ValueError("Formal core recovery table lacks the expected 60-subject rarefied branch")
    core = core[core_columns].rename(
        columns={
            "endpoint_mean_bray_to_baseline_centroid": "endpoint_bray_distance_rarefied_primary",
            "endpoint_mean_jaccard_to_baseline_consensus": "endpoint_jaccard_distance_rarefied_primary",
            "recovery_auc_bray_actual_day": "bray_auc_distance_rarefied_primary",
            "recovery_auc_jaccard_actual_day": "jaccard_auc_distance_rarefied_primary",
            "baseline_richness": "richness_rarefied_primary",
            "baseline_bacteroides_fraction": "bacteroides_fraction_rarefied_primary",
            "baseline_phocaeicola_fraction": "phocaeicola_fraction_rarefied_primary",
            "baseline_bacteroides_plus_phocaeicola_fraction": "bact_phoc_fraction_rarefied_primary",
            "original_endpoint_group": "core_original_endpoint_group",
            "original_endpoint_rank": "core_original_endpoint_rank",
        }
    )
    out = out.join(core, how="inner")
    if out.shape[0] != 60 or out.isna().any().any():
        raise ValueError("Subject covariate/recovery table is incomplete")
    out["bray_recovery"] = 1.0 - out["endpoint_bray_distance_rarefied_primary"]
    out["jaccard_recovery"] = 1.0 - out["endpoint_jaccard_distance_rarefied_primary"]
    out["bray_auc_recovery"] = 1.0 - out["bray_auc_distance_rarefied_primary"]
    out["jaccard_auc_recovery"] = 1.0 - out["jaccard_auc_distance_rarefied_primary"]
    out = out.rename(
        columns={
            "bacteroides_fraction": "bacteroides_fraction_unrarefied",
            "phocaeicola_fraction": "phocaeicola_fraction_unrarefied",
            "bact_phoc_fraction": "bact_phoc_fraction_unrarefied",
            "richness": "richness_unrarefied",
        }
    )
    out["bacteroides_fraction"] = out["bacteroides_fraction_rarefied_primary"]
    out["phocaeicola_fraction"] = out["phocaeicola_fraction_rarefied_primary"]
    out["bact_phoc_fraction"] = out["bact_phoc_fraction_rarefied_primary"]
    out["richness"] = out["richness_rarefied_primary"]
    out["bray_recovery_z"] = zscore(out["bray_recovery"])
    out["jaccard_recovery_z"] = zscore(out["jaccard_recovery"])
    out["bray_auc_recovery_z"] = zscore(out["bray_auc_recovery"])
    out["jaccard_auc_recovery_z"] = zscore(out["jaccard_auc_recovery"])
    out["bact_z"] = zscore(np.log1p(out["bacteroides_fraction"]).rename("bact_z"))
    out["phoc_z"] = zscore(np.log1p(out["phocaeicola_fraction"]).rename("phoc_z"))
    out["bp_z"] = zscore(np.log1p(out["bact_phoc_fraction"]).rename("bp_z"))
    out["richness_z"] = zscore(out["richness"].rename("richness_z"))
    out["log_depth_z"] = zscore(np.log(out["bacterial_reads"]).rename("log_depth_z"))
    out["bp_unrarefied_z"] = zscore(
        np.log1p(out["bact_phoc_fraction_unrarefied"]).rename("bp_unrarefied_z")
    )
    out["richness_unrarefied_z"] = zscore(out["richness_unrarefied"].rename("richness_unrarefied_z"))

    # Species epithets whose current Bracken genus is Phocaeicola.  HUMAnN's
    # stratified labels retain older Bacteroides names, so this set supports a
    # transparent taxonomy-harmonized split while leaving the combined share unchanged.
    phoc_epithets: set[str] = set()
    for lineage in species.index.astype(str):
        parts = lineage.split("|")
        if "g__Phocaeicola" not in parts:
            continue
        terminal = parts[-1].removeprefix("s__")
        if terminal.startswith("Phocaeicola_"):
            phoc_epithets.add(terminal.removeprefix("Phocaeicola_"))
    return out, phoc_epithets


def _parse_stratified_taxon(row_name: str, phoc_epithets: set[str]) -> tuple[str, str, str, str]:
    taxon = row_name.split("|", 1)[1] if "|" in row_name else "unclassified"
    if taxon == "unclassified":
        return taxon, "unclassified", "unclassified", "unclassified row"
    genus = "unclassified"
    species = ""
    for part in taxon.split("."):
        if part.startswith("g__"):
            genus = part.removeprefix("g__")
        elif part.startswith("s__"):
            species = part.removeprefix("s__")
    harmonized = genus
    reason = "literal genus"
    if genus == "Bacteroides" and species.startswith("Bacteroides_"):
        epithet = species.removeprefix("Bacteroides_")
        if epithet in phoc_epithets:
            harmonized = "Phocaeicola"
            reason = "species-epithet crosswalk to current Bracken genus"
    return taxon, genus, harmonized, reason


def aggregate_pathways(
    baseline_meta: pd.DataFrame,
    unstrat: pd.DataFrame,
    strat: pd.DataFrame,
    phoc_epithets: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_subject = baseline_meta.set_index("ids")["subject"]
    subject_rows: list[pd.DataFrame] = []
    contribution_rows: list[pd.DataFrame] = []
    map_rows: list[dict[str, str]] = []
    for pathway in PATHWAYS:
        values = unstrat[pathway].rename("cpm")
        frame = pd.DataFrame({"subject": sample_subject.loc[values.index], "cpm": values})
        subject = frame.groupby("subject", sort=True).agg(
            baseline_cpm=("cpm", "mean"),
            positive_timepoints=("cpm", lambda x: int((x > 0).sum())),
            positive_sample_mean_cpm=("cpm", lambda x: float(x.loc[x > 0].mean()) if (x > 0).any() else np.nan),
        )
        subject["presence_any"] = (subject["positive_timepoints"] >= 1).astype(int)
        subject["presence_majority"] = (subject["positive_timepoints"] >= 2).astype(int)
        subject["pathway"] = pathway
        subject_rows.append(subject.reset_index())

        rows = strat.loc[
            strat.index.to_series().astype(str).map(
                lambda value: value.startswith(pathway + ":") or value.startswith(pathway + "|")
            ).to_numpy()
        ]
        classes: dict[str, str] = {}
        for row_name in rows.index.astype(str):
            taxon, literal, harmonized, reason = _parse_stratified_taxon(row_name, phoc_epithets)
            if harmonized == "Bacteroides":
                contribution_class = "Bacteroides"
            elif harmonized == "Phocaeicola":
                contribution_class = "Phocaeicola"
            elif harmonized == "unclassified":
                contribution_class = "unclassified_stratum"
            else:
                contribution_class = "other_resolved"
            classes[row_name] = contribution_class
            map_rows.append(
                {
                    "pathway": pathway,
                    "stratified_row": row_name,
                    "taxon_label": taxon,
                    "literal_genus": literal,
                    "harmonized_genus": harmonized,
                    "contribution_class": contribution_class,
                    "harmonization_reason": reason,
                }
            )
        by_class: dict[str, pd.Series] = {}
        for name in ("Bacteroides", "Phocaeicola", "other_resolved", "unclassified_stratum"):
            selected = [row for row, category in classes.items() if category == name]
            by_class[name] = rows.loc[selected].sum(axis=0) if selected else pd.Series(0.0, index=unstrat.index)
        sample_contrib = pd.DataFrame(by_class, index=unstrat.index)
        sample_contrib["stratified_sum"] = rows.sum(axis=0)
        sample_contrib["total_cpm"] = unstrat[pathway]
        remainder = sample_contrib["total_cpm"] - sample_contrib["stratified_sum"]
        if (remainder < -1e-6).any():
            raise ValueError(f"Stratified abundance exceeds total for {pathway}")
        sample_contrib["unresolved_remainder"] = remainder.clip(lower=0.0)
        sample_contrib["subject"] = sample_subject.loc[sample_contrib.index]
        summed = sample_contrib.groupby("subject", sort=True).mean(numeric_only=True)
        summed["pathway"] = pathway
        summed["bact_phoc_cpm"] = summed["Bacteroides"] + summed["Phocaeicola"]
        summed["cpm_minus_bact_phoc"] = (summed["total_cpm"] - summed["bact_phoc_cpm"]).clip(lower=0.0)
        denominator = summed["total_cpm"].replace(0, np.nan)
        for name in (
            "Bacteroides",
            "Phocaeicola",
            "other_resolved",
            "unclassified_stratum",
            "unresolved_remainder",
            "stratified_sum",
            "bact_phoc_cpm",
        ):
            summed[name + "_share"] = summed[name] / denominator
        contribution_rows.append(summed.reset_index())

    subject_pathways = pd.concat(subject_rows, ignore_index=True)
    contributions = pd.concat(contribution_rows, ignore_index=True)
    taxon_map = pd.DataFrame(map_rows)
    if subject_pathways.groupby("pathway")["subject"].nunique().ne(60).any():
        raise ValueError("Incomplete subject-level pathway aggregation")
    return subject_pathways, contributions, taxon_map


def bootstrap_interval(values: Sequence[float], rng: np.random.Generator, n_boot: int = BOOTSTRAPS) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return np.nan, np.nan
    estimates = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        estimates[i] = np.median(rng.choice(array, size=len(array), replace=True))
    return tuple(np.quantile(estimates, [0.025, 0.975]).tolist())


def descriptive_tables(
    baseline_meta: pd.DataFrame,
    unstrat: pd.DataFrame,
    subject_pathways: pd.DataFrame,
    contributions: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_subject = baseline_meta.set_index("ids")["subject"]
    zero_rows: list[dict[str, float | int | str]] = []
    contribution_summary: list[dict[str, float | int | str]] = []
    rng = np.random.default_rng(seed)
    for pathway in PATHWAYS:
        sample_values = unstrat[pathway]
        subj = subject_pathways.loc[subject_pathways["pathway"] == pathway]
        positive_sample = sample_values.loc[sample_values > 0]
        positive_subject = subj.loc[subj["baseline_cpm"] > 0, "baseline_cpm"]
        zero_rows.append(
            {
                "pathway": pathway,
                "pathway_label": PATHWAY_SHORT[pathway],
                "samples": int(len(sample_values)),
                "subjects": int(subj["subject"].nunique()),
                "sample_presence_n": int((sample_values > 0).sum()),
                "sample_presence_rate": float((sample_values > 0).mean()),
                "subject_presence_any_n": int(subj["presence_any"].sum()),
                "subject_presence_any_rate": float(subj["presence_any"].mean()),
                "subject_presence_majority_n": int(subj["presence_majority"].sum()),
                "subject_presence_majority_rate": float(subj["presence_majority"].mean()),
                "positive_sample_median_cpm": float(positive_sample.median()),
                "positive_sample_q1_cpm": float(positive_sample.quantile(0.25)),
                "positive_sample_q3_cpm": float(positive_sample.quantile(0.75)),
                "positive_subject_median_cpm": float(positive_subject.median()),
                "positive_subject_q1_cpm": float(positive_subject.quantile(0.25)),
                "positive_subject_q3_cpm": float(positive_subject.quantile(0.75)),
                "aggregation": "A1-A4 arithmetic mean per subject",
            }
        )
        cont = contributions.loc[contributions["pathway"] == pathway]
        for component, column in (
            ("Bacteroides", "Bacteroides_share"),
            ("Phocaeicola", "Phocaeicola_share"),
            ("Bacteroides_or_Phocaeicola", "bact_phoc_cpm_share"),
            ("other_resolved_taxa", "other_resolved_share"),
            ("unclassified_stratum", "unclassified_stratum_share"),
            ("unresolved_remainder", "unresolved_remainder_share"),
            ("all_taxonomically_stratified", "stratified_sum_share"),
        ):
            values = cont[column].dropna().clip(lower=0.0, upper=1.0)
            low, high = bootstrap_interval(values.to_numpy(), rng)
            contribution_summary.append(
                {
                    "pathway": pathway,
                    "component": component,
                    "n_positive_subjects": int(len(values)),
                    "mean_share": float(values.mean()),
                    "median_share": float(values.median()),
                    "q1_share": float(values.quantile(0.25)),
                    "q3_share": float(values.quantile(0.75)),
                    "median_bootstrap_ci_low": low,
                    "median_bootstrap_ci_high": high,
                    "bootstrap_replicates": BOOTSTRAPS,
                    "bootstrap_seed": seed,
                }
            )
    return pd.DataFrame(zero_rows), pd.DataFrame(contribution_summary)


def gene_level_contribution_audit(
    paths: Paths,
    baseline_meta: pd.DataFrame,
    phoc_epithets: set[str],
    covariates: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Audit Bacteroides/Phocaeicola attribution for nadD, nudC and pncC KOs."""
    unstrat = pd.read_csv(paths.raw / "ko_cpm_unstratified.tsv", sep="\t", index_col=0)
    strat = pd.read_csv(paths.raw / "ko_cpm_stratified.tsv", sep="\t", index_col=0)
    for frame in (unstrat, strat):
        frame.columns = (
            frame.columns.astype(str)
            .str.removesuffix("_Abundance-CPM")
            .str.removesuffix("-CPM")
        )
    sample_subject = baseline_meta.set_index("ids")["subject"]
    if set(sample_subject.index) != set(unstrat.columns):
        raise ValueError("KO table is not the same A1--A4 sample set as the pathway table")
    detail_rows: list[pd.DataFrame] = []
    map_rows: list[dict[str, str]] = []
    for gene, kos in NAD_GENE_KOS.items():
        missing = [ko for ko in kos if ko not in unstrat.index]
        if missing:
            raise ValueError(f"Missing KO rows for {gene}: {missing}")
        total = unstrat.loc[list(kos)].sum(axis=0).astype(float)
        row_mask = strat.index.to_series().astype(str).map(
            lambda value: any(value.startswith(ko + "|") for ko in kos)
        )
        rows = strat.loc[row_mask.to_numpy()].astype(float)
        classes: dict[str, str] = {}
        for row_name in rows.index.astype(str):
            taxon, literal, harmonized, reason = _parse_stratified_taxon(row_name, phoc_epithets)
            if harmonized == "Bacteroides":
                category = "Bacteroides"
            elif harmonized == "Phocaeicola":
                category = "Phocaeicola"
            elif harmonized == "unclassified":
                category = "unclassified_stratum"
            else:
                category = "other_resolved"
            classes[row_name] = category
            map_rows.append(
                {
                    "gene": gene,
                    "ko": row_name.split("|", 1)[0],
                    "stratified_row": row_name,
                    "taxon_label": taxon,
                    "literal_genus": literal,
                    "harmonized_genus": harmonized,
                    "contribution_class": category,
                    "harmonization_reason": reason,
                }
            )
        sample = pd.DataFrame(index=total.index)
        for category in ("Bacteroides", "Phocaeicola", "other_resolved", "unclassified_stratum"):
            selected = [row for row, value in classes.items() if value == category]
            sample[category] = rows.loc[selected].sum(axis=0) if selected else 0.0
        sample["stratified_sum"] = rows.sum(axis=0)
        sample["total_cpm"] = total
        remainder = sample["total_cpm"] - sample["stratified_sum"]
        if remainder.min() < -0.01:
            raise ValueError(f"KO stratified abundance materially exceeds total for {gene}")
        sample["unresolved_remainder"] = remainder.clip(lower=0.0)
        sample["subject"] = sample_subject.loc[sample.index]
        subject = sample.groupby("subject", sort=True).mean(numeric_only=True)
        subject["gene"] = gene
        subject["ko"] = "+".join(kos)
        subject["bact_phoc_cpm"] = subject["Bacteroides"] + subject["Phocaeicola"]
        total_denom = subject["total_cpm"].replace(0, np.nan)
        assigned_denom = subject["stratified_sum"].replace(0, np.nan)
        subject["bacteroides_share_total"] = subject["Bacteroides"] / total_denom
        subject["phocaeicola_share_total"] = subject["Phocaeicola"] / total_denom
        subject["bact_phoc_share_total"] = subject["bact_phoc_cpm"] / total_denom
        subject["bact_phoc_share_assigned"] = subject["bact_phoc_cpm"] / assigned_denom
        subject["resolved_share_total"] = subject["stratified_sum"] / total_denom
        subject["unresolved_remainder_share"] = subject["unresolved_remainder"] / total_denom
        subject = subject.join(covariates[["hdr_ldr_group"]], how="left")
        detail_rows.append(subject.reset_index())
    detail = pd.concat(detail_rows, ignore_index=True)

    rng = np.random.default_rng(seed + 101)
    summary_rows: list[dict[str, object]] = []
    for gene in NAD_GENE_KOS:
        gene_data = detail.loc[detail["gene"] == gene]
        for group in ("All", "LDR", "HDR", "Middle"):
            part = gene_data if group == "All" else gene_data.loc[gene_data["hdr_ldr_group"] == group]
            for measure in (
                "bacteroides_share_total",
                "phocaeicola_share_total",
                "bact_phoc_share_total",
                "bact_phoc_share_assigned",
                "resolved_share_total",
                "unresolved_remainder_share",
            ):
                values = part[measure].dropna().clip(0, 1)
                low, high = bootstrap_interval(values, rng)
                summary_rows.append(
                    {
                        "gene": gene,
                        "ko": "+".join(NAD_GENE_KOS[gene]),
                        "group": group,
                        "measure": measure,
                        "n_subjects": int(len(values)),
                        "mean": float(values.mean()),
                        "median": float(values.median()),
                        "q1": float(values.quantile(0.25)),
                        "q3": float(values.quantile(0.75)),
                        "median_bootstrap_ci_low": low,
                        "median_bootstrap_ci_high": high,
                        "bootstrap_replicates": BOOTSTRAPS,
                        "bootstrap_seed": seed + 101,
                    }
                )
    return detail, pd.DataFrame(summary_rows), pd.DataFrame(map_rows)


MODEL_COVARIATES: Mapping[str, tuple[str, ...]] = {
    "unadjusted": (),
    "bacteroides_adjusted": ("bact_z", "richness_z", "log_depth_z"),
    "phocaeicola_adjusted": ("phoc_z", "richness_z", "log_depth_z"),
    "taxa_separate_adjusted": ("bact_z", "phoc_z", "richness_z", "log_depth_z"),
    "combined_bp_adjusted": ("bp_z", "richness_z", "log_depth_z"),
    "combined_bp_unrarefied_adjusted_sensitivity": (
        "bp_unrarefied_z",
        "richness_unrarefied_z",
        "log_depth_z",
    ),
}


def _base_result(
    *,
    analysis: str,
    recovery_metric: str,
    data_variant: str,
    model: str,
    pathway: str,
    component: str,
    method: str,
    predictor: str,
) -> dict[str, object]:
    return {
        "analysis": analysis,
        "recovery_metric": recovery_metric,
        "data_variant": data_variant,
        "model": model,
        "pathway": pathway,
        "component": component,
        "method": method,
        "predictor": predictor,
        "n": 0,
        "events": np.nan,
        "beta": np.nan,
        "se_hc3": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p_value": np.nan,
        "q_value": np.nan,
        "ratio": np.nan,
        "ratio_ci_low": np.nan,
        "ratio_ci_high": np.nan,
        "effect_scale": "",
        "status": "not_run",
        "note": "",
    }


def fit_ols_hc3(
    data: pd.DataFrame,
    outcome: str,
    predictor: str,
    covariates: Sequence[str],
    base: dict[str, object],
    effect_scale: str,
) -> dict[str, object]:
    columns = [outcome, predictor, *covariates]
    frame = data.loc[:, columns].replace([np.inf, -np.inf], np.nan).dropna()
    result = dict(base)
    result["n"] = int(len(frame))
    result["effect_scale"] = effect_scale
    if len(frame) < len(covariates) + 6 or frame[predictor].nunique() < 2 or frame[outcome].nunique() < 2:
        result["status"] = "non_estimable"
        result["note"] = "insufficient sample size or variation"
        return result
    design = sm.add_constant(frame[[predictor, *covariates]], has_constant="add")
    try:
        fit = sm.OLS(frame[outcome], design).fit(cov_type="HC3")
        beta = float(fit.params[predictor])
        se = float(fit.bse[predictor])
        ci = fit.conf_int(alpha=0.05).loc[predictor].to_numpy(dtype=float)
        p_value = float(fit.pvalues[predictor])
    except (ValueError, np.linalg.LinAlgError) as exc:
        result["status"] = "fit_failed"
        result["note"] = str(exc)
        return result
    result.update(
        beta=beta,
        se_hc3=se,
        ci_low=float(ci[0]),
        ci_high=float(ci[1]),
        p_value=p_value,
        status="ok",
    )
    if effect_scale.startswith("log1p"):
        result.update(
            ratio=float(math.exp(beta)),
            ratio_ci_low=float(math.exp(ci[0])),
            ratio_ci_high=float(math.exp(ci[1])),
        )
    return result


def fit_logit_hc3(
    data: pd.DataFrame,
    outcome: str,
    predictor: str,
    covariates: Sequence[str],
    base: dict[str, object],
) -> dict[str, object]:
    columns = [outcome, predictor, *covariates]
    frame = data.loc[:, columns].replace([np.inf, -np.inf], np.nan).dropna()
    result = dict(base)
    result["n"] = int(len(frame))
    result["events"] = int(frame[outcome].sum()) if len(frame) else 0
    result["effect_scale"] = "log odds of pathway detection"
    if frame[outcome].nunique() < 2:
        result["status"] = "non_estimable_no_presence_variation"
        result["note"] = "all subjects have the same detection status"
        return result
    if len(frame) < len(covariates) + 8 or frame[predictor].nunique() < 2:
        result["status"] = "non_estimable"
        result["note"] = "insufficient sample size or predictor variation"
        return result
    design = sm.add_constant(frame[[predictor, *covariates]], has_constant="add")
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", PerfectSeparationWarning)
            fit = sm.GLM(frame[outcome], design, family=sm.families.Binomial()).fit(
                cov_type="HC3", maxiter=200
            )
        beta = float(fit.params[predictor])
        se = float(fit.bse[predictor])
        ci = fit.conf_int(alpha=0.05).loc[predictor].to_numpy(dtype=float)
        p_value = float(fit.pvalues[predictor])
        unstable = any(issubclass(item.category, PerfectSeparationWarning) for item in caught)
        unstable = unstable or not np.isfinite([beta, se, *ci, p_value]).all() or abs(beta) > 20 or se > 20
    except (ValueError, np.linalg.LinAlgError) as exc:
        result["status"] = "fit_failed"
        result["note"] = str(exc)
        return result
    result.update(
        beta=beta,
        se_hc3=se,
        ci_low=float(ci[0]),
        ci_high=float(ci[1]),
        # Wald p values and odds-ratio intervals are not defensible under
        # complete/quasi-complete separation. Preserve the raw coefficient as
        # a diagnostic, but exclude the inferential fields from BH families.
        p_value=np.nan if unstable else p_value,
        ratio=np.nan if unstable else float(math.exp(np.clip(beta, -700, 700))),
        ratio_ci_low=np.nan if unstable else float(math.exp(np.clip(ci[0], -700, 700))),
        ratio_ci_high=np.nan if unstable else float(math.exp(np.clip(ci[1], -700, 700))),
        status="unstable_separation" if unstable else "ok",
        note=(
            "possible complete/quasi-complete separation; Wald p value and odds-ratio CI suppressed"
            if unstable
            else ""
        ),
    )
    return result


def run_hurdle_models(subject_data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    variants = {
        "original_total": "baseline_cpm",
        "minus_stratified_bact_phoc": "cpm_minus_bact_phoc",
    }
    analyses = [
        ("continuous_recovery", "bray", "bray_recovery_z", subject_data.index == subject_data.index),
        ("continuous_recovery", "jaccard_sensitivity", "jaccard_recovery_z", subject_data.index == subject_data.index),
        ("continuous_recovery", "bray_auc_sensitivity", "bray_auc_recovery_z", subject_data.index == subject_data.index),
        ("continuous_recovery", "jaccard_auc_sensitivity", "jaccard_auc_recovery_z", subject_data.index == subject_data.index),
        (
            "endpoint_defined_hdr_ldr_descriptive",
            "not_applicable",
            "hdr_indicator",
            subject_data["hdr_ldr_group"].isin(["HDR", "LDR"]),
        ),
    ]
    for analysis, recovery_metric, predictor, mask in analyses:
        data_analysis = subject_data.loc[mask].copy()
        if predictor == "hdr_indicator":
            data_analysis[predictor] = (data_analysis["hdr_ldr_group"] == "HDR").astype(int)
        for pathway in PATHWAYS:
            path_data = data_analysis.loc[data_analysis["pathway"] == pathway].copy()
            for variant, abundance_col in variants.items():
                path_data["presence_model"] = (path_data[abundance_col] > 0).astype(int)
                path_data["log1p_positive"] = np.log1p(path_data[abundance_col].where(path_data[abundance_col] > 0))
                path_data["rank_positive"] = rank_normal(path_data[abundance_col].where(path_data[abundance_col] > 0))
                for model, covariates in MODEL_COVARIATES.items():
                    logistic_base = _base_result(
                        analysis=analysis,
                        recovery_metric=recovery_metric,
                        data_variant=variant,
                        model=model,
                        pathway=pathway,
                        component="presence",
                        method="logistic_glm_hc3",
                        predictor=predictor,
                    )
                    rows.append(
                        fit_logit_hc3(
                            path_data,
                            "presence_model",
                            predictor,
                            covariates,
                            logistic_base,
                        )
                    )
                    if variant == "original_total":
                        majority_base = _base_result(
                            analysis=analysis,
                            recovery_metric=recovery_metric,
                            data_variant="original_total_majority_detection",
                            model=model,
                            pathway=pathway,
                            component="presence",
                            method="logistic_glm_hc3_majority_sensitivity",
                            predictor=predictor,
                        )
                        rows.append(
                            fit_logit_hc3(
                                path_data,
                                "presence_majority",
                                predictor,
                                covariates,
                                majority_base,
                            )
                        )
                    positive = path_data.loc[path_data[abundance_col] > 0].copy()
                    abundance_base = _base_result(
                        analysis=analysis,
                        recovery_metric=recovery_metric,
                        data_variant=variant,
                        model=model,
                        pathway=pathway,
                        component="positive_abundance",
                        method="ols_log1p_hc3",
                        predictor=predictor,
                    )
                    rows.append(
                        fit_ols_hc3(
                            positive,
                            "log1p_positive",
                            predictor,
                            covariates,
                            abundance_base,
                            "log1p(CPM) difference per predictor unit",
                        )
                    )
                    rank_base = _base_result(
                        analysis=analysis,
                        recovery_metric=recovery_metric,
                        data_variant=variant,
                        model=model,
                        pathway=pathway,
                        component="positive_abundance",
                        method="ols_rank_normal_hc3_sensitivity",
                        predictor=predictor,
                    )
                    rows.append(
                        fit_ols_hc3(
                            positive,
                            "rank_positive",
                            predictor,
                            covariates,
                            rank_base,
                            "rank-normal SD difference per predictor unit",
                        )
                    )
    results = pd.DataFrame(rows)
    # Multiplicity families combine estimable presence and positive-abundance
    # parts across the three pathways. Rank-normal models are separate sensitivity families.
    family_branch = np.select(
        [
            results["method"].eq("ols_rank_normal_hc3_sensitivity"),
            results["method"].eq("logistic_glm_hc3_majority_sensitivity"),
        ],
        ["rank_sensitivity", "majority_detection_sensitivity"],
        default="primary_hurdle",
    )
    results["family_id"] = (
        results["analysis"].astype(str)
        + "|"
        + results["recovery_metric"].astype(str)
        + "|"
        + results["data_variant"].astype(str)
        + "|"
        + results["model"].astype(str)
        + "|"
        + family_branch
    )
    results["q_value"] = results.groupby("family_id", sort=False)["p_value"].transform(bh_adjust)
    results["selection_coupled"] = results["analysis"].eq("endpoint_defined_hdr_ldr_descriptive")
    results["analysis_tier"] = "sensitivity"
    results.loc[results["selection_coupled"], "analysis_tier"] = "descriptive_selection_coupled"
    primary = (
        results["analysis"].eq("continuous_recovery")
        & results["recovery_metric"].eq("bray")
        & results["data_variant"].eq("original_total")
        & results["model"].eq("combined_bp_adjusted")
        & results["method"].isin(["logistic_glm_hc3", "ols_log1p_hc3"])
    )
    results.loc[primary, "analysis_tier"] = "primary"
    return results


def build_subject_table(
    covariates: pd.DataFrame,
    subject_pathways: pd.DataFrame,
    contributions: pd.DataFrame,
) -> pd.DataFrame:
    merged = subject_pathways.merge(contributions, on=["subject", "pathway"], validate="one_to_one")
    merged = merged.merge(covariates.reset_index(), on="subject", validate="many_to_one")
    if merged.shape[0] != 180:
        raise ValueError("Expected 60 subjects x 3 pathways")
    return merged.sort_values(["pathway", "subject"], kind="mergesort").reset_index(drop=True)


def run_baseline_aggregation_sensitivity(
    baseline_meta: pd.DataFrame,
    unstrat: pd.DataFrame,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    """Check that A1--A4 mean aggregation is not a single-summary artefact."""
    meta_lookup = baseline_meta.set_index("ids")[["subject", "canon"]]
    rows: list[dict[str, object]] = []
    schemes = {
        "mean_A1_A4_primary": (set(BASELINE), "mean"),
        "median_A1_A4": (set(BASELINE), "median"),
        "mean_A2_A4": ({"A2", "A3", "A4"}, "mean"),
        "A4_only": ({"A4"}, "mean"),
    }
    for pathway in PATHWAYS:
        long = pd.DataFrame(
            {
                "cpm": unstrat[pathway],
                "subject": meta_lookup.loc[unstrat.index, "subject"],
                "canon": meta_lookup.loc[unstrat.index, "canon"],
            }
        )
        for scheme, (visits, reducer) in schemes.items():
            selected = long.loc[long["canon"].isin(visits)]
            if reducer == "median":
                aggregated = selected.groupby("subject")["cpm"].median()
            else:
                aggregated = selected.groupby("subject")["cpm"].mean()
            frame = covariates.join(aggregated.rename("baseline_cpm"), how="inner").copy()
            frame["presence"] = (frame["baseline_cpm"] > 0).astype(int)
            frame["log1p_positive"] = np.log1p(frame["baseline_cpm"].where(frame["baseline_cpm"] > 0))
            frame["rank_positive"] = rank_normal(frame["baseline_cpm"].where(frame["baseline_cpm"] > 0))
            predictor = "bray_recovery_z"
            covar = MODEL_COVARIATES["combined_bp_adjusted"]
            base = _base_result(
                analysis="baseline_aggregation_sensitivity",
                recovery_metric="bray",
                data_variant=scheme,
                model="combined_bp_adjusted",
                pathway=pathway,
                component="presence",
                method="logistic_glm_hc3",
                predictor=predictor,
            )
            rows.append(fit_logit_hc3(frame, "presence", predictor, covar, base))
            positive = frame.loc[frame["baseline_cpm"] > 0].copy()
            base = _base_result(
                analysis="baseline_aggregation_sensitivity",
                recovery_metric="bray",
                data_variant=scheme,
                model="combined_bp_adjusted",
                pathway=pathway,
                component="positive_abundance",
                method="ols_log1p_hc3",
                predictor=predictor,
            )
            rows.append(
                fit_ols_hc3(
                    positive,
                    "log1p_positive",
                    predictor,
                    covar,
                    base,
                    "log1p(CPM) difference per 1 SD better Bray recovery",
                )
            )
            base = _base_result(
                analysis="baseline_aggregation_sensitivity",
                recovery_metric="bray",
                data_variant=scheme,
                model="combined_bp_adjusted",
                pathway=pathway,
                component="positive_abundance",
                method="ols_rank_normal_hc3_sensitivity",
                predictor=predictor,
            )
            rows.append(
                fit_ols_hc3(
                    positive,
                    "rank_positive",
                    predictor,
                    covar,
                    base,
                    "rank-normal SD difference per 1 SD better Bray recovery",
                )
            )
    results = pd.DataFrame(rows)
    results["family_id"] = (
        "baseline_aggregation_sensitivity|bray|"
        + results["data_variant"].astype(str)
        + "|combined_bp_adjusted|"
        + np.where(results["method"].eq("ols_rank_normal_hc3_sensitivity"), "rank_sensitivity", "hurdle")
    )
    results["q_value"] = results.groupby("family_id", sort=False)["p_value"].transform(bh_adjust)
    results["selection_coupled"] = False
    results["analysis_tier"] = "sensitivity"
    return results


def classify_attribution(models: pd.DataFrame, contributions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    continuous = models.loc[
        (models["analysis"] == "continuous_recovery")
        & (models["recovery_metric"] == "bray")
        & (models["method"] == "ols_log1p_hc3")
        & (models["component"] == "positive_abundance")
        & (models["data_variant"] == "original_total")
    ]
    selected = models.loc[
        (models["analysis"] == "endpoint_defined_hdr_ldr_descriptive")
        & (models["recovery_metric"] == "not_applicable")
        & (models["method"] == "ols_log1p_hc3")
        & (models["component"] == "positive_abundance")
        & (models["data_variant"] == "original_total")
    ]
    for pathway in PATHWAYS:
        cont = contributions.loc[(contributions["pathway"] == pathway) & (contributions["total_cpm"] > 0)]
        share = cont["bact_phoc_cpm_share"].clip(0, 1)
        unresolved = cont["unresolved_remainder_share"].clip(0, 1)
        c_primary = continuous.loc[
            (continuous["pathway"] == pathway) & (continuous["model"] == "combined_bp_adjusted")
        ]
        group_unadjusted = selected.loc[(selected["pathway"] == pathway) & (selected["model"] == "unadjusted")]
        group_adjusted = selected.loc[
            (selected["pathway"] == pathway) & (selected["model"] == "combined_bp_adjusted")
        ]
        continuous_beta = float(c_primary["beta"].iloc[0]) if len(c_primary) else np.nan
        continuous_q = float(c_primary["q_value"].iloc[0]) if len(c_primary) else np.nan
        group_beta_unadjusted = float(group_unadjusted["beta"].iloc[0]) if len(group_unadjusted) else np.nan
        group_q_unadjusted = float(group_unadjusted["q_value"].iloc[0]) if len(group_unadjusted) else np.nan
        group_beta_adjusted = float(group_adjusted["beta"].iloc[0]) if len(group_adjusted) else np.nan
        group_q_adjusted = float(group_adjusted["q_value"].iloc[0]) if len(group_adjusted) else np.nan
        attenuation = np.nan
        if np.isfinite(group_beta_unadjusted) and abs(group_beta_unadjusted) > 1e-12 and np.isfinite(group_beta_adjusted):
            attenuation = 1.0 - abs(group_beta_adjusted) / abs(group_beta_unadjusted)
        median_share = float(share.median()) if len(share) else np.nan
        median_unresolved = float(unresolved.median()) if len(unresolved) else np.nan
        if np.isfinite(continuous_q) and continuous_q >= 0.05:
            continuous_judgment = "no_supported_primary_continuous_association_to_attribute"
        else:
            continuous_judgment = "association_present_but_taxonomic_driver_inconclusive"
        if (
            np.isfinite(group_q_unadjusted)
            and group_q_unadjusted < 0.05
            and np.isfinite(group_q_adjusted)
            and group_q_adjusted >= 0.05
            and np.isfinite(attenuation)
            and attenuation >= 0.50
        ):
            selected_judgment = "selection_coupled_contrast_attenuated_after_joint_adjustment"
        else:
            selected_judgment = "selection_coupled_driver_inconclusive"
        rows.append(
            {
                "pathway": pathway,
                "median_bact_phoc_share": median_share,
                "median_unresolved_remainder_share": median_unresolved,
                "primary_continuous_beta_adjusted": continuous_beta,
                "primary_continuous_q": continuous_q,
                "primary_continuous_judgment": continuous_judgment,
                "selection_coupled_beta_unadjusted": group_beta_unadjusted,
                "selection_coupled_q_unadjusted": group_q_unadjusted,
                "selection_coupled_beta_joint_adjusted": group_beta_adjusted,
                "selection_coupled_q_joint_adjusted": group_q_adjusted,
                "selection_coupled_absolute_effect_attenuation_fraction": attenuation,
                "selection_coupled_judgment": selected_judgment,
                "caveat": (
                    "Joint adjustment cannot isolate Bacteroides/Phocaeicola from richness/depth; "
                    "taxonomically unresolved CPM prevents a definitive driver assignment. "
                    "This is an attribution audit, not evidence about activity, flux, or causality."
                ),
            }
        )
    return pd.DataFrame(rows)


def _forest_rows(models: pd.DataFrame, recovery_metric: str, method: str) -> pd.DataFrame:
    return models.loc[
        (models["analysis"] == "continuous_recovery")
        & (models["recovery_metric"] == recovery_metric)
        & (models["model"] == "combined_bp_adjusted")
        & (models["method"] == method)
        & (models["component"] == "positive_abundance")
        & (models["status"] == "ok")
    ].copy()


def set_figure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.4,
            "ytick.labelsize": 6.4,
            "legend.fontsize": 6.2,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.13, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def render_summary_figure(
    subject_data: pd.DataFrame,
    zero_summary: pd.DataFrame,
    contribution_summary: pd.DataFrame,
    models: pd.DataFrame,
    paths: Paths,
    seed: int,
) -> dict[str, object]:
    """Render a quantitative-grid summary using Python/matplotlib only."""
    set_figure_style()
    rng = np.random.default_rng(seed)
    fig = plt.figure(figsize=(7.2, 6.8), constrained_layout=False)
    grid = fig.add_gridspec(2, 2, left=0.13, right=0.985, bottom=0.18, top=0.93, wspace=0.38, hspace=0.42)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # a: subject-level baseline distributions; the endpoint-defined groups are
    # deliberately shown as a descriptive selection-coupled contrast.
    x_positions = np.arange(len(PATHWAYS), dtype=float)
    for i, pathway in enumerate(PATHWAYS):
        part = subject_data.loc[
            (subject_data["pathway"] == pathway) & subject_data["hdr_ldr_group"].isin(["LDR", "HDR"])
        ]
        for group, offset in (("LDR", -0.13), ("HDR", 0.13)):
            values = np.log10(part.loc[part["hdr_ldr_group"] == group, "baseline_cpm"] + 1.0).to_numpy()
            jitter = rng.uniform(-0.045, 0.045, len(values))
            ax_a.scatter(
                np.full(len(values), i + offset) + jitter,
                values,
                s=15,
                color=GROUP_COLORS[group],
                alpha=0.78,
                linewidth=0.25,
                edgecolor="white",
                zorder=3,
            )
            if len(values):
                median = np.median(values)
                ax_a.plot([i + offset - 0.07, i + offset + 0.07], [median, median], color="#222222", lw=1.2, zorder=4)
    ax_a.set_xticks(x_positions, [PATHWAY_SHORT[p] for p in PATHWAYS], rotation=18, ha="right")
    ax_a.set_ylabel(r"Baseline inferred potential, log$_{10}$(CPM + 1)")
    ax_a.set_title("Endpoint-defined groups (descriptive only)", loc="left", fontweight="bold")
    ax_a.text(
        0.01,
        0.98,
        "A1–A4 mean; n=15/group\nselection-coupled HDR/LDR",
        transform=ax_a.transAxes,
        va="top",
        ha="left",
        fontsize=6.1,
        color="#555555",
    )
    ax_a.legend(
        handles=[
            Line2D([0], [0], marker="o", lw=0, label="LDR", color=GROUP_COLORS["LDR"], markersize=4),
            Line2D([0], [0], marker="o", lw=0, label="HDR", color=GROUP_COLORS["HDR"], markersize=4),
        ],
        frameon=False,
        loc="upper right",
        ncol=2,
        handletextpad=0.35,
        columnspacing=0.8,
    )
    panel_label(ax_a, "a")

    # b: prevalence and positive-value summary.  This visibly shows why the
    # presence component is non-estimable for the two ubiquitous pathways.
    y = np.arange(len(PATHWAYS))[::-1]
    sample_rate = zero_summary.set_index("pathway").loc[list(PATHWAYS), "sample_presence_rate"].to_numpy()
    subject_rate = zero_summary.set_index("pathway").loc[list(PATHWAYS), "subject_presence_any_rate"].to_numpy()
    ax_b.hlines(y, sample_rate, subject_rate, color="#C9CFD6", lw=2.2, zorder=1)
    ax_b.scatter(sample_rate, y, s=26, facecolor="white", edgecolor="#68737D", lw=0.9, label="Samples (n=240)", zorder=3)
    ax_b.scatter(subject_rate, y, s=27, color="#2D6A8A", edgecolor="white", lw=0.35, label="Subjects, any A1–A4 (n=60)", zorder=4)
    ax_b.set_yticks(y, [PATHWAY_SHORT[p] for p in PATHWAYS])
    ax_b.set_xlim(-0.02, 1.04)
    ax_b.set_xlabel("Detection proportion")
    ax_b.set_title("Zero inflation is pathway-specific", loc="left", fontweight="bold")
    ax_b.axvline(1.0, color="#DDDDDD", lw=0.7, ls="--", zorder=0)
    ax_b.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.01, 0.97))
    panel_label(ax_b, "b")

    # c: hero forest plot, continuous within-subject recovery, original versus
    # taxonomically depleted pathway estimates.
    forest = _forest_rows(models, "bray", "ols_log1p_hc3")
    labels = list(PATHWAYS)
    y_base = np.arange(len(labels))[::-1]
    for variant, offset, marker, color, label in (
        ("original_total", 0.10, "o", "#3F6C8E", "Total pathway CPM"),
        ("minus_stratified_bact_phoc", -0.10, "s", "#B0734F", "After subtracting Bact./Phoc. strata"),
    ):
        part = forest.loc[forest["data_variant"] == variant].set_index("pathway")
        for idx, pathway in enumerate(labels):
            if pathway not in part.index:
                continue
            row = part.loc[pathway]
            ax_c.errorbar(
                float(row["beta"]),
                y_base[idx] + offset,
                xerr=np.array([[float(row["beta"] - row["ci_low"])], [float(row["ci_high"] - row["beta"])]]),
                fmt=marker,
                ms=4.2,
                color=color,
                mec="white",
                mew=0.35,
                capsize=2,
                lw=1.0,
                label=label if idx == 0 else None,
            )
    ax_c.axvline(0, color="#555555", lw=0.75, ls="--")
    ax_c.set_yticks(y_base, [PATHWAY_SHORT[p] for p in labels])
    ax_c.set_xlabel("Adjusted β in log1p(CPM) per 1 SD better Bray recovery\n(95% CI; HC3)")
    ax_c.set_title("Full-cohort continuous recovery (n=60; rarefied)", loc="left", fontweight="bold")
    c_handles, c_labels = ax_c.get_legend_handles_labels()
    fig.legend(
        c_handles,
        c_labels,
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.13, 0.052),
        ncol=1,
        fontsize=5.8,
    )
    panel_label(ax_c, "c")

    # d: attribution audit. Use median shares, retaining the unresolved portion
    # rather than renormalising only the identifiable strata.
    summary = contribution_summary.pivot(index="pathway", columns="component", values="mean_share")
    order = list(PATHWAYS)
    bottoms = np.zeros(len(order), dtype=float)
    components = [
        ("Bacteroides", "#6A8FB3", "Bacteroides"),
        ("Phocaeicola", "#B7CCE0", "Phocaeicola\n(harmonized)"),
        ("other_resolved_taxa", "#9ABF9A", "Other resolved"),
        ("unclassified_stratum", "#D4B27A", "Unclassified stratum"),
        ("unresolved_remainder", "#E4E4E4", "Unresolved remainder"),
    ]
    for component, color, label in components:
        values = summary.reindex(order)[component].fillna(0).to_numpy(dtype=float)
        ax_d.barh(y_base, values, left=bottoms, height=0.58, color=color, edgecolor="white", lw=0.4, label=label)
        bottoms += values
    ax_d.set_yticks(y_base, [PATHWAY_SHORT[p] for p in order])
    ax_d.set_xlim(0, 1)
    ax_d.set_xlabel("Mean share of total baseline pathway CPM")
    ax_d.set_title("Taxonomic attribution of inferred potential", loc="left", fontweight="bold")
    handles, legend_labels = ax_d.get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.70, 0.055),
        ncol=3,
        fontsize=5.5,
        columnspacing=0.7,
    )
    panel_label(ax_d, "d")

    fig.suptitle(
        "Subject-level audit of baseline inferred microbial functional potential",
        x=0.13,
        y=0.985,
        ha="left",
        fontsize=10.5,
        fontweight="bold",
    )
    fig.text(
        0.985,
        0.014,
        "CPM reflects inferred gene-content potential—not NAD concentration, pathway activity, flux, or mechanism.",
        ha="right",
        va="bottom",
        fontsize=6.2,
        color="#555555",
    )

    paths.figures.mkdir(parents=True, exist_ok=True)
    pdf = paths.figures / "human_pathway_reanalysis_summary.pdf"
    tif = paths.figures / "human_pathway_reanalysis_summary.tif"
    fig.savefig(pdf, bbox_inches=None, facecolor="white")
    fig.savefig(tif, dpi=600, bbox_inches=None, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    with Image.open(tif) as opened:
        if opened.mode != "RGB":
            rgb = opened.convert("RGB")
            rgb.save(tif, dpi=(600, 600), compression="tiff_lzw")
    return {
        "pdf": str(pdf.relative_to(paths.project_root)),
        "tif": str(tif.relative_to(paths.project_root)),
        "width_in": 7.2,
        "height_in": 6.8,
        "tiff_dpi": 600,
        "backend": "Python/matplotlib",
        "archetype": "quantitative grid with a continuous-recovery forest hero panel",
    }


def build_gene_attribution_figure_source(
    gene_contributions: pd.DataFrame,
    seed: int,
) -> pd.DataFrame:
    """Create a plot-only descriptive summary from the audited subject table.

    This does not fit a new inferential model.  Each component is divided by
    the participant's total unstratified KO CPM before the subject-level median
    and its non-parametric bootstrap interval are calculated.
    """
    rng = np.random.default_rng(seed + 202)
    specs = (
        ("bacteroides", "Bacteroides", lambda frame: frame["Bacteroides"]),
        ("phocaeicola", "Phocaeicola", lambda frame: frame["Phocaeicola"]),
        ("other_resolved", "Other resolved taxa", lambda frame: frame["other_resolved"]),
        (
            "unclassified_unresolved",
            "Unclassified / unresolved",
            lambda frame: frame["unclassified_stratum"] + frame["unresolved_remainder"],
        ),
        (
            "bact_phoc_combined",
            "Bacteroides + Phocaeicola",
            lambda frame: frame["Bacteroides"] + frame["Phocaeicola"],
        ),
    )
    rows: list[dict[str, object]] = []
    for gene, kos in NAD_GENE_KOS.items():
        part = gene_contributions.loc[gene_contributions["gene"] == gene].copy()
        denominator = part["total_cpm"].where(part["total_cpm"] > 0)
        for component, label, numerator in specs:
            values = (numerator(part) / denominator).dropna().clip(0, 1)
            low, high = bootstrap_interval(values.to_numpy(), rng)
            rows.append(
                {
                    "gene": gene,
                    "ko": "+".join(kos),
                    "group": "All",
                    "component": component,
                    "component_label": label,
                    "n_subjects": int(len(values)),
                    "median_share": float(values.median()),
                    "q1_share": float(values.quantile(0.25)),
                    "q3_share": float(values.quantile(0.75)),
                    "median_bootstrap_ci_low": low,
                    "median_bootstrap_ci_high": high,
                    "bootstrap_replicates": BOOTSTRAPS,
                    "bootstrap_seed": seed + 202,
                    "denominator": "subject-level total unstratified KO CPM",
                    "baseline_aggregation": "A1-A4 arithmetic mean per subject",
                    "source_table": "nad_gene_stratified_contributions_by_subject.csv",
                }
            )
    return pd.DataFrame(rows)


def render_gene_attribution_figure(
    gene_contributions: pd.DataFrame,
    figure_source: pd.DataFrame,
    paths: Paths,
    seed: int,
) -> dict[str, object]:
    """Render a compact reviewer-facing KO attribution figure in Python."""
    set_figure_style()
    genes = list(NAD_GENE_KOS)
    centers = dict(zip(genes, np.arange(len(genes), dtype=float)[::-1]))
    gene_labels = {gene: f"{gene}\n{'+'.join(NAD_GENE_KOS[gene])}" for gene in genes}
    component_style = {
        "bacteroides": ("#365F7D", "o"),
        "phocaeicola": ("#79A9C5", "s"),
        "other_resolved": ("#6F9C73", "^"),
        "unclassified_unresolved": ("#A79C8E", "D"),
    }
    offsets = {
        "bacteroides": 0.24,
        "phocaeicola": 0.08,
        "other_resolved": -0.08,
        "unclassified_unresolved": -0.24,
    }

    fig = plt.figure(figsize=(7.2, 3.8), constrained_layout=False)
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.68, 1.0),
        left=0.12,
        right=0.985,
        bottom=0.25,
        top=0.72,
        wspace=0.24,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])

    # a: separately estimated component medians.  They are not stacked because
    # component-wise medians need not sum to exactly 100%.
    component_rows = figure_source.loc[
        figure_source["component"].isin(component_style)
    ].copy()
    for component, (color, marker) in component_style.items():
        part = component_rows.loc[component_rows["component"] == component].set_index("gene")
        for gene in genes:
            row = part.loc[gene]
            median = 100 * float(row["median_share"])
            low = 100 * float(row["median_bootstrap_ci_low"])
            high = 100 * float(row["median_bootstrap_ci_high"])
            ypos = centers[gene] + offsets[component]
            ax_a.errorbar(
                median,
                ypos,
                xerr=np.array([[median - low], [high - median]]),
                fmt=marker,
                color=color,
                mec="white",
                mew=0.45,
                ms=5.2,
                capsize=2.0,
                elinewidth=1.1,
                zorder=3,
            )
            ax_a.text(
                min(high + 1.2, 77.0),
                ypos,
                f"{median:.1f}",
                va="center",
                ha="left",
                fontsize=5.8,
                color="#333333",
            )
    for separator in (0.5, 1.5):
        ax_a.axhline(separator, color="#ECEFF1", lw=0.7, zorder=0)
    ax_a.set_yticks([centers[gene] for gene in genes], [gene_labels[gene] for gene in genes])
    ax_a.set_xlim(-1, 80)
    ax_a.set_xticks([0, 20, 40, 60, 80])
    ax_a.set_ylim(-0.42, 2.42)
    ax_a.set_xlabel("Share of total unstratified KO CPM (%)")
    ax_a.set_title("Source-resolved attribution", loc="left", fontweight="bold")
    ax_a.grid(axis="x", color="#E7EAED", lw=0.65)
    ax_a.set_axisbelow(True)
    ax_a.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker=component_style[key][1],
                color=component_style[key][0],
                markeredgecolor="white",
                markeredgewidth=0.4,
                lw=1.0,
                label=label,
            )
            for key, label in (
                ("bacteroides", "Bacteroides"),
                ("phocaeicola", "Phocaeicola"),
                ("other_resolved", "Other resolved"),
                ("unclassified_unresolved", "Unclassified / unresolved"),
            )
        ],
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(-0.02, 1.08),
        ncol=2,
        columnspacing=0.9,
        handlelength=1.4,
        handletextpad=0.35,
        borderaxespad=0,
        fontsize=5.9,
    )
    panel_label(ax_a, "a")

    # b: retain the individual-subject distribution for the combined genus
    # contribution, with the same subject-level denominator as panel a.
    rng = np.random.default_rng(seed + 303)
    combined = figure_source.loc[
        figure_source["component"] == "bact_phoc_combined"
    ].set_index("gene")
    for gene in genes:
        part = gene_contributions.loc[gene_contributions["gene"] == gene]
        values = 100 * (part["Bacteroides"] + part["Phocaeicola"]) / part["total_cpm"]
        values = values.replace([np.inf, -np.inf], np.nan).dropna().clip(0, 100)
        ypos = centers[gene]
        jitter = rng.uniform(-0.16, 0.16, len(values))
        ax_b.scatter(
            values,
            np.full(len(values), ypos) + jitter,
            s=9,
            color="#AAB4BC",
            alpha=0.48,
            linewidth=0,
            zorder=1,
        )
        row = combined.loc[gene]
        median = 100 * float(row["median_share"])
        low = 100 * float(row["median_bootstrap_ci_low"])
        high = 100 * float(row["median_bootstrap_ci_high"])
        ax_b.errorbar(
            median,
            ypos,
            xerr=np.array([[median - low], [high - median]]),
            fmt="D",
            color="#8E4F3D",
            mec="white",
            mew=0.5,
            ms=5.5,
            capsize=2.4,
            elinewidth=1.35,
            zorder=4,
        )
        ax_b.text(
            min(high + 2.2, 92),
            ypos,
            f"{median:.1f}%",
            va="center",
            ha="left",
            fontsize=6.2,
            fontweight="bold",
            color="#704033",
        )
    ax_b.axvline(50, color="#D9DDE0", lw=0.7, ls="--", zorder=0)
    ax_b.set_yticks([centers[gene] for gene in genes], [gene_labels[gene] for gene in genes])
    ax_b.set_xlim(-2, 103)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_ylim(-0.42, 2.42)
    ax_b.set_xlabel("Combined share of total KO CPM (%)")
    ax_b.set_title("Bacteroides + Phocaeicola", loc="left", fontweight="bold")
    ax_b.grid(axis="x", color="#EEF0F2", lw=0.6)
    ax_b.set_axisbelow(True)
    panel_label(ax_b, "b")

    fig.suptitle(
        "Taxonomic attribution of NAD-related inferred functional potential",
        x=0.12,
        y=0.975,
        ha="left",
        fontsize=10.3,
        fontweight="bold",
    )
    fig.text(
        0.12,
        0.895,
        "All participants; A1–A4 arithmetic mean per subject (n=60 per gene)",
        ha="left",
        va="center",
        fontsize=6.8,
        color="#4E5961",
    )
    fig.text(
        0.12,
        0.035,
        "Median and 95% bootstrap CI (5,000 resamples); grey points show participants. "
        "Denominator: each participant's total unstratified KO CPM.\n"
        "Unclassified / unresolved = unclassified stratum + total-minus-stratified remainder. "
        "Component medians need not sum to 100%.\n"
        "CPM denotes inferred gene-content potential—not activity, metabolite concentration, or flux.",
        ha="left",
        va="bottom",
        fontsize=5.5,
        color="#515960",
        linespacing=1.35,
    )

    paths.figures.mkdir(parents=True, exist_ok=True)
    base = paths.figures / "nad_gene_taxonomic_attribution"
    pdf = base.with_suffix(".pdf")
    svg = base.with_suffix(".svg")
    tif = base.with_suffix(".tif")
    fig.savefig(pdf, bbox_inches=None, facecolor="white")
    fig.savefig(svg, bbox_inches=None, facecolor="white")
    fig.savefig(tif, dpi=600, bbox_inches=None, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    with Image.open(tif) as opened:
        if opened.mode != "RGB":
            rgb = opened.convert("RGB")
            rgb.save(tif, dpi=(600, 600), compression="tiff_lzw")
    return {
        "pdf": str(pdf.relative_to(paths.project_root)),
        "svg": str(svg.relative_to(paths.project_root)),
        "tif": str(tif.relative_to(paths.project_root)),
        "source_table": str(
            (paths.tables / "nad_gene_attribution_figure_source.csv").relative_to(paths.project_root)
        ),
        "width_in": 7.2,
        "height_in": 3.8,
        "tiff_dpi": 600,
        "backend": "Python/matplotlib",
        "archetype": "compact quantitative grid with component intervals and subject distribution",
        "group": "All participants",
        "n_subjects_per_gene": 60,
        "center": "subject-level median",
        "interval": "95% non-parametric bootstrap CI; 5,000 resamples",
        "denominator": "subject-level total unstratified KO CPM",
    }


def write_readme(paths: Paths) -> None:
    text = """# Human inferred-pathway reviewer reanalysis

This directory implements the R1.8 / R2.1--R2.4 human functional-pathway audit.
The input pathway table contains only A1--A4, so the sampling unit is the
participant (60 ciprofloxacin-exposed participants), never the 240 repeated
baseline samples.

## Estimand and model contract

- **Primary association:** baseline inferred functional potential versus the
  continuous C4--C6 Bray--Curtis return to each participant's A1--A4 taxonomic
  centroid in the Methods-faithful rarefied Bracken branch. A larger recovery
  score means closer return to baseline.
- **Hurdle strategy:** pathway detection is analysed with logistic GLM and
  HC3 covariance; abundance among pathway-positive participants is analysed
  as `log1p(CPM)` with OLS-HC3. A rank-normal OLS-HC3 model is the tail-robust
  sensitivity analysis. These are association models with pathway status/
  abundance as the response; they do not establish prospective prediction or
  causal direction. A component is explicitly non-estimable if all participants
  have identical detection status.
- **Single primary specification:** rarefied endpoint Bray recovery, original
  total pathway CPM, and combined Bacteroides/Phocaeicola + richness + raw
  bacterial-read-depth adjustment. All alternative adjustment sets, recovery
  metrics, taxonomic subtraction and rank models are labelled sensitivity
  analyses; they are not used to select a favourable result.
- **Adjustment:** baseline Bacteroides and Phocaeicola are tested separately,
  jointly, and as a combined abundance, always alongside species richness and
  log bacterial Bracken read depth. Continuous covariates are standardized.
- **Multiplicity:** BH correction is applied within each declared model family
  across the three pathways and the estimable hurdle components. Rank models
  and sensitivity outcomes form separate families.
- **Sensitivity:** Jaccard return, actual-day-normalized Bray and Jaccard AUC,
  majority detection (>=2/4 baseline visits), rank-normal positive abundance,
  mean/median/A2--A4/A4-only baseline aggregation, and subtraction of
  stratified Bacteroides/Phocaeicola contributions.
- **HDR/LDR:** bottom/top 15 endpoint-diversity groups are reported only as
  selection-coupled descriptions; they are not independent validation.

## Interpretation boundary

HUMAnN-style CPM is inferred gene-content potential. It is not a measurement
of NAD, NMN, nicotinamide, tryptophan, kynurenine, pathway activity, metabolic
flux, or causal mediation. PWY-5265 is not treated as a validated negative
control. If it associates with recovery, the result is interpreted as possible
shared richness/depth or compositional dependence.

The pathway-stratified audit is complemented by a direct gene-level check of
the Sankey claim for `nadD` (K00969), `nudC` (K03426), and `pncC` (the sum of
K03742 and K03743). These mappings follow the official KEGG nicotinate and
nicotinamide metabolism entries: https://www.kegg.jp/entry/ko00760 and
https://www.kegg.jp/entry/K03743. Total-share, assigned-share, and unresolved
remainder are kept separate so an apparently large fraction is not created by
renormalising only taxonomically assigned reads.

## Run

```bash
MPLCONFIGDIR=/tmp/pathway_mpl python reviewer_response/02_analysis/pathway_reanalysis/run_pathway_reanalysis.py
pytest -q reviewer_response/02_analysis/pathway_reanalysis/test_pathway_reanalysis.py
```

Outputs are written to `reviewer_response/03_tables/pathway_reanalysis/` and
`reviewer_response/04_figures/pathway_reanalysis/`.

The compact reviewer-facing KO attribution figure is exported as
`nad_gene_taxonomic_attribution.pdf`, `.svg`, and a 600-dpi LZW RGB `.tif`.
Its source-data table is `nad_gene_attribution_figure_source.csv`, derived
only from the audited subject-level KO contribution table. It reports all 60
participants; endpoint-defined HDR/LDR facets are deliberately omitted here
because that split is selection-coupled rather than an independent subgroup
validation.

## Figure contract

- Core conclusion: zero inflation and taxonomic attribution must be modelled
  explicitly before interpreting baseline pathway--recovery associations.
- Archetype: quantitative grid; the adjusted continuous-recovery forest is the
  hero evidence.
- Backend: Python only (matplotlib), including PDF/TIFF export and QA.
- Output: 7.2 x 6.8 inches, editable-text PDF and 600-dpi LZW RGB TIFF.
- Reviewer risk: endpoint-defined grouping, repeated-measure pseudoreplication,
  zero inflation, compositional confounding, unresolved taxonomy, and causal
  overstatement are all surfaced rather than hidden.

For the gene-level attribution figure, the denominator is each participant's
total unstratified KO CPM after A1--A4 arithmetic-mean aggregation. Component
points are subject-level medians with 95% non-parametric bootstrap confidence
intervals (5,000 resamples); grey points retain the participant distribution
for the combined Bacteroides + Phocaeicola share. Because component medians are
calculated separately, they are not forced to sum to 100%. The figure is a
7.2 x 3.8 inch quantitative grid with editable-text PDF/SVG and a 600-dpi LZW
RGB TIFF.

## Software-method citation

The statistical and figure workflow was informed by:

Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). Scientific Agent
Skills: A Library of Procedural Knowledge for Research Agents.
arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065
"""
    readme = paths.project_root / "reviewer_response" / "02_analysis" / "pathway_reanalysis" / "README.md"
    readme.write_text(text, encoding="utf-8")


def _jsonable(value: object) -> object:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_manifest(
    paths: Paths,
    figure_meta: dict[str, object],
    gene_figure_meta: dict[str, object],
    outputs: Sequence[Path],
    seed: int,
) -> None:
    inputs = [
        paths.raw / "meta.csv",
        paths.raw / "pathway_cpm_unstratified.tsv",
        paths.raw / "pathway_cpm_stratified.tsv",
        paths.raw / "ko_cpm_unstratified.tsv",
        paths.raw / "ko_cpm_stratified.tsv",
        paths.raw / "combine_brk.mpa",
        paths.project_root / "species.csv",
        paths.project_root / "reviewer_response" / "03_tables" / "core_subject_recovery_metrics.csv",
    ]
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed": seed,
        "sampling_unit": "subject",
        "baseline_visits": list(BASELINE),
        "endpoint_visits": list(ENDPOINT),
        "subjects": 60,
        "pathways": PATHWAYS,
        "claim_boundary": "inferred functional potential only; no NAD concentration, activity, flux, mediation, or mechanism claim",
        "primary_recovery": "Methods-faithful rarefied branch: 1 - mean C4-C6 Bray-Curtis distance to subject-specific A1-A4 baseline centroid",
        "primary_model": "association model with pathway response: two-part presence logistic GLM-HC3 plus positive log1p(CPM) OLS-HC3; rarefied Bacteroides/Phocaeicola and richness plus raw bacterial-read depth adjusted",
        "input_files": {str(p.relative_to(paths.project_root)): sha256(p) for p in inputs},
        "output_files": {str(p.relative_to(paths.project_root)): sha256(p) for p in outputs if p.exists()},
        "figure": figure_meta,
        "figures": {
            "pathway_summary": figure_meta,
            "nad_gene_taxonomic_attribution": gene_figure_meta,
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": importlib.metadata.version("statsmodels"),
            "matplotlib": mpl.__version__,
            "pillow": importlib.metadata.version("Pillow"),
        },
    }
    target = paths.tables / "analysis_manifest.json"
    target.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=_jsonable), encoding="utf-8")


def run(paths: Paths, seed: int = SEED) -> dict[str, Path]:
    paths.tables.mkdir(parents=True, exist_ok=True)
    paths.figures.mkdir(parents=True, exist_ok=True)
    baseline_meta, unstrat, strat = load_pathways(paths)
    full_meta = pd.read_csv(paths.raw / "meta.csv", dtype=str)
    covariates, phoc_epithets = load_taxonomy_and_recovery(paths, full_meta)
    subject_pathways, contributions, taxon_map = aggregate_pathways(
        baseline_meta, unstrat, strat, phoc_epithets
    )
    subject_data = build_subject_table(covariates, subject_pathways, contributions)
    zero_summary, contribution_summary = descriptive_tables(
        baseline_meta, unstrat, subject_pathways, contributions, seed
    )
    gene_contributions, gene_contribution_summary, gene_taxon_map = gene_level_contribution_audit(
        paths, baseline_meta, phoc_epithets, covariates, seed
    )
    gene_figure_source = build_gene_attribution_figure_source(gene_contributions, seed)
    models = run_hurdle_models(subject_data)
    aggregation_sensitivity = run_baseline_aggregation_sensitivity(
        baseline_meta, unstrat, covariates
    )
    models = pd.concat([models, aggregation_sensitivity], ignore_index=True)
    attribution = classify_attribution(models, contributions)

    outputs = {
        "subject_data": paths.tables / "subject_level_baseline_pathways.csv",
        "zero_summary": paths.tables / "pathway_zero_inflation_summary.csv",
        "models": paths.tables / "pathway_hurdle_models.csv",
        "aggregation_sensitivity": paths.tables / "baseline_aggregation_sensitivity.csv",
        "contributions": paths.tables / "stratified_contributions_by_subject.csv",
        "contribution_summary": paths.tables / "stratified_contribution_summary.csv",
        "taxon_map": paths.tables / "stratified_taxon_harmonization.csv",
        "gene_contributions": paths.tables / "nad_gene_stratified_contributions_by_subject.csv",
        "gene_contribution_summary": paths.tables / "nad_gene_stratified_contribution_summary.csv",
        "gene_taxon_map": paths.tables / "nad_gene_taxon_harmonization.csv",
        "gene_figure_source": paths.tables / "nad_gene_attribution_figure_source.csv",
        "attribution": paths.tables / "bacteroides_phocaeicola_attribution_judgment.csv",
    }
    subject_data.to_csv(outputs["subject_data"], index=False)
    zero_summary.to_csv(outputs["zero_summary"], index=False)
    models.to_csv(outputs["models"], index=False)
    aggregation_sensitivity.to_csv(outputs["aggregation_sensitivity"], index=False)
    contributions.to_csv(outputs["contributions"], index=False)
    contribution_summary.to_csv(outputs["contribution_summary"], index=False)
    taxon_map.to_csv(outputs["taxon_map"], index=False)
    gene_contributions.to_csv(outputs["gene_contributions"], index=False)
    gene_contribution_summary.to_csv(outputs["gene_contribution_summary"], index=False)
    gene_taxon_map.to_csv(outputs["gene_taxon_map"], index=False)
    gene_figure_source.to_csv(outputs["gene_figure_source"], index=False)
    attribution.to_csv(outputs["attribution"], index=False)

    figure_meta = render_summary_figure(subject_data, zero_summary, contribution_summary, models, paths, seed)
    gene_figure_meta = render_gene_attribution_figure(
        gene_contributions, gene_figure_source, paths, seed
    )
    write_readme(paths)
    output_paths = list(outputs.values()) + [
        paths.figures / "human_pathway_reanalysis_summary.pdf",
        paths.figures / "human_pathway_reanalysis_summary.tif",
        paths.figures / "nad_gene_taxonomic_attribution.pdf",
        paths.figures / "nad_gene_taxonomic_attribution.svg",
        paths.figures / "nad_gene_taxonomic_attribution.tif",
    ]
    write_manifest(paths, figure_meta, gene_figure_meta, output_paths, seed)
    outputs["manifest"] = paths.tables / "analysis_manifest.json"
    outputs["figure_pdf"] = paths.figures / "human_pathway_reanalysis_summary.pdf"
    outputs["figure_tif"] = paths.figures / "human_pathway_reanalysis_summary.tif"
    outputs["gene_figure_pdf"] = paths.figures / "nad_gene_taxonomic_attribution.pdf"
    outputs["gene_figure_svg"] = paths.figures / "nad_gene_taxonomic_attribution.svg"
    outputs["gene_figure_tif"] = paths.figures / "nad_gene_taxonomic_attribution.tif"
    return outputs


def default_paths(project_root: Path) -> Paths:
    response = project_root / "reviewer_response"
    return Paths(
        project_root=project_root,
        tables=response / "03_tables" / "pathway_reanalysis",
        figures=response / "04_figures" / "pathway_reanalysis",
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script = Path(__file__).resolve()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=script.parents[3])
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    project_root = args.project_root.resolve()
    outputs = run(default_paths(project_root), seed=args.seed)
    for name, path in outputs.items():
        print(f"{name}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
