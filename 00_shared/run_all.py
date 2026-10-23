#!/usr/bin/env python3
"""Reproduce the human taxonomic analyses supported by the supplied tables.

This script deliberately does not fabricate the missing HUMAnN/MetaCyc, MDSINE2,
FastSpar-bootstrap, or mouse data.  Approximate panels are labelled as such in the
figures and in ``outputs/reproduction_manifest.tsv``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from scipy.spatial.distance import pdist, squareform
from scipy.stats import (
    mannwhitneyu,
    pearsonr,
    rankdata,
    spearmanr,
    t as student_t,
    ttest_ind,
    ttest_rel,
)
from sklearn.cluster import KMeans


SEED = 20261009
TIME_ORDER = [f"A{i}" for i in range(1, 5)] + [f"B{i}" for i in range(1, 7)] + [f"C{i}" for i in range(1, 7)]
BASELINE = ["A1", "A2", "A3", "A4"]
EXPOSURE = ["B1", "B2", "B3", "B4", "B5"]
RECOVERY = ["B6", "C1", "C2", "C3", "C4", "C5", "C6"]
ENDPOINT = ["C4", "C5", "C6"]
PHASE_COLORS = {"baseline": "#DCEAF7", "exposure": "#F8E4DF", "recovery": "#E4F3E9"}
GROUP_COLORS = {"LDR": "#0072B2", "HDR": "#D55E00", "Middle": "#999999"}
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


@dataclass
class Dataset:
    root: Path
    meta: pd.DataFrame
    abundance: pd.DataFrame
    taxa: pd.DataFrame
    metrics: pd.DataFrame
    genus: pd.DataFrame
    species_named: pd.DataFrame
    group_table: pd.DataFrame


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjustment preserving NaNs."""
    values = np.asarray(list(pvalues), dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[valid] = restored
    return out


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan, np.nan
    result = spearmanr(x[mask], y[mask])
    return float(result.statistic), float(result.pvalue)


def strip_prefix(series: pd.Series) -> pd.Series:
    return series.fillna("Unclassified").astype(str).str.replace(r"^[a-z]__", "", regex=True).str.replace("_", " ")


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def add_phase_spans(ax: mpl.axes.Axes) -> None:
    ax.axvspan(-0.5, 3.5, color=PHASE_COLORS["baseline"], zorder=-10)
    ax.axvspan(3.5, 8.5, color=PHASE_COLORS["exposure"], zorder=-10)
    ax.axvspan(8.5, 15.5, color=PHASE_COLORS["recovery"], zorder=-10)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.10, 1.05, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def export_figure(fig: mpl.figure.Figure, stem: Path, dpi: int = 300) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white", transparent=False)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    # Matplotlib commonly writes an RGBA PNG even with an opaque white canvas.
    # Composite explicitly so the delivered raster has no latent alpha channel.
    with Image.open(png_path) as opened:
        rgba = opened.convert("RGBA")
        background = Image.new("RGBA", rgba.size, (255, 255, 255, 255))
        background.alpha_composite(rgba)
        background.convert("RGB").save(png_path, dpi=(dpi, dpi), optimize=True)


def load_dataset(root: Path) -> Dataset:
    meta = pd.read_csv(root / "meta.csv", dtype={"ids": str, "library": str, "subject": str, "canon": str})
    abundance = pd.read_csv(root / "species.csv", index_col=0)
    taxa = pd.read_csv(root / "taxa.csv", dtype=str).set_index("taxon")

    if meta["ids"].duplicated().any():
        raise ValueError("meta.csv contains duplicate ids")
    if meta["subject"].nunique() != 65 or not (meta.groupby("subject").size() == 16).all():
        raise ValueError("Expected 65 subjects with 16 samples each")
    if set(meta["ids"]) != set(abundance.columns):
        raise ValueError("meta ids and abundance columns differ")
    if set(abundance.index) != set(taxa.index):
        raise ValueError("abundance rows and taxonomy rows differ")
    if not np.isfinite(abundance.to_numpy()).all() or (abundance.to_numpy() < 0).any():
        raise ValueError("abundance matrix contains invalid values")

    abundance = abundance.loc[taxa.index, meta["ids"]]
    col_sums = abundance.sum(axis=0)
    prop = abundance.div(col_sums.replace(0, np.nan), axis=1)
    log_prop = np.log(prop.where(prop > 0))
    shannon = -(prop * log_prop).sum(axis=0)
    richness = (abundance > 0).sum(axis=0).astype(float)
    simpson = 1.0 - prop.pow(2).sum(axis=0)
    evenness = shannon / np.log(richness.where(richness > 1))

    metrics = meta.copy()
    metrics["day"] = pd.to_numeric(metrics["day"], errors="coerce")
    metrics["canons"] = pd.to_numeric(metrics["canons"], errors="coerce")
    metrics["shannon"] = metrics["ids"].map(shannon)
    metrics["richness"] = metrics["ids"].map(richness)
    metrics["simpson"] = metrics["ids"].map(simpson)
    metrics["evenness"] = metrics["ids"].map(evenness)
    metrics["canon"] = pd.Categorical(metrics["canon"], TIME_ORDER, ordered=True)

    # Paired presence/absence Jaccard distance to the same person's A1 sample.
    presence = abundance.to_numpy() > 0
    sample_index = {sample: i for i, sample in enumerate(abundance.columns)}
    paired_jaccard: dict[str, float] = {}
    bray_to_baseline_centroid: dict[str, float] = {}
    mean_jaccard_to_baseline: dict[str, float] = {}
    for _, subject_rows in metrics.groupby("subject", observed=True):
        a1_id = str(subject_rows.loc[subject_rows["canon"] == "A1", "ids"].iloc[0])
        a1 = presence[:, sample_index[a1_id]]
        baseline_ids = subject_rows.loc[subject_rows["canon"].isin(BASELINE), "ids"].tolist()
        baseline_prop = prop[baseline_ids].mean(axis=1).to_numpy()
        baseline_binary = [presence[:, sample_index[x]] for x in baseline_ids]
        for sample in subject_rows["ids"]:
            current_presence = presence[:, sample_index[sample]]
            union = np.logical_or(a1, current_presence).sum()
            paired_jaccard[sample] = 0.0 if union == 0 else 1.0 - np.logical_and(a1, current_presence).sum() / union
            current_prop = prop[sample].to_numpy()
            denom = np.abs(current_prop + baseline_prop).sum()
            bray_to_baseline_centroid[sample] = 0.0 if denom == 0 else np.abs(current_prop - baseline_prop).sum() / denom
            jd = []
            for baseline_presence in baseline_binary:
                u = np.logical_or(baseline_presence, current_presence).sum()
                jd.append(0.0 if u == 0 else 1.0 - np.logical_and(baseline_presence, current_presence).sum() / u)
            mean_jaccard_to_baseline[sample] = float(np.mean(jd))
    metrics["jaccard_to_A1"] = metrics["ids"].map(paired_jaccard)
    metrics["bray_to_baseline_centroid"] = metrics["ids"].map(bray_to_baseline_centroid)
    metrics["mean_jaccard_to_baseline"] = metrics["ids"].map(mean_jaccard_to_baseline)

    # Aggregate leaf taxa into named species and genera; re-normalize because names collapse strains.
    genus_labels = strip_prefix(taxa["g"])
    species_labels = strip_prefix(taxa["s"])
    genus = abundance.groupby(genus_labels, sort=False).sum()
    species_named = abundance.groupby(species_labels, sort=False).sum()

    exposed = metrics.loc[metrics["abx"] == "Y"]
    endpoint = exposed.loc[exposed["canon"].isin(ENDPOINT)].groupby("subject")["shannon"].mean().sort_values(kind="mergesort")
    baseline_mean = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject")["shannon"].mean()
    group_table = pd.DataFrame({"baseline_shannon": baseline_mean, "endpoint_shannon": endpoint})
    group_table["rank"] = group_table["endpoint_shannon"].rank(method="first").astype(int)
    group_table["recovery_group"] = "Middle"
    group_table.loc[endpoint.index[:15], "recovery_group"] = "LDR"
    group_table.loc[endpoint.index[-15:], "recovery_group"] = "HDR"
    group_table.index.name = "subject"
    metrics = metrics.merge(group_table[["recovery_group"]], left_on="subject", right_index=True, how="left")
    metrics["recovery_group"] = metrics["recovery_group"].fillna("Control")

    return Dataset(root, meta, abundance, taxa, metrics, genus, species_named, group_table)


def write_data_audit(ds: Dataset, tables: Path) -> None:
    tables.mkdir(parents=True, exist_ok=True)
    sample_sums = ds.abundance.sum(axis=0)
    audit = {
        "samples": int(ds.meta.shape[0]),
        "subjects": int(ds.meta["subject"].nunique()),
        "cip_subjects": int(ds.meta.loc[ds.meta["abx"] == "Y", "subject"].nunique()),
        "control_subjects": int(ds.meta.loc[ds.meta["abx"] == "N", "subject"].nunique()),
        "leaf_taxa": int(ds.abundance.shape[0]),
        "named_species": int(ds.species_named.shape[0]),
        "named_genera": int(ds.genus.shape[0]),
        "sample_sum_min": float(sample_sums.min()),
        "sample_sum_median": float(sample_sums.median()),
        "sample_sum_max": float(sample_sums.max()),
        "nonzero_fraction": float((ds.abundance.to_numpy() > 0).mean()),
        "max_possible_richness_from_supplied_table": int(ds.abundance.shape[0]),
        "observed_richness_min": int(ds.metrics["richness"].min()),
        "observed_richness_median": float(ds.metrics["richness"].median()),
        "observed_richness_max": int(ds.metrics["richness"].max()),
        "input_sha256": {
            name: hashlib.sha256((ds.root / name).read_bytes()).hexdigest()
            for name in ["meta.csv", "species.csv", "taxa.csv"]
        },
    }
    (tables / "data_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    ds.metrics.sort_values(["subject", "canons"]).to_csv(tables / "sample_metrics.csv", index=False)
    ds.group_table.sort_values("endpoint_shannon").to_csv(tables / "hdr_ldr_assignments.csv")

    day_rows = []
    for abx, group in ds.metrics.groupby("abx"):
        for canon, part in group.groupby("canon", observed=True):
            modes = part["day"].mode()
            day_rows.append(
                {
                    "abx": abx,
                    "canon": str(canon),
                    "n": len(part),
                    "day_min": part["day"].min(),
                    "day_median": part["day"].median(),
                    "day_mode": modes.iloc[0] if len(modes) else np.nan,
                    "day_max": part["day"].max(),
                }
            )
    pd.DataFrame(day_rows).to_csv(tables / "canonical_day_mapping.csv", index=False)

    # Presence-threshold sensitivity for the paired Jaccard trajectory.
    jaccard_sensitivity = []
    for threshold in [0.0, 0.001, 0.01]:
        presence = ds.abundance.to_numpy() > threshold
        sample_index = {sample: i for i, sample in enumerate(ds.abundance.columns)}
        for _, rows_for_subject in ds.metrics.loc[ds.metrics["abx"] == "Y"].groupby("subject", observed=True):
            a1_id = rows_for_subject.loc[rows_for_subject["canon"] == "A1", "ids"].iloc[0]
            a1 = presence[:, sample_index[a1_id]]
            for row in rows_for_subject.itertuples():
                current = presence[:, sample_index[row.ids]]
                union = np.logical_or(a1, current).sum()
                distance = 0.0 if union == 0 else 1.0 - np.logical_and(a1, current).sum() / union
                jaccard_sensitivity.append({"threshold_percent": threshold, "subject": row.subject, "canon": str(row.canon), "jaccard_to_A1": distance})
    sensitivity = pd.DataFrame(jaccard_sensitivity)
    sensitivity.groupby(["threshold_percent", "canon"])["jaccard_to_A1"].agg(["count", "mean", "std", "median"]).reset_index().to_csv(
        tables / "jaccard_presence_threshold_sensitivity.csv", index=False
    )

    rows = []
    for group in ["LDR", "HDR"]:
        part = ds.group_table.query("recovery_group == @group")
        sr = spearmanr(part["baseline_shannon"], part["endpoint_shannon"])
        pr = pearsonr(part["baseline_shannon"], part["endpoint_shannon"])
        rows.append(
            {
                "group": group,
                "n": len(part),
                "spearman_rho": sr.statistic,
                "spearman_p": sr.pvalue,
                "pearson_r": pr.statistic,
                "pearson_p": pr.pvalue,
            }
        )
    pd.DataFrame(rows).to_csv(tables / "baseline_endpoint_correlations.csv", index=False)


def paired_time_stats(ds: Dataset, metric: str, groups: tuple[str, ...] = ("all_exposed",)) -> pd.DataFrame:
    results = []
    selections: dict[str, pd.DataFrame] = {"all_exposed": ds.metrics.loc[ds.metrics["abx"] == "Y"]}
    selections.update({g: ds.metrics.loc[ds.metrics["recovery_group"] == g] for g in ("LDR", "HDR")})
    for group in groups:
        pivot = selections[group].pivot(index="subject", columns="canon", values=metric)
        pvals = []
        temp = []
        for time in TIME_ORDER:
            values = pivot[time].dropna()
            if time == "A1":
                stat, p = np.nan, np.nan
            else:
                aligned = pivot[["A1", time]].dropna()
                test = ttest_rel(aligned[time], aligned["A1"])
                stat, p = float(test.statistic), float(test.pvalue)
            pvals.append(p)
            temp.append(
                {
                    "group": group,
                    "metric": metric,
                    "canon": time,
                    "n": int(len(values)),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)),
                    "sem": float(values.sem()),
                    "t_vs_A1": stat,
                    "p_vs_A1": p,
                }
            )
        qvals = bh_adjust(pvals)
        for row, q in zip(temp, qvals):
            row["q_BH_vs_A1"] = q
            results.append(row)
    return pd.DataFrame(results)


def draw_boxplot(ax: mpl.axes.Axes, data: list[np.ndarray], colors: list[str], ylabel: str) -> None:
    bp = ax.boxplot(data, positions=np.arange(len(data)), widths=0.72, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor("#444444")
        patch.set_alpha(0.85)
    for key in ("whiskers", "caps", "medians"):
        for item in bp[key]:
            item.set_color("#444444")
            item.set_linewidth(0.8)
    ax.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.5)


def annotate_qvalues(ax: mpl.axes.Axes, stats: pd.DataFrame, y_frac: float = 0.97) -> None:
    y0, y1 = ax.get_ylim()
    y = y0 + y_frac * (y1 - y0)
    for i, row in stats.set_index("canon").reindex(TIME_ORDER).iterrows():
        q = row["q_BH_vs_A1"]
        if not np.isfinite(q) or i == "A1":
            continue
        label = "q<.001" if q < 0.001 else f"q={q:.3f}"
        ax.text(TIME_ORDER.index(i), y, label, rotation=90, va="top", ha="center", fontsize=5.2)


def figure1(ds: Dataset, figures: Path, tables: Path) -> None:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    s_stats = paired_time_stats(ds, "shannon")
    s_stats.to_csv(tables / "figure1_shannon_paired_tests.csv", index=False)

    # Jaccard inferential comparison is not identifiable from the manuscript wording.
    j_rows = []
    for time in TIME_ORDER:
        values = exposed.loc[exposed["canon"] == time, "jaccard_to_A1"]
        j_rows.append({"canon": time, "n": len(values), "mean": values.mean(), "sd": values.std(ddof=1), "median": values.median()})
    pd.DataFrame(j_rows).to_csv(tables / "figure1_jaccard_descriptive.csv", index=False)

    fig = plt.figure(figsize=(11.2, 8.0), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.42])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[:, 1])

    phase_box_colors = ["#76A5D3"] * 4 + ["#D96B5F"] * 5 + ["#A7D7A0"] * 7
    shannon_data = [exposed.loc[exposed["canon"] == t, "shannon"].to_numpy() for t in TIME_ORDER]
    draw_boxplot(ax_b, shannon_data, phase_box_colors, "Shannon index")
    means = exposed.groupby("canon", observed=True)["shannon"].mean().reindex(TIME_ORDER)
    ax_b.plot(range(16), means, color="#A52A2A", linestyle="--", marker="o", markersize=2.5, label="Mean trend")
    annotate_qvalues(ax_b, s_stats)
    ax_b.set_title("CIP-exposed cohort (n=60); paired t tests vs A1, BH corrected")
    panel_label(ax_b, "B")

    jac_data = [exposed.loc[exposed["canon"] == t, "jaccard_to_A1"].to_numpy() for t in TIME_ORDER]
    draw_boxplot(ax_c, jac_data, phase_box_colors, "Paired Jaccard distance to A1")
    jmeans = exposed.groupby("canon", observed=True)["jaccard_to_A1"].mean().reindex(TIME_ORDER)
    ax_c.plot(range(16), jmeans, color="#A52A2A", linestyle="--", marker="o", markersize=2.5)
    ax_c.set_title("Descriptive only: the article does not identify the Jaccard test contrast")
    panel_label(ax_c, "C")

    sample_ids = exposed.sort_values("canons").groupby("canon", observed=True)["ids"].apply(list)
    species_mean = pd.DataFrame({t: ds.species_named[sample_ids[t]].mean(axis=1) for t in TIME_ORDER})
    top = species_mean.mean(axis=1).nlargest(30).index
    plot_table = species_mean.loc[top]
    others = 100.0 - plot_table.sum(axis=0)
    others = others.clip(lower=0)
    plot_table.loc["Others"] = others
    plot_table.T.to_csv(tables / "figure1_top30_species_mean_percent.csv")
    cmap = plt.colormaps["tab20"].resampled(len(plot_table))
    colors = [cmap(i) for i in range(len(plot_table) - 1)] + ["#C8C8C8"]
    ax_d.stackplot(range(16), plot_table.to_numpy(), labels=plot_table.index, colors=colors, linewidth=0)
    ax_d.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
    ax_d.set_ylim(0, 100)
    ax_d.set_ylabel("Mean relative abundance (%)")
    ax_d.set_title("Top 30 named species; remaining abundance grouped as Others")
    ax_d.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), ncol=1, frameon=False, fontsize=5.5)
    panel_label(ax_d, "D")
    fig.suptitle("Figure 1 reproduction — panels supported by supplied human abundance data", fontsize=13)
    export_figure(fig, figures / "figure1_reproduced")


def differential_proxy(ds: Dataset, tables: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Time-stratified HDR/LDR differential abundance proxy.

    LEfSe's exact implementation is not operational in the supplied environment and
    its preprocessing is undocumented.  We therefore report signed log2 fold changes
    and Mann-Whitney P values with within-time BH adjustment, never calling these LDA.
    """
    results = []
    selected_subjects = ds.group_table.query("recovery_group in ['LDR','HDR']").index
    metric_subset = ds.metrics.loc[ds.metrics["subject"].isin(selected_subjects)]
    pseudocount = 1e-3
    selected_ids = metric_subset["ids"].tolist()
    prevalence = (ds.species_named[selected_ids] > 0).mean(axis=1)
    mean_abundance = ds.species_named[selected_ids].mean(axis=1)
    eligible = prevalence.index[(prevalence >= 0.20) & (mean_abundance >= 0.01)]
    for time in TIME_ORDER:
        rows = metric_subset.loc[metric_subset["canon"] == time]
        ldr_ids = rows.loc[rows["recovery_group"] == "LDR", "ids"].tolist()
        hdr_ids = rows.loc[rows["recovery_group"] == "HDR", "ids"].tolist()
        pvals, temp = [], []
        for feature in eligible:
            ldr = ds.species_named.loc[feature, ldr_ids].to_numpy(dtype=float)
            hdr = ds.species_named.loc[feature, hdr_ids].to_numpy(dtype=float)
            if np.all(ldr == 0) and np.all(hdr == 0):
                continue
            test = mannwhitneyu(hdr, ldr, alternative="two-sided", method="asymptotic")
            log2fc = math.log2((hdr.mean() + pseudocount) / (ldr.mean() + pseudocount))
            rank_biserial = 2.0 * float(test.statistic) / (len(hdr) * len(ldr)) - 1.0
            pvals.append(test.pvalue)
            temp.append({"canon": time, "species": feature, "log2FC_HDR_over_LDR": log2fc, "rank_biserial_HDR_over_LDR": rank_biserial, "p": test.pvalue})
        qvals = bh_adjust(pvals)
        for row, q in zip(temp, qvals):
            row["q_BH"] = q
            results.append(row)
    full = pd.DataFrame(results)
    full.to_csv(tables / "figure2_differential_species_proxy_all.csv", index=False)
    ranking = full.groupby("species")["rank_biserial_HDR_over_LDR"].apply(lambda x: x.abs().mean()).sort_values(ascending=False)
    top = ranking.head(30).index
    matrix = full.loc[full["species"].isin(top)].pivot(index="canon", columns="species", values="rank_biserial_HDR_over_LDR").reindex(TIME_ORDER)
    matrix = matrix.loc[:, ranking.loc[top].sort_values(ascending=False).index]
    return full, matrix


def figure2(ds: Dataset, figures: Path, tables: Path) -> None:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"]
    pivot = exposed.pivot(index="subject", columns="canon", values="shannon").reindex(columns=TIME_ORDER)
    corr = pivot.corr(method="spearman")
    corr.to_csv(tables / "figure2_timepoint_spearman.csv")
    km = KMeans(n_clusters=3, random_state=SEED, n_init=50).fit(1 - corr.to_numpy())
    pd.DataFrame({"canon": TIME_ORDER, "kmeans_cluster": km.labels_}).to_csv(tables / "figure2_timepoint_kmeans.csv", index=False)

    diff_full, diff_matrix = differential_proxy(ds, tables)
    sh_stats = []
    rich_stats = []
    for metric, sink in [("shannon", sh_stats), ("richness", rich_stats)]:
        for time in TIME_ORDER:
            ldr = ds.metrics.loc[(ds.metrics["recovery_group"] == "LDR") & (ds.metrics["canon"] == time), metric]
            hdr = ds.metrics.loc[(ds.metrics["recovery_group"] == "HDR") & (ds.metrics["canon"] == time), metric]
            test = ttest_ind(hdr, ldr, equal_var=False)
            sink.append({"metric": metric, "canon": time, "mean_LDR": ldr.mean(), "mean_HDR": hdr.mean(), "difference_HDR_minus_LDR": hdr.mean() - ldr.mean(), "welch_t": test.statistic, "p": test.pvalue})
        q = bh_adjust([x["p"] for x in sink])
        for row, value in zip(sink, q):
            row["q_BH"] = value
    pd.DataFrame(sh_stats + rich_stats).to_csv(tables / "figure2_group_comparisons.csv", index=False)

    fig = plt.figure(figsize=(12.2, 12.0), layout="constrained")
    gs = fig.add_gridspec(3, 6, height_ratios=[1.05, 0.9, 1.35])
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:])
    ax_f = fig.add_subplot(gs[2, :])

    im = ax_a.imshow(corr, vmin=0.25, vmax=1, cmap="RdBu_r", interpolation="nearest")
    ax_a.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
    ax_a.set_yticks(range(16), TIME_ORDER)
    ax_a.add_patch(mpl.patches.Rectangle((12.5, 12.5), 3, 3, fill=False, linestyle="--", linewidth=1.5, edgecolor="black"))
    fig.colorbar(im, ax=ax_a, label="Spearman rho", shrink=0.78)
    ax_a.set_title("Across-subject Shannon correlation")
    panel_label(ax_a, "A")

    rank_table = ds.group_table.sort_values("endpoint_shannon")
    endpoint_samples = exposed.loc[exposed["canon"].isin(ENDPOINT)].pivot(index="subject", columns="canon", values="shannon").loc[rank_table.index]
    x = np.arange(len(rank_table))
    width = 0.18
    endpoint_colors = ["#9BD59B", "#50B865", "#147A3B"]
    for j, time in enumerate(ENDPOINT):
        ax_b.bar(x + (j - 1.5) * width, endpoint_samples[time], width, color=endpoint_colors[j], label=time)
    ax_b.bar(x + 1.5 * width, rank_table["endpoint_shannon"], width, facecolor="white", edgecolor="black", linewidth=0.45, label="Mean")
    ax_b.axvspan(-0.8, 14.8, fill=False, edgecolor=GROUP_COLORS["LDR"], linestyle="--", linewidth=1.5)
    ax_b.axvspan(44.2, 59.8, fill=False, edgecolor=GROUP_COLORS["HDR"], linestyle="--", linewidth=1.5)
    ax_b.text(7, ax_b.get_ylim()[1] * 0.95, "LDR", color=GROUP_COLORS["LDR"], ha="center", fontweight="bold")
    ax_b.text(52, ax_b.get_ylim()[1] * 0.95, "HDR", color=GROUP_COLORS["HDR"], ha="center", fontweight="bold")
    ax_b.set(xlabel="Subjects ranked by C4–C6 mean", ylabel="Shannon index")
    ax_b.set_xticks([])
    ax_b.legend(ncol=4, frameon=False, loc="upper left")
    panel_label(ax_b, "B")

    correlation_rows = []
    for group in ["LDR", "HDR"]:
        part = rank_table.query("recovery_group == @group")
        color = GROUP_COLORS[group]
        ax_c.scatter(part["baseline_shannon"], part["endpoint_shannon"], color=color, label=group, s=24, alpha=0.9)
        coef = np.polyfit(part["baseline_shannon"], part["endpoint_shannon"], 1)
        xx = np.linspace(part["baseline_shannon"].min(), part["baseline_shannon"].max(), 50)
        ax_c.plot(xx, np.polyval(coef, xx), color=color, linestyle="--")
        rho, p = safe_spearman(part["baseline_shannon"].to_numpy(), part["endpoint_shannon"].to_numpy())
        correlation_rows.append({"group": group, "n": len(part), "rho": rho, "p": p})
    pd.DataFrame(correlation_rows).to_csv(tables / "figure2_panelC_correlations.csv", index=False)
    text_y = 0.98
    for row in correlation_rows:
        ax_c.text(0.02, text_y, f"{row['group']}: rho={row['rho']:.2f}, p={row['p']:.4g}", transform=ax_c.transAxes, color=GROUP_COLORS[row["group"]], va="top")
        text_y -= 0.09
    ax_c.set(xlabel="Baseline mean Shannon (A1–A4)", ylabel="Endpoint mean Shannon (C4–C6)")
    ax_c.legend(frameon=False, loc="lower right")
    panel_label(ax_c, "C")

    def trajectory(ax: mpl.axes.Axes, metric: str, title: str, ylabel: str) -> None:
        add_phase_spans(ax)
        for group in ["LDR", "HDR"]:
            d = ds.metrics.loc[ds.metrics["recovery_group"] == group]
            summary = d.groupby("canon", observed=True)[metric].agg(["mean", "sem"]).reindex(TIME_ORDER)
            xx = np.arange(16)
            ax.plot(xx, summary["mean"], marker="o", markersize=3, color=GROUP_COLORS[group], label=group)
            ax.fill_between(xx, summary["mean"] - summary["sem"], summary["mean"] + summary["sem"], color=GROUP_COLORS[group], alpha=0.18, linewidth=0)
        ax.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)

    trajectory(ax_d, "richness", "Observed leaf taxa (not original richness scale)", "Observed leaf taxa")
    panel_label(ax_d, "D")
    trajectory(ax_e, "shannon", "Shannon (mean ± SEM)", "Shannon index")
    panel_label(ax_e, "E")

    imf = ax_f.imshow(diff_matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax_f.set_yticks(range(16), TIME_ORDER)
    ax_f.set_xticks(range(len(diff_matrix.columns)), diff_matrix.columns, rotation=55, ha="right", fontsize=6)
    significant = diff_full.loc[(diff_full["q_BH"] < 0.05) & diff_full["species"].isin(diff_matrix.columns)]
    colpos = {c: i for i, c in enumerate(diff_matrix.columns)}
    rowpos = {c: i for i, c in enumerate(TIME_ORDER)}
    for row in significant.itertuples():
        ax_f.text(colpos[row.species], rowpos[row.canon], "*", ha="center", va="center", color="black", fontsize=6)
    fig.colorbar(imf, ax=ax_f, label="Rank-biserial effect (HDR > 0)", shrink=0.8)
    ax_f.set_title("Differential-abundance proxy: rank-biserial effect; Mann–Whitney + within-time BH (not LEfSe LDA)")
    panel_label(ax_f, "F")

    fig.suptitle("Figure 2 reproduction — current supplied table", fontsize=13)
    export_figure(fig, figures / "figure2_reproduced")


def feature_spearman(matrix: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Row-wise Spearman correlation with a common response."""
    ranked_x = rankdata(matrix, axis=1)
    ranked_y = rankdata(y)
    xc = ranked_x - ranked_x.mean(axis=1, keepdims=True)
    yc = ranked_y - ranked_y.mean()
    denom = np.sqrt(np.sum(xc * xc, axis=1) * np.sum(yc * yc))
    rho = np.divide(xc @ yc, denom, out=np.full(matrix.shape[0], np.nan), where=denom > 0)
    df = len(y) - 2
    tval = rho * np.sqrt(np.divide(df, np.maximum(1e-15, 1 - rho * rho)))
    p = 2 * student_t.sf(np.abs(tval), df)
    p[~np.isfinite(rho)] = np.nan
    return rho, p


def aggregate_metric_table(ds: Dataset, level: str) -> pd.DataFrame:
    if level == "genus":
        return ds.genus
    if level == "species":
        return ds.species_named
    raise ValueError(level)


def correlation_persistence(ds: Dataset, level: str, tables: Path) -> pd.DataFrame:
    feature_table = aggregate_metric_table(ds, level)
    rows = []
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])]
    for group in ["LDR", "HDR"]:
        for time in TIME_ORDER:
            meta = selected.loc[(selected["recovery_group"] == group) & (selected["canon"] == time)].sort_values("subject")
            matrix = feature_table[meta["ids"]].to_numpy(dtype=float)
            rho, p = feature_spearman(matrix, meta["shannon"].to_numpy(dtype=float))
            for feature, rv, pv in zip(feature_table.index, rho, p):
                rows.append({"level": level, "group": group, "canon": time, "feature": feature, "rho": rv, "p": pv, "positive_significant": bool(np.isfinite(pv) and pv < 0.05 and rv > 0)})
    out = pd.DataFrame(rows)
    out.to_csv(tables / f"suppfigure3_{level}_shannon_correlations.csv", index=False)
    return out


def phase_name(time: str) -> str:
    if time in BASELINE:
        return "A"
    if time in EXPOSURE:
        return "B"
    return "C"


def supplementary_figure3(ds: Dataset, figures: Path, tables: Path) -> None:
    corr_tables = {level: correlation_persistence(ds, level, tables) for level in ["genus", "species"]}
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 11.0), layout="constrained")

    for col, level in enumerate(["genus", "species"]):
        data = corr_tables[level].copy()
        data["phase"] = data["canon"].map(phase_name)
        phase_hit = data.groupby(["group", "feature", "phase"])["positive_significant"].any().unstack(fill_value=False)
        for phase in ["A", "B", "C"]:
            if phase not in phase_hit:
                phase_hit[phase] = False
        labels = []
        counts = {g: [] for g in ["LDR", "HDR"]}
        combos = [(True, False, False), (False, True, False), (False, False, True), (True, True, False), (True, False, True), (False, True, True), (True, True, True)]
        for combo in combos:
            labels.append("&".join([p for p, flag in zip(["A", "B", "C"], combo) if flag]))
            for group in ["LDR", "HDR"]:
                sub = phase_hit.loc[group]
                mask = np.ones(len(sub), dtype=bool)
                for phase, flag in zip(["A", "B", "C"], combo):
                    mask &= sub[phase].to_numpy() == flag
                counts[group].append(int(mask.sum()))
        ax = axes[0, col]
        x = np.arange(len(labels))
        ax.bar(x - 0.18, counts["LDR"], 0.36, color=GROUP_COLORS["LDR"], label="LDR")
        ax.bar(x + 0.18, counts["HDR"], 0.36, color=GROUP_COLORS["HDR"], label="HDR")
        ax.set_xticks(x, labels, rotation=45, ha="right")
        ax.set_ylabel(f"Number of {level} features")
        ax.set_title(f"Phase overlap — {level}")
        ax.legend(frameon=False)
        panel_label(ax, "A" if level == "genus" else "B")

        feature_counts = data.groupby(["group", "feature"])["positive_significant"].sum().reset_index()
        ax = axes[1, col]
        bins = np.arange(0.5, 16.6, 1)
        for group in ["LDR", "HDR"]:
            vals = feature_counts.loc[feature_counts["group"] == group, "positive_significant"]
            vals = vals[vals > 0]
            hist, _ = np.histogram(vals, bins=bins)
            ax.plot(range(1, 17), hist, marker="o", markersize=3, color=GROUP_COLORS[group], label=group)
        ax.set(xlabel="Number of significant-positive time points", ylabel=f"Number of {level} features", xticks=range(1, 17))
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False)
        panel_label(ax, "C" if level == "genus" else "D")

        sig = data.loc[data["positive_significant"]]
        patterns = sig.groupby(["group", "feature"])["canon"].apply(lambda x: "·".join(sorted(set(x), key=TIME_ORDER.index))).reset_index(name="pattern")
        pattern_counts = patterns.groupby(["group", "pattern"]).size().rename("n").reset_index()
        top_patterns = pattern_counts.groupby("pattern")["n"].sum().nlargest(18).index
        plot = pattern_counts.loc[pattern_counts["pattern"].isin(top_patterns)].pivot(index="pattern", columns="group", values="n").fillna(0)
        plot["total"] = plot.sum(axis=1)
        plot = plot.sort_values("total").drop(columns="total")
        ax = axes[2, col]
        y = np.arange(len(plot))
        ax.barh(y - 0.18, plot.get("LDR", 0), 0.36, color=GROUP_COLORS["LDR"], label="LDR")
        ax.barh(y + 0.18, plot.get("HDR", 0), 0.36, color=GROUP_COLORS["HDR"], label="HDR")
        ax.set_yticks(y, plot.index, fontsize=5.5)
        ax.set_xlabel(f"Number of {level} features")
        ax.set_title("Most frequent exact time-point combinations")
        ax.legend(frameon=False)
        panel_label(ax, "E" if level == "genus" else "F")

    fig.suptitle("Supplementary Figure 3 reproduction — raw P<0.05 and rho>0, as stated", fontsize=13)
    export_figure(fig, figures / "supplementary_figure3_reproduced")


def pcoa(distance_matrix: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray, float]:
    n = distance_matrix.shape[0]
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ (distance_matrix ** 2) @ h
    values, vectors = np.linalg.eigh(b)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    positive = values > 0
    coords = vectors[:, positive][:, :n_components] * np.sqrt(values[positive][:n_components])
    explained = values[positive][:n_components] / values[positive].sum()
    negative_fraction = float(np.abs(values[values < 0]).sum() / values[positive].sum()) if positive.any() else np.nan
    return coords, explained, negative_fraction


def permanova(
    distance_matrix: np.ndarray,
    labels: np.ndarray,
    permutations: int = 999,
    seed: int = SEED,
    strata: np.ndarray | None = None,
) -> tuple[float, float]:
    """One-factor PERMANOVA using Anderson's distance-based sums of squares."""
    labels = np.asarray(labels)
    n = len(labels)
    groups = np.unique(labels)
    tri = np.triu_indices(n, 1)
    ss_total = np.sum(distance_matrix[tri] ** 2) / n

    def pseudo_f(current: np.ndarray) -> float:
        within = 0.0
        for group in groups:
            idx = np.flatnonzero(current == group)
            if len(idx) > 1:
                sub = distance_matrix[np.ix_(idx, idx)]
                within += np.sum(sub[np.triu_indices(len(idx), 1)] ** 2) / len(idx)
        among = ss_total - within
        return (among / (len(groups) - 1)) / (within / (n - len(groups)))

    observed = pseudo_f(labels)
    rng = np.random.default_rng(seed)
    if strata is None:
        shuffled = lambda: rng.permutation(labels)
    else:
        strata = np.asarray(strata)
        unique_strata = np.unique(strata)

        def shuffled() -> np.ndarray:
            result = labels.copy()
            for stratum in unique_strata:
                idx = np.flatnonzero(strata == stratum)
                result[idx] = rng.permutation(result[idx])
            return result

    permuted = np.array([pseudo_f(shuffled()) for _ in range(permutations)])
    p = (1 + np.sum(permuted >= observed)) / (permutations + 1)
    return float(observed), float(p)


def supplementary_figure1(ds: Dataset, figures: Path, tables: Path) -> None:
    control = ds.metrics.loc[ds.metrics["abx"] == "N"]
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"]
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.8), layout="constrained")

    phase_colors = ["#76A5D3"] * 4 + ["#D96B5F"] * 5 + ["#A7D7A0"] * 7
    draw_boxplot(axes[0, 0], [control.loc[control["canon"] == t, "shannon"].to_numpy() for t in TIME_ORDER], phase_colors, "Shannon index")
    axes[0, 0].set_title("Unexposed controls (n=5); descriptive")
    panel_label(axes[0, 0], "A")

    draw_boxplot(axes[0, 1], [control.loc[control["canon"] == t, "jaccard_to_A1"].to_numpy() for t in TIME_ORDER], phase_colors, "Paired Jaccard distance to A1")
    axes[0, 1].set_title("Unexposed controls (n=5); descriptive")
    panel_label(axes[0, 1], "B")

    ids_by_time = exposed.groupby("canon", observed=True)["ids"].apply(list)
    prevalence = pd.DataFrame({t: (ds.species_named[ids_by_time[t]] > 0.01).mean(axis=1) for t in TIME_ORDER})
    overall_prev = (ds.species_named[exposed["ids"]] > 0.01).mean(axis=1)
    core = overall_prev.loc[overall_prev > 0.5].sort_values(ascending=False).head(30).index
    prevalence.loc[core].to_csv(tables / "suppfigure1_core_species_prevalence.csv")
    im = axes[1, 0].imshow(prevalence.loc[core].to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    axes[1, 0].set_yticks(range(len(core)), core, fontsize=5.5)
    axes[1, 0].set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
    axes[1, 0].set_title("Core species prevalence (RA>0.01% in >50% of samples)")
    fig.colorbar(im, ax=axes[1, 0], label="Prevalence", shrink=0.75)
    panel_label(axes[1, 0], "C")

    baseline_ids = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject")["ids"].apply(list)
    recovery_ids = exposed.loc[exposed["canon"].isin(RECOVERY)].groupby("subject")["ids"].apply(list)
    paired_rows = []
    for feature in ds.species_named.index:
        before = np.array([ds.species_named.loc[feature, ids].mean() for ids in baseline_ids])
        after = np.array([ds.species_named.loc[feature, recovery_ids.loc[s]].mean() for s in baseline_ids.index])
        test = ttest_rel(after, before)
        paired_rows.append({"species": feature, "baseline_mean": before.mean(), "recovery_mean": after.mean(), "difference": after.mean() - before.mean(), "p": test.pvalue})
    paired = pd.DataFrame(paired_rows)
    paired["q_BH"] = bh_adjust(paired["p"])
    paired.to_csv(tables / "suppfigure1_phaseA_vs_phaseC_species.csv", index=False)
    article_species = [
        "Alistipes onderdonkii",
        "Bifidobacterium adolescentis",
        "Blautia wexlerae",
        "Faecalibacterium prausnitzii",
        "Parabacteroides distasonis",
        "Eubacterium ramulus",
        "Ruminococcus gnavus",
        "Bacteroides ovatus",
        "Bacteroides xylanisolvens",
    ]
    present_article_species = [x for x in article_species if x in set(paired["species"])]
    if len(present_article_species) == 9:
        top = paired.set_index("species").loc[present_article_species].reset_index()
    else:
        top = paired.assign(absdiff=lambda x: x["difference"].abs()).sort_values(["q_BH", "absdiff"], ascending=[True, False]).head(9)
    y = np.arange(len(top))
    axes[1, 1].barh(y - 0.18, top["baseline_mean"], 0.36, color="#4C78A8", label="Baseline A1–A4")
    axes[1, 1].barh(y + 0.18, top["recovery_mean"], 0.36, color="#59A14F", label="Recovery B6–C6")
    axes[1, 1].set_yticks(y, top["species"], fontsize=6)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("Subject-mean relative abundance (%)")
    axes[1, 1].set_title("Largest paired phase shifts (subject is the unit)")
    axes[1, 1].legend(frameon=False)
    panel_label(axes[1, 1], "D")

    fig.suptitle("Supplementary Figure 1 reproduction", fontsize=13)
    export_figure(fig, figures / "supplementary_figure1_reproduced")


def supplementary_figure2(ds: Dataset, figures: Path, tables: Path) -> None:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])].copy()
    fig = plt.figure(figsize=(12.0, 10.2), layout="constrained")
    gs = fig.add_gridspec(3, 4)
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])
    ax_e = fig.add_subplot(gs[2, :])

    phase_colors = ["#76A5D3"] * 4 + ["#D96B5F"] * 5 + ["#A7D7A0"] * 7
    draw_boxplot(ax_a, [exposed.loc[exposed["canon"] == t, "jaccard_to_A1"].to_numpy() for t in TIME_ORDER], phase_colors, "Paired Jaccard distance to A1")
    ax_a.set_title("All CIP-exposed participants (n=60)")
    panel_label(ax_a, "A")

    baseline = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject")["shannon"].mean()
    recovery_pivot = exposed.loc[exposed["canon"].isin(RECOVERY)].pivot(index="subject", columns="canon", values="shannon").reindex(columns=RECOVERY)
    deviations = recovery_pivot.sub(baseline, axis=0)
    classes = pd.Series(index=deviations.index, dtype=str)
    classes[(deviations > 0).all(axis=1)] = "Always above baseline"
    classes[(deviations < 0).all(axis=1)] = "Always below baseline"
    classes = classes.fillna("Crosses baseline")
    dev_out = deviations.copy()
    dev_out.insert(0, "trajectory_class", classes)
    dev_out.to_csv(tables / "suppfigure2_recovery_shannon_deviations.csv")
    class_colors = {"Always above baseline": "#009E73", "Always below baseline": "#D55E00", "Crosses baseline": "#0072B2"}
    for subject, row in deviations.iterrows():
        ax_b.plot(range(7), row, color=class_colors[classes[subject]], alpha=0.28, linewidth=0.8)
    for name, color in class_colors.items():
        mean = deviations.loc[classes == name].mean()
        ax_b.plot(range(7), mean, color=color, linewidth=2.2, label=f"{name} (n={(classes == name).sum()})")
    ax_b.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax_b.set_xticks(range(7), RECOVERY, rotation=45, ha="right")
    ax_b.set_ylabel("Shannon minus subject baseline mean")
    ax_b.legend(frameon=False, fontsize=6)
    ax_b.set_title("Recovery trajectories; operational three-class rule documented")
    panel_label(ax_b, "B")

    # Jaccard PCoA for all 960 exposed samples.
    ids = exposed.sort_values(["canons", "subject"])["ids"].tolist()
    pmat = (ds.abundance[ids].to_numpy().T > 0).astype(bool)
    distances = squareform(pdist(pmat, metric="jaccard"))
    coords, explained, negative_fraction = pcoa(distances)
    pcoa_meta = exposed.set_index("ids").loc[ids].reset_index()
    pcoa_meta["PCo1"] = coords[:, 0]
    pcoa_meta["PCo2"] = coords[:, 1]
    pcoa_meta["phase"] = pcoa_meta["canon"].map(lambda x: phase_name(str(x)))
    pcoa_meta.to_csv(tables / "suppfigure2_jaccard_pcoa_scores.csv", index=False)
    phase_palette = {"A": "#4C78A8", "B": "#E45756", "C": "#59A14F"}
    for phase in ["A", "B", "C"]:
        part = pcoa_meta.loc[pcoa_meta["phase"] == phase]
        ax_c.scatter(part["PCo1"], part["PCo2"], s=7, alpha=0.35, color=phase_palette[phase], label=phase)
        ax_c.scatter(part["PCo1"].mean(), part["PCo2"].mean(), s=70, marker="X", color=phase_palette[phase], edgecolor="black")
    fstat, pp = permanova(
        distances,
        pcoa_meta["phase"].to_numpy(),
        permutations=999,
        strata=pcoa_meta["subject"].to_numpy(),
    )
    ax_c.set(xlabel=f"PCo1 ({explained[0]*100:.1f}%)", ylabel=f"PCo2 ({explained[1]*100:.1f}%)")
    ax_c.set_title(f"Jaccard PCoA; subject-restricted PERMANOVA F={fstat:.2f}, p={pp:.3f}")
    ax_c.legend(frameon=False)
    panel_label(ax_c, "C")

    selected_ids = selected.sort_values(["canons", "subject"])["ids"].tolist()
    selected_presence = (ds.abundance[selected_ids].to_numpy().T > 0).astype(bool)
    selected_dist = squareform(pdist(selected_presence, metric="jaccard"))
    selected_coords, selected_exp, selected_negative_fraction = pcoa(selected_dist)
    selected_meta = selected.set_index("ids").loc[selected_ids].reset_index()
    selected_meta["PCo1"] = selected_coords[:, 0]
    selected_meta["PCo2"] = selected_coords[:, 1]
    centroids = selected_meta.groupby(["recovery_group", "canon"], observed=True)[["PCo1", "PCo2"]].mean().reset_index()
    for group in ["LDR", "HDR"]:
        part = centroids.loc[centroids["recovery_group"] == group].set_index("canon").reindex(TIME_ORDER)
        ax_d.plot(part["PCo1"], part["PCo2"], marker="o", markersize=4, color=GROUP_COLORS[group], label=group)
        for time in ["A1", "B1", "B6", "C6"]:
            ax_d.text(part.loc[time, "PCo1"], part.loc[time, "PCo2"], time, fontsize=6)
    ax_d.set(xlabel=f"PCo1 ({selected_exp[0]*100:.1f}%)", ylabel=f"PCo2 ({selected_exp[1]*100:.1f}%)")
    ax_d.set_title("HDR/LDR centroid trajectories")
    ax_d.legend(frameon=False)
    panel_label(ax_d, "D")

    target = "Faecalibacterium prausnitzii"
    if target not in ds.species_named.index:
        candidates = [x for x in ds.species_named.index if "Faecalibacterium prausnitzii" in x]
        target = candidates[0] if candidates else ds.species_named.mean(axis=1).idxmax()
    for group in ["LDR", "HDR"]:
        part = selected.loc[selected["recovery_group"] == group]
        values = []
        sems = []
        for time in TIME_ORDER:
            sample_ids = part.loc[part["canon"] == time, "ids"]
            val = ds.species_named.loc[target, sample_ids]
            values.append(val.mean())
            sems.append(val.sem())
        x = np.arange(16)
        values = np.asarray(values)
        sems = np.asarray(sems)
        ax_e.plot(x, values, marker="o", color=GROUP_COLORS[group], label=group)
        ax_e.fill_between(x, values - sems, values + sems, color=GROUP_COLORS[group], alpha=0.18)
    add_phase_spans(ax_e)
    ax_e.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
    ax_e.set_ylabel("Relative abundance (%)")
    ax_e.set_title(target)
    ax_e.legend(frameon=False)
    panel_label(ax_e, "E")
    (tables / "suppfigure2_pcoa_diagnostics.json").write_text(
        json.dumps(
            {
                "all_exposed_negative_eigenvalue_fraction_of_positive_sum": negative_fraction,
                "HDR_LDR_negative_eigenvalue_fraction_of_positive_sum": selected_negative_fraction,
                "variance_percent_uses_positive_eigenvalues_only": True,
                "PERMANOVA_permutations": 999,
                "PERMANOVA_restriction": "phase labels permuted within subject",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    fig.suptitle("Supplementary Figure 2 reproduction", fontsize=13)
    export_figure(fig, figures / "supplementary_figure2_reproduced")


def network_proxy(ds: Dataset, group: str, stage: str, max_features: int = 40) -> tuple[nx.Graph, pd.DataFrame]:
    subset = ds.metrics.loc[ds.metrics["recovery_group"] == group]
    if stage == "Baseline":
        sample_rows = subset.loc[subset["canon"].isin(BASELINE)]
    elif stage == "End":
        sample_rows = subset.loc[subset["canon"].isin(ENDPOINT)]
    else:
        sample_rows = subset.loc[subset["canon"] == stage]
    subject_tables = []
    for subject, rows in sample_rows.groupby("subject"):
        subject_tables.append(ds.genus[rows["ids"]].mean(axis=1).rename(subject))
    data = pd.concat(subject_tables, axis=1)
    features = data.mean(axis=1).nlargest(max_features).index
    data = data.loc[features]
    rho = data.T.corr(method="spearman")
    p_rows = []
    graph = nx.Graph()
    graph.add_nodes_from(features)
    for i, first in enumerate(features):
        for second in features[i + 1 :]:
            r, p = safe_spearman(data.loc[first].to_numpy(), data.loc[second].to_numpy())
            p_rows.append({"group": group, "stage": stage, "genus1": first, "genus2": second, "rho": r, "p": p})
            if np.isfinite(p) and abs(r) > 0.6 and p < 0.05:
                graph.add_edge(first, second, weight=abs(r), sign=np.sign(r))
    graph.remove_nodes_from(list(nx.isolates(graph)))
    return graph, pd.DataFrame(p_rows)


def figure4(ds: Dataset, figures: Path, tables: Path) -> None:
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])]
    fig = plt.figure(figsize=(13.0, 17.5), layout="constrained")
    gs = fig.add_gridspec(5, 6, height_ratios=[1.0, 0.72, 0.82, 0.82, 1.05])
    ax_a1 = fig.add_subplot(gs[0, :3])
    ax_a2 = fig.add_subplot(gs[0, 3:])
    ax_b = fig.add_subplot(gs[1, :2])
    ax_nodes = fig.add_subplot(gs[1, 2:4])
    ax_edges = fig.add_subplot(gs[1, 4:])
    network_axes = [fig.add_subplot(gs[2, i]) for i in range(6)] + [fig.add_subplot(gs[3, i]) for i in range(6)]
    ax_d1 = fig.add_subplot(gs[4, :3])
    ax_d2 = fig.add_subplot(gs[4, 3:])

    # A: stacked genus compositions (honest approximation to alluvial flow).
    group_compositions = {}
    selected_genera = set()
    for group in ["LDR", "HDR"]:
        group_rows = selected.loc[selected["recovery_group"] == group]
        comp = pd.DataFrame({t: ds.genus[group_rows.loc[group_rows["canon"] == t, "ids"]].mean(axis=1) for t in TIME_ORDER})
        group_compositions[group] = comp
        selected_genera.update(comp.mean(axis=1).nlargest(16).index)
    selected_genera = sorted(selected_genera, key=lambda g: -np.mean([group_compositions[x].loc[g].mean() if g in group_compositions[x].index else 0 for x in ["LDR", "HDR"]]))[:16]
    palette = [plt.colormaps["tab20"](i) for i in range(16)]
    for ax, group in [(ax_a1, "LDR"), (ax_a2, "HDR")]:
        comp = group_compositions[group]
        plotted = comp.reindex(selected_genera).fillna(0)
        plotted.loc["Others"] = (100 - plotted.sum(axis=0)).clip(lower=0)
        plotted.T.to_csv(tables / f"figure4_genus_composition_{group}.csv")
        ax.stackplot(range(16), plotted.to_numpy(), labels=plotted.index, colors=palette + ["#C8C8C8"], linewidth=0)
        ax.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
        ax.set_ylim(0, 100)
        ax.set_title(f"{group}: mean genus composition")
        ax.set_ylabel("Relative abundance (%)")
        if group == "HDR":
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=5.5)
    panel_label(ax_a1, "A")

    # B: metric heat maps, columns HDR/LDR, rows time.
    matrices = []
    for metric in ["shannon", "richness", "simpson"]:
        matrix = selected.groupby(["canon", "recovery_group"], observed=True)[metric].mean().unstack().reindex(TIME_ORDER)[["HDR", "LDR"]]
        scaled = (matrix - matrix.min().min()) / (matrix.max().max() - matrix.min().min())
        matrices.append(scaled.to_numpy())
    heat = np.concatenate(matrices, axis=1)
    ax_b.imshow(heat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax_b.set_yticks(range(16), TIME_ORDER)
    ax_b.set_xticks(
        range(6),
        ["Shannon\nHDR", "Shannon\nLDR", "Richness*\nHDR", "Richness*\nLDR", "Simpson\nHDR", "Simpson\nLDR"],
        rotation=45,
        ha="right",
    )
    for x in [1.5, 3.5]:
        ax_b.axvline(x, color="white", linewidth=2)
    ax_b.set_title("Mean diversity; within-metric 0–1 scale\n*richness cannot match original")
    panel_label(ax_b, "B")

    # C: six stage networks per group using a documented Spearman proxy.
    stages = ["Baseline", "B6", "C1", "C2", "C3", "End"]
    network_stats = []
    pair_stats = []
    for row_group, group in enumerate(["HDR", "LDR"]):
        for stage_idx, stage in enumerate(stages):
            axis_index = stage_idx if row_group == 0 else 6 + stage_idx
            ax = network_axes[axis_index]
            graph, pairs = network_proxy(ds, group, stage)
            pair_stats.append(pairs)
            pos = nx.spring_layout(graph, seed=SEED, weight="weight")
            degrees = dict(graph.degree())
            nx.draw_networkx_edges(graph, pos, ax=ax, width=0.4, alpha=0.35, edge_color=["#555555" if graph.edges[e]["sign"] > 0 else "#B2182B" for e in graph.edges])
            nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=[8 + 3 * degrees[n] for n in graph.nodes], node_color="#56B4E9" if group == "LDR" else "#D55E00", alpha=0.8, linewidths=0)
            ax.set_title(f"{group} {stage}", fontsize=7)
            ax.set_axis_off()
            network_stats.append({"group": group, "stage": stage, "nodes": graph.number_of_nodes(), "edges": graph.number_of_edges(), "method": "top40 genus subject-mean Spearman proxy; |rho|>0.6, raw p<0.05"})
    pd.DataFrame(network_stats).to_csv(tables / "figure4_network_proxy_summary.csv", index=False)
    pd.concat(pair_stats, ignore_index=True).to_csv(tables / "figure4_network_proxy_all_pairs.csv", index=False)
    network_summary = pd.DataFrame(network_stats)
    for group in ["LDR", "HDR"]:
        part = network_summary.loc[network_summary["group"] == group].set_index("stage").reindex(stages)
        ax_nodes.plot(range(6), part["nodes"], marker="o", color=GROUP_COLORS[group], label=group)
        ax_edges.plot(range(6), part["edges"], marker="o", color=GROUP_COLORS[group], label=group)
    for ax, ylabel in [(ax_nodes, "Number of nodes"), (ax_edges, "Number of edges")]:
        ax.set_xticks(range(6), stages, rotation=45, ha="right")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False)
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.5)
    panel_label(network_axes[0], "C")

    # D: Bacteroides/Faecalibacterium correlations with three diversity indices.
    correlation_rows = []
    for ax, group in [(ax_d1, "LDR"), (ax_d2, "HDR")]:
        group_meta = selected.loc[selected["recovery_group"] == group]
        for genus, color, marker in [("Bacteroides", "#0072B2", "o"), ("Faecalibacterium", "#D55E00", "s")]:
            values = []
            for time in TIME_ORDER:
                rows = group_meta.loc[group_meta["canon"] == time].sort_values("subject")
                abundance = ds.genus.loc[genus, rows["ids"]].to_numpy() if genus in ds.genus.index else np.zeros(len(rows))
                rho, p = safe_spearman(abundance, rows["shannon"].to_numpy())
                values.append(rho)
                correlation_rows.append({"group": group, "canon": time, "genus": genus, "metric": "shannon", "rho": rho, "p": p})
            ax.plot(range(16), values, color=color, marker=marker, markersize=3, label=genus)
        add_phase_spans(ax)
        ax.axhline(0, color="#777777", linewidth=0.6)
        ax.set_xticks(range(16), TIME_ORDER, rotation=45, ha="right")
        ax.set_ylim(-1, 1)
        ax.set_ylabel("Spearman rho with Shannon")
        ax.set_title(group)
        ax.legend(frameon=False)
    pd.DataFrame(correlation_rows).to_csv(tables / "figure4_dominant_genus_correlations.csv", index=False)
    panel_label(ax_d1, "D")
    fig.suptitle("Figure 4 reproduction — FastSpar networks replaced by labelled Spearman proxies", fontsize=13)
    export_figure(fig, figures / "figure4_reproduced")


def reviewer_reanalysis(ds: Dataset, figures: Path, tables: Path) -> None:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    base = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject").agg(
        baseline_shannon=("shannon", "mean"),
        baseline_shannon_sd=("shannon", "std"),
    )
    end = exposed.loc[exposed["canon"].isin(ENDPOINT)].groupby("subject").agg(
        endpoint_shannon=("shannon", "mean"),
        endpoint_mean_bray_to_baseline_centroid=("bray_to_baseline_centroid", "mean"),
        endpoint_mean_jaccard_to_baseline_samples=("mean_jaccard_to_baseline", "mean"),
    )
    subject = base.join(end)
    subject["signed_shannon_change"] = subject["endpoint_shannon"] - subject["baseline_shannon"]
    subject["absolute_shannon_displacement"] = subject["signed_shannon_change"].abs()

    rec = exposed.loc[exposed["canon"].isin(RECOVERY)].pivot(index="subject", columns="canon", values="shannon").reindex(columns=RECOVERY)
    dev = rec.sub(subject["baseline_shannon"], axis=0)
    subject["mean_abs_recovery_displacement"] = dev.abs().mean(axis=1)
    subject["auc_abs_recovery_displacement_canonical"] = np.trapz(dev.abs().to_numpy(), x=np.arange(len(RECOVERY)), axis=1)

    # Composition-centroid distances, actual-day AUCs and an explicit return rule.
    prop = ds.abundance.div(ds.abundance.sum(axis=0), axis=1)
    structural_rows = []
    for person, person_rows in exposed.groupby("subject"):
        baseline_rows = person_rows.loc[person_rows["canon"].isin(BASELINE)]
        endpoint_rows = person_rows.loc[person_rows["canon"].isin(ENDPOINT)]
        recovery_rows = person_rows.loc[person_rows["canon"].isin(RECOVERY)].sort_values("day")
        baseline_centroid = prop[baseline_rows["ids"]].mean(axis=1).to_numpy()
        endpoint_centroid = prop[endpoint_rows["ids"]].mean(axis=1).to_numpy()
        baseline_centroid /= baseline_centroid.sum()
        endpoint_centroid /= endpoint_centroid.sum()
        centroid_bray = np.abs(baseline_centroid - endpoint_centroid).sum() / np.abs(baseline_centroid + endpoint_centroid).sum()

        baseline_binary = (ds.abundance[baseline_rows["ids"]].to_numpy() > 0).sum(axis=1) >= 2
        endpoint_binary = (ds.abundance[endpoint_rows["ids"]].to_numpy() > 0).sum(axis=1) >= 2
        union = np.logical_or(baseline_binary, endpoint_binary).sum()
        centroid_jaccard = 0.0 if union == 0 else 1.0 - np.logical_and(baseline_binary, endpoint_binary).sum() / union

        actual_days = recovery_rows["day"].to_numpy(dtype=float)
        shannon_abs_dev = np.abs(recovery_rows["shannon"].to_numpy(dtype=float) - subject.loc[person, "baseline_shannon"])
        bray_dev = recovery_rows["bray_to_baseline_centroid"].to_numpy(dtype=float)
        time_span = actual_days[-1] - actual_days[0]
        shannon_auc = np.trapz(shannon_abs_dev, x=actual_days) / time_span if time_span > 0 else np.nan
        bray_auc = np.trapz(bray_dev, x=actual_days) / time_span if time_span > 0 else np.nan

        band = subject.loc[person, "baseline_shannon_sd"]
        return_index = np.nan
        return_day = np.nan
        for i in range(len(recovery_rows)):
            if (shannon_abs_dev[i:] <= band).all():
                return_index = i
                return_day = actual_days[i]
                break
        structural_rows.append(
            {
                "subject": person,
                "endpoint_centroid_bray": centroid_bray,
                "endpoint_binary_centroid_jaccard": centroid_jaccard,
                "shannon_auc_actual_day_normalized": shannon_auc,
                "bray_auc_actual_day_normalized": bray_auc,
                "persistent_return_stage_index": return_index,
                "persistent_return_actual_day": return_day,
            }
        )
    subject = subject.join(pd.DataFrame(structural_rows).set_index("subject"))

    # Baseline dominant genus; keep strict and renamed-inclusive Bacteroides definitions.
    dominant = {}
    bact_fraction = {}
    bp_fraction = {}
    bp_dominant = {}
    for person, rows in exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject"):
        mean_genus = ds.genus[rows["ids"]].mean(axis=1)
        dominant[person] = mean_genus.idxmax()
        bact_fraction[person] = float(mean_genus.get("Bacteroides", 0.0))
        bp = float(mean_genus.get("Bacteroides", 0.0) + mean_genus.get("Phocaeicola", 0.0))
        bp_fraction[person] = bp
        competitors = mean_genus.drop(labels=[x for x in ["Bacteroides", "Phocaeicola"] if x in mean_genus.index])
        bp_dominant[person] = bool(bp > competitors.max())
    subject["dominant_genus"] = pd.Series(dominant)
    subject["bacteroides_dominant"] = subject["dominant_genus"].eq("Bacteroides")
    subject["baseline_bacteroides_percent"] = pd.Series(bact_fraction)
    subject["baseline_bacteroides_plus_phocaeicola_percent"] = pd.Series(bp_fraction)
    subject["bacteroides_plus_phocaeicola_dominant"] = pd.Series(bp_dominant)
    subject = subject.join(ds.group_table[["recovery_group", "rank"]])
    subject.to_csv(tables / "reviewer_subject_level_recovery_metrics.csv")

    # Continuous associations: Spearman plus OLS/HC3 implemented with closed-form matrix algebra.
    outcomes = [
        "endpoint_shannon",
        "absolute_shannon_displacement",
        "endpoint_centroid_bray",
        "endpoint_binary_centroid_jaccard",
        "shannon_auc_actual_day_normalized",
        "bray_auc_actual_day_normalized",
        "mean_abs_recovery_displacement",
    ]
    model_rows = []
    for outcome in outcomes:
        for model_name, columns in [
            ("unadjusted", ["baseline_shannon"]),
            ("adjusted_Bacteroides", ["baseline_shannon", "baseline_bacteroides_percent"]),
            ("adjusted_BacteroidesPlusPhocaeicola", ["baseline_shannon", "baseline_bacteroides_plus_phocaeicola_percent"]),
        ]:
            analysis = subject[[outcome] + columns].dropna()
            y = analysis[outcome].to_numpy(dtype=float)
            x = np.column_stack([np.ones(len(analysis))] + [analysis[c].to_numpy(dtype=float) for c in columns])
            beta = np.linalg.lstsq(x, y, rcond=None)[0]
            resid = y - x @ beta
            bread = np.linalg.inv(x.T @ x)
            leverage = np.sum(x * (x @ bread), axis=1)
            meat = x.T @ np.diag((resid / np.maximum(1e-8, 1 - leverage)) ** 2) @ x
            cov = bread @ meat @ bread
            se = np.sqrt(np.diag(cov))
            df = len(y) - x.shape[1]
            for j, term in enumerate(["intercept"] + columns):
                tval = beta[j] / se[j]
                pval = 2 * student_t.sf(abs(tval), df)
                model_rows.append({"outcome": outcome, "model": model_name, "term": term, "n": len(y), "estimate": beta[j], "HC3_SE": se[j], "t": tval, "df": df, "p": pval, "CI95_low": beta[j] - student_t.ppf(0.975, df) * se[j], "CI95_high": beta[j] + student_t.ppf(0.975, df) * se[j]})
    models = pd.DataFrame(model_rows)
    models.to_csv(tables / "reviewer_continuous_models_HC3.csv", index=False)

    # Extreme-group cutoff sensitivity.
    cutoff_rows = []
    sorted_subject = subject.sort_values("endpoint_shannon")
    for n_extreme in [10, 15, 20]:
        low = sorted_subject.head(n_extreme)
        high = sorted_subject.tail(n_extreme)
        for outcome in ["baseline_shannon", "absolute_shannon_displacement", "endpoint_centroid_bray", "endpoint_binary_centroid_jaccard", "bray_auc_actual_day_normalized"]:
            test = ttest_ind(high[outcome], low[outcome], equal_var=False)
            cutoff_rows.append({"extreme_n_per_group": n_extreme, "outcome": outcome, "mean_low": low[outcome].mean(), "mean_high": high[outcome].mean(), "difference_high_minus_low": high[outcome].mean() - low[outcome].mean(), "welch_t": test.statistic, "p": test.pvalue})
    cutoff = pd.DataFrame(cutoff_rows)
    cutoff["q_BH_within_outcome"] = cutoff.groupby("outcome")["p"].transform(bh_adjust)
    cutoff.to_csv(tables / "reviewer_extreme_group_sensitivity.csv", index=False)

    # How much do genuine recovery rankings overlap the endpoint-diversity HDR/LDR labels?
    original_hdr = set(subject.index[subject["recovery_group"] == "HDR"])
    original_ldr = set(subject.index[subject["recovery_group"] == "LDR"])
    overlap_rows = []
    for metric in ["absolute_shannon_displacement", "endpoint_centroid_bray", "endpoint_binary_centroid_jaccard", "shannon_auc_actual_day_normalized", "bray_auc_actual_day_normalized"]:
        ordered = subject[metric].sort_values()
        high_recovery = set(ordered.head(15).index)
        low_recovery = set(ordered.tail(15).index)
        overlap_rows.append(
            {
                "recovery_metric_lower_is_better": metric,
                "high_recovery_overlap_original_HDR_n_of_15": len(high_recovery & original_hdr),
                "low_recovery_overlap_original_LDR_n_of_15": len(low_recovery & original_ldr),
            }
        )
    pd.DataFrame(overlap_rows).to_csv(tables / "reviewer_recovery_ranking_overlap.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.4), layout="constrained")
    outcomes_plot = [
        ("absolute_shannon_displacement", "|Endpoint − baseline| Shannon"),
        ("endpoint_centroid_bray", "Baseline–endpoint centroid Bray–Curtis"),
        ("endpoint_binary_centroid_jaccard", "Binary baseline–endpoint centroid Jaccard"),
    ]
    for ax, (outcome, ylabel) in zip(axes.flat[:3], outcomes_plot):
        colors = np.where(subject["bacteroides_plus_phocaeicola_dominant"], "#D55E00", "#0072B2")
        ax.scatter(subject["baseline_shannon"], subject[outcome], c=colors, s=30, alpha=0.8)
        coef = np.polyfit(subject["baseline_shannon"], subject[outcome], 1)
        xx = np.linspace(subject["baseline_shannon"].min(), subject["baseline_shannon"].max(), 100)
        ax.plot(xx, np.polyval(coef, xx), color="black", linestyle="--", linewidth=1)
        rho, p = safe_spearman(subject["baseline_shannon"].to_numpy(), subject[outcome].to_numpy())
        ax.set(xlabel="Baseline mean Shannon", ylabel=ylabel, title=f"Full cohort n=60: rho={rho:.2f}, p={p:.3g}")
    legend_handles = [
        mpl.lines.Line2D([], [], linestyle="none", marker="o", color="#D55E00", label="Bacteroides+Phocaeicola dominant"),
        mpl.lines.Line2D([], [], linestyle="none", marker="o", color="#0072B2", label="Other dominant genus"),
    ]
    axes[0, 0].legend(handles=legend_handles, frameon=False, fontsize=6, loc="upper left")
    panel_label(axes[0, 0], "A")
    panel_label(axes[0, 1], "B")
    panel_label(axes[1, 0], "C")

    ax = axes[1, 1]
    for outcome, color, marker in [("endpoint_centroid_bray", "#0072B2", "o"), ("endpoint_binary_centroid_jaccard", "#D55E00", "s")]:
        part = cutoff.loc[cutoff["outcome"] == outcome]
        ax.plot(part["extreme_n_per_group"], part["difference_high_minus_low"], color=color, marker=marker, label=outcome)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set(xlabel="Extreme subjects per group", ylabel="High endpoint group − low endpoint group", xticks=[10, 15, 20])
    ax.set_title("Sensitivity to arbitrary extreme-group cutoff")
    ax.legend(frameon=False, fontsize=6)
    panel_label(ax, "D")
    fig.suptitle("Reviewer-focused recovery analysis independent of absolute endpoint diversity", fontsize=12)
    export_figure(fig, figures / "reviewer_reanalysis")


def write_manifest(root: Path, output: Path) -> None:
    rows = [
        ("Figure 1B", "reproduced with audit caveat", "Shannon boxes reproduced; paired-vs-A1 tests are explicit, but the manuscript's unspecified test does not reproduce all printed P values"),
        ("Figure 1C", "descriptive reproduction", "Jaccard values reproduced; manuscript does not identify the inferential contrast"),
        ("Figure 1D", "reproduced", "top-30 named species after taxonomy aggregation"),
        ("Figure 2A-E", "reproduced with audit caveat", "group correlations and k-means boundaries disagree with the manuscript; richness scale cannot match because the table has only 2,035 rows"),
        ("Figure 2F", "method proxy", "Mann-Whitney/rank-biserial heat map; exact LEfSe preprocessing unavailable"),
        ("Figure 3", "not reproduced", "MDSINE2 parameters, count/absolute abundance input, seeds and posterior outputs absent"),
        ("Figure 4A/B/D", "reproduced", "stacked composition used for alluvial panel; supplied-table richness clearly labelled"),
        ("Figure 4C", "method proxy", "subject-mean Spearman networks replace FastSpar/SparCC bootstrap"),
        ("Figure 5 / Supplementary 4", "not reproducible", "HUMAnN/MetaCyc/KO/EC and stratified contribution tables absent"),
        ("Figure 6 / Supplementary 5", "not reproducible", "mouse qPCR, 16S and shotgun data absent"),
        ("Supplementary Figure 1", "reproduced", "panel D implemented as the plotted phase A-vs-C comparison, not the erroneous legend"),
        ("Supplementary Figure 2", "reproduced", "three-class trajectory rule explicit; Jaccard PCoA and subject-restricted PERMANOVA use fixed seed"),
        ("Supplementary Figure 3", "reproduced", "raw P<0.05 and positive rho exactly as stated; no hidden multiplicity correction"),
        ("Reviewer reanalysis", "new analysis", "continuous subject-level recovery metrics, actual-day AUC, HC3 models, Bacteroides+Phocaeicola adjustment, cutoff sensitivity"),
    ]
    pd.DataFrame(rows, columns=["item", "status", "scope_or_limitation"]).to_csv(output / "reproduction_manifest.tsv", sep="\t", index=False)
    versions = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": mpl.__version__,
        "networkx": nx.__version__,
        "scikit_learn": importlib.metadata.version("scikit-learn"),
        "Pillow": importlib.metadata.version("Pillow"),
        "random_seed": SEED,
        "input_root": str(root.resolve()),
    }
    (output / "software_versions.json").write_text(json.dumps(versions, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "outputs")
    args = parser.parse_args()
    figures = args.output / "figures"
    tables = args.output / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    set_style()
    ds = load_dataset(args.data_root)
    write_data_audit(ds, tables)
    figure1(ds, figures, tables)
    figure2(ds, figures, tables)
    supplementary_figure1(ds, figures, tables)
    supplementary_figure2(ds, figures, tables)
    supplementary_figure3(ds, figures, tables)
    figure4(ds, figures, tables)
    reviewer_reanalysis(ds, figures, tables)
    write_manifest(args.data_root, args.output)
    print(f"Completed reproduction. Outputs: {args.output.resolve()}")


if __name__ == "__main__":
    main()
