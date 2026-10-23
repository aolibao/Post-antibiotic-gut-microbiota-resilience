#!/usr/bin/env python3
"""Render a source-calibrated reconstruction of manuscript Figure 3.

The renderer consumes the audit tables made by ``summarize_posteriors.py``.
For layout development only it can read the two historical short-chain export
directories with ``--development-fallback``.  That mode is watermarked in the
manifest and must never be described as the final no-qPCR result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
from PIL import Image
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

# Required publication settings: editable vector text and a robust sans stack.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

FIGSIZE = (7.805, 11.162)
TIFF_DPI = 1000
VECTOR_RASTER_DPI = 600
TIME_ORDER = [f"A{i}" for i in range(1, 5)] + [f"B{i}" for i in range(1, 7)] + [f"C{i}" for i in range(1, 7)]
LDR_BLUE = "#0072B5"
HDR_RED = "#BC3C29"
GREY = "#B9B9B9"
GENUS_COLORS = {
    "Faecalibacterium": "#2E9E44",
    "Bacteroides": "#8E2A8A",
    "Alistipes": "#F2A51A",
}
DEV_EXPORTS = {
    "LDR": Path(
        "/home/qqt/MDSINE2_Server_Package/MDSINE2_Paper/datasets/p1/output/"
        "mdsine2/inference_fast_from_raw_tables/p1-seed0-raw/figures_p1"
    ),
    "HDR": Path(
        "/home/qqt/MDSINE2_Server_Package/MDSINE2_Paper/datasets/p4/output/"
        "mdsine2/inference_fast_raw/p4-seed1-raw/posterior_P4"
    ),
}


@dataclass
class GroupData:
    group: str
    cocluster: pd.DataFrame
    modules: pd.Series
    perturbation: pd.DataFrame
    growth: pd.Series
    zscores: pd.DataFrame
    network: pd.DataFrame
    source_mode: str


def publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.size": 6.4,
            "axes.labelsize": 7.0,
            "axes.titlesize": 7.2,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.fontsize": 5.8,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_square(path: Path) -> pd.DataFrame:
    """Read a labelled square matrix, tolerating comma files named .tsv."""
    first = path.open("r", encoding="utf-8").readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    frame = pd.read_csv(path, sep=sep, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    common = frame.index.intersection(frame.columns, sort=False)
    if len(common) < 2:
        raise ValueError(f"Not a labelled square matrix: {path}")
    frame = frame.loc[common, common].apply(pd.to_numeric, errors="coerce")
    values = frame.to_numpy(float)
    values = np.nan_to_num(values, nan=1.0)
    values = np.clip((values + values.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(values, 1.0)
    return pd.DataFrame(values, index=common, columns=common)


def _read_module_table(path: Path) -> pd.Series:
    frame = pd.read_csv(path, sep=None, engine="python")
    lower = {str(c).lower().replace(" ", "_"): c for c in frame.columns}
    taxon_col = lower.get("taxon") or lower.get("name") or frame.columns[0]
    module_col = (
        lower.get("module")
        or lower.get("cluster_assignment")
        or lower.get("cluster")
        or frame.columns[-1]
    )
    out = frame.set_index(taxon_col)[module_col]
    out.index = out.index.astype(str)
    # Preserve explicit numeric/M<number> labels so independently exported
    # network edges keep the same node identity.  Only opaque labels are
    # compacted, with a deterministic natural-sort mapping.
    labels = out.astype(str)
    parsed = labels.str.extract(r"^(?:M)?(\d+)$", expand=False)
    if parsed.notna().all():
        values = parsed.astype(int)
    else:
        unique = sorted(
            labels.unique(),
            key=lambda value: tuple(
                (0, int(part)) if part.isdigit() else (1, part)
                for part in re.split(r"(\d+)", value)
                if part
            ),
        )
        mapping = {value: i for i, value in enumerate(unique)}
        values = labels.map(mapping).astype(int)
    return pd.Series(values.to_numpy(), index=out.index, name="module", dtype=int)


def _read_effect_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    frame.index = frame.index.astype(str)
    return frame


def _old_group_exports(group: str) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]:
    root = DEV_EXPORTS[group]
    cocluster = _read_square(root / "clustering" / "coclusters.tsv")
    modules = _read_module_table(root / "clustering" / "clusterassignments.tsv")
    effect = pd.read_csv(root / "antibiotics" / "values.tsv", sep="\t", index_col=0)
    bf = pd.read_csv(root / "antibiotics" / "bayes_factors.tsv", sep="\t", index_col=0).iloc[0]
    effect.index = effect.index.astype(str)
    effect["bf"] = pd.to_numeric(bf.reindex(effect.index), errors="coerce")
    effect = effect.rename(
        columns={
            "25th percentile": "q025",
            "75th percentile": "q975",
        }
    )
    growth_frame = pd.read_csv(root / "growth" / "values.tsv", sep="\t")
    # The final token in the verbose name is the stable s__ identifier.
    growth_frame["taxon"] = growth_frame["name"].astype(str).str.extract(r"(s__\S+)$", expand=False)
    growth = growth_frame.set_index("taxon")["median"]
    return cocluster, modules, effect, growth


def _load_standard_group(group: str, tables: Path) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series] | None:
    required = {
        "cocluster": tables / f"{group}_cocluster_consensus.csv",
        "modules": tables / f"{group}_module_assignments.csv",
        "effect": tables / "perturbation_effects_combined.csv",
        "growth": tables / "growth_effects_combined.csv",
    }
    if not all(p.exists() for p in required.values()):
        return None
    cocluster = _read_square(required["cocluster"])
    modules = _read_module_table(required["modules"])
    effects = pd.read_csv(required["effect"])
    taxon_col = "taxon" if "taxon" in effects else effects.columns[0]
    effects = effects.set_index(taxon_col)
    prefix = f"{group}_"
    rename = {c: c[len(prefix):] for c in effects.columns if c.startswith(prefix)}
    effect = effects[list(rename)].rename(columns=rename)
    if "bf" not in effect:
        for candidate in ("bayes_factor", "BF"):
            if candidate in effect:
                effect["bf"] = effect[candidate]
                break
    if "q025" not in effect and "ci_low" in effect:
        effect["q025"] = effect["ci_low"]
    if "q975" not in effect and "ci_high" in effect:
        effect["q975"] = effect["ci_high"]
    growths = pd.read_csv(required["growth"])
    taxon_col = "taxon" if "taxon" in growths else growths.columns[0]
    growths = growths.set_index(taxon_col)
    candidates = [f"{group}_median", f"{group}_growth_median", "median"]
    growth_col = next((c for c in candidates if c in growths), None)
    if growth_col is None:
        raise ValueError(f"No growth median column for {group}: {required['growth']}")
    return cocluster, modules, effect, growths[growth_col]


def _terminal_species(path: str) -> str:
    return str(path).rsplit("|", 1)[-1]


def _load_raw_stage_zscores(group: str, taxa: Iterable[str], project_root: Path) -> pd.DataFrame:
    """Calculate taxon-wise stage z-scores from the supplied Bracken reads."""
    meta = pd.read_csv(project_root / "raw" / "meta.csv", dtype={"ids": str, "subject": str})
    assignments = pd.read_csv(
        project_root / "reproduction" / "outputs" / "integrated_fig124_raw" / "tables"
        / "raw_bracken_hdr_ldr_assignments.csv",
        dtype={"subject": str},
    )
    subjects = assignments.loc[assignments["recovery_group"] == group, "subject"].tolist()
    meta = meta.loc[meta["subject"].isin(subjects) & meta["canon"].isin(TIME_ORDER)].copy()
    raw = pd.read_csv(project_root / "raw" / "combine_brk.mpa", sep="\t", index_col=0)
    raw.columns = raw.columns.str.removesuffix(".brk.report")
    rows = raw.index.to_series().astype(str)
    mask = rows.str.startswith("k__Bacteria|") & rows.str.rsplit("|", n=1).str[-1].str.startswith("s__")
    species = raw.loc[mask, meta["ids"].tolist()].copy()
    species.index = [_terminal_species(x) for x in species.index]
    species = species.groupby(level=0, sort=False).sum()
    taxa = pd.Index([str(x) for x in taxa])
    species = species.reindex(taxa).fillna(0.0)
    # Match the model feature universe: proportions among the 259 retained taxa.
    rel = species.div(species.sum(axis=0).replace(0, np.nan), axis=1)
    by_stage = pd.DataFrame(
        {
            stage: rel[meta.loc[meta["canon"] == stage, "ids"]].mean(axis=1)
            for stage in TIME_ORDER
        }
    )
    center = by_stage.mean(axis=1)
    scale = by_stage.std(axis=1, ddof=0).replace(0, np.nan)
    z = by_stage.sub(center, axis=0).div(scale, axis=0).fillna(0.0)
    return z.clip(-3.0, 3.0)


def _load_prepared_stage_zscores(group: str, taxa: Iterable[str], inputs: Path) -> pd.DataFrame | None:
    group_dir = inputs / group
    reads_path, meta_path = group_dir / "reads.tsv", group_dir / "metadata.tsv"
    if not reads_path.exists() or not meta_path.exists():
        return None
    reads = pd.read_csv(reads_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t")
    sample_col = next(
        (c for c in ["sampleID", "sample", "sample_id", "ids", "name"] if c in meta),
        meta.columns[0],
    )
    stage_col = next((c for c in ["canon", "stage", "timepoint"] if c in meta), None)
    if stage_col is None:
        return None
    samples = meta[sample_col].astype(str)
    reads.index = reads.index.astype(str)
    reads.columns = reads.columns.astype(str)
    if set(samples).issubset(reads.index) and not set(samples).issubset(reads.columns):
        reads = reads.T
    if not set(samples).issubset(reads.columns):
        return None
    reads = reads.reindex(index=pd.Index([str(x) for x in taxa]), fill_value=0)
    rel = reads.div(reads.sum(axis=0).replace(0, np.nan), axis=1)
    by_stage = pd.DataFrame(
        {
            stage: rel[meta.loc[meta[stage_col] == stage, sample_col].astype(str)].mean(axis=1)
            for stage in TIME_ORDER
        }
    )
    z = by_stage.sub(by_stage.mean(axis=1), axis=0).div(by_stage.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return z.fillna(0).clip(-3, 3)


def _bh(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def _module_network(zscores: pd.DataFrame, modules: pd.Series) -> pd.DataFrame:
    common = zscores.index.intersection(modules.index, sort=False)
    profiles = zscores.loc[common].groupby(modules.loc[common]).mean()
    records: list[dict[str, float | int]] = []
    for i, source in enumerate(profiles.index):
        for target in profiles.index[i + 1 :]:
            rho, p = spearmanr(profiles.loc[source], profiles.loc[target])
            records.append({"source": int(source), "target": int(target), "rho": rho, "p": p})
    edges = pd.DataFrame(records)
    if edges.empty:
        return pd.DataFrame(columns=["source", "target", "rho", "p", "q"])
    edges["q"] = _bh(edges["p"].fillna(1).to_numpy())
    selected = edges.loc[(edges["rho"].abs() >= 0.55) & (edges["q"] < 0.10)].copy()
    if selected.empty:
        selected = edges.reindex(edges["rho"].abs().sort_values(ascending=False).index).head(min(12, len(edges))).copy()
    return selected


def _read_network(group: str, tables: Path, zscores: pd.DataFrame, modules: pd.Series) -> pd.DataFrame:
    path = tables / f"module_spearman_edges_{group}.csv"
    if path.exists():
        frame = pd.read_csv(path)
        if {"source", "target", "rho"}.issubset(frame.columns):
            if "selected_q05" in frame and frame["selected_q05"].astype(bool).any():
                return frame.loc[frame["selected_q05"].astype(bool)].copy()
            if "q" in frame and (pd.to_numeric(frame["q"], errors="coerce") < 0.05).any():
                return frame.loc[pd.to_numeric(frame["q"], errors="coerce") < 0.05].copy()
            if "selected_nominal_p05" in frame and frame["selected_nominal_p05"].astype(bool).any():
                return frame.loc[frame["selected_nominal_p05"].astype(bool)].copy()
            # An empty statistically selected network is a valid result.  Do
            # not manufacture visually convenient edges.
            return frame.iloc[0:0].copy()
    return _module_network(zscores, modules)


def load_group(
    group: str,
    tables: Path,
    inputs: Path,
    project_root: Path,
    development_fallback: bool,
) -> GroupData:
    standard = _load_standard_group(group, tables)
    if standard is None:
        if not development_fallback:
            raise FileNotFoundError(
                f"Final summary tables for {group} are incomplete. Run summarize_posteriors.py, "
                "or explicitly use --development-fallback for a non-final visual test."
            )
        cocluster, modules, perturbation, growth = _old_group_exports(group)
        mode = "development_only_historical_short_chain"
    else:
        cocluster, modules, perturbation, growth = standard
        summary_manifest = tables / "posterior_summary_manifest.json"
        tier = "unknown"
        if summary_manifest.exists():
            try:
                tier = str(json.loads(summary_manifest.read_text(encoding="utf-8"))["groups"][group]["evidence_tier"])
            except (KeyError, TypeError, json.JSONDecodeError):
                tier = "unknown"
        mode = "formal_no_qpcr_summary" if tier.startswith("formal") else f"{tier}_posterior_summary"
    taxa = cocluster.index.intersection(modules.index, sort=False)
    cocluster = cocluster.loc[taxa, taxa]
    modules = modules.reindex(taxa)
    perturbation = perturbation.reindex(taxa)
    growth = pd.to_numeric(growth.reindex(taxa), errors="coerce")
    zscores = _load_prepared_stage_zscores(group, taxa, inputs)
    if zscores is None:
        zscores = _load_raw_stage_zscores(group, taxa, project_root)
    network = _read_network(group, tables, zscores, modules)
    return GroupData(group, cocluster, modules, perturbation, growth, zscores, network, mode)


def _ordered_cocluster(data: GroupData) -> tuple[pd.DataFrame, np.ndarray]:
    dist = np.clip(1.0 - data.cocluster.to_numpy(float), 0, 1)
    np.fill_diagonal(dist, 0)
    linkage = hierarchy.linkage(squareform(dist, checks=False), method="average", optimal_ordering=True)
    order = hierarchy.leaves_list(linkage)
    return data.cocluster.iloc[order, order], linkage


def draw_cluster_panel(fig: plt.Figure, slot, data: GroupData, color: str, label: str) -> None:
    ordered, link = _ordered_cocluster(data)
    gs = slot.subgridspec(
        2,
        4,
        width_ratios=[0.052, 0.035, 0.145, 0.768],
        height_ratios=[0.18, 0.82],
        wspace=0.02,
        hspace=0.015,
    )
    ax_top = fig.add_subplot(gs[0, 3])
    ax_left = fig.add_subplot(gs[1, 2])
    ax_heat = fig.add_subplot(gs[1, 3])
    ax_cbar = fig.add_subplot(gs[1, 0])
    hierarchy.dendrogram(link, ax=ax_top, no_labels=True, color_threshold=0, above_threshold_color="#505050")
    hierarchy.dendrogram(
        link,
        ax=ax_left,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#505050",
    )
    for ax in (ax_top, ax_left):
        ax.set_axis_off()
        for collection in ax.collections:
            collection.set_linewidth(0.38)
    cmap = LinearSegmentedColormap.from_list(f"white_{data.group}", ["#FFFFFF", color])
    image = ax_heat.imshow(ordered, cmap=cmap, vmin=0, vmax=1, interpolation="nearest", rasterized=True)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("#222222")
    cb = fig.colorbar(image, cax=ax_cbar, orientation="vertical", ticks=np.linspace(0, 1, 6))
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.2, width=0.5, length=2)
    ax_top.text(-0.29, 1.12, label, transform=ax_top.transAxes, fontsize=11, va="top", ha="left")
    ax_top.text(
        0.01,
        1.08,
        f"{data.group} · {data.modules.nunique()} modules",
        transform=ax_top.transAxes,
        fontsize=6.2,
        fontweight="bold",
        color=color,
        va="top",
    )


def _pretty_taxon(taxon: str) -> str:
    text = str(taxon).removeprefix("s__").replace("_", " ")
    return text.replace("[", "").replace("]", "")


def _combined_effects(ldr: GroupData, hdr: GroupData) -> pd.DataFrame:
    taxa = ldr.perturbation.index.intersection(hdr.perturbation.index, sort=False)
    out = pd.DataFrame(index=taxa)
    for group, data in [("LDR", ldr), ("HDR", hdr)]:
        median_col = "median" if "median" in data.perturbation else "mean"
        out[f"{group}_median"] = pd.to_numeric(data.perturbation.loc[taxa, median_col], errors="coerce")
        out[f"{group}_bf"] = pd.to_numeric(data.perturbation.loc[taxa].get("bf", np.nan), errors="coerce")
    out["contrast"] = (out["LDR_median"] - out["HDR_median"]).abs()
    out["evidence"] = np.log10(1 + out[["LDR_bf", "HDR_bf"]].replace(np.inf, 1e6).max(axis=1).fillna(0))
    out["rank_score"] = out["contrast"] * (1 + 0.12 * out["evidence"])
    return out


def _source_calibrated_panel_c(effects: pd.DataFrame) -> pd.DataFrame:
    """Resolve the 25 visible source slots without inventing ambiguous IDs.

    The raster labels six Faecalibacterium and two Alistipes entries only as
    ``sp.``.  Those slots are resolved deterministically to the strongest
    negative LDR-minus-HDR contrasts within the matching genus.  Full IDs and
    the resolution rule are retained in the exported source-data table.
    """

    anchors = {
        "blautia_argi": "s__Blautia_argi",
        "blautia_obeum": "s__Blautia_obeum",
        "faecium": "s__Bacteroides_faecium",
        "blautia_hydrogenotrophica": "s__Blautia_hydrogenotrophica",
        "bacteroides_humanifaecis": "s__Bacteroides_humanifaecis",
        "bacteroides_eggerthii": "s__Bacteroides_eggerthii",
        "alistipes_dispar": "s__Alistipes_dispar",
        "faecalibacterium_prausnitzii": "s__Faecalibacterium_prausnitzii",
        "alistipes_shahii": "s__Alistipes_shahii",
        "alistipes_megaguti": "s__Alistipes_megaguti",
        "alistipes_communis": "s__Alistipes_communis",
        "faecalibacterium_taiwanense": "s__Faecalibacterium_taiwanense",
        "alistipes_indistinctus": "s__Alistipes_indistinctus",
        "alistipes_ihumii": "s__Alistipes_ihumii",
        "faecalibacterium_duncaniae": "s__Faecalibacterium_duncaniae",
        "faecalibacterium_wellingii": "s__Faecalibacterium_wellingii",
    }
    negative = effects.assign(signed=effects["LDR_median"] - effects["HDR_median"])
    faec_sp = negative.loc[
        negative.index.str.startswith("s__Faecalibacterium_sp."),
    ].sort_values(["signed", "rank_score"], ascending=[True, False]).index.tolist()[:6]
    alistipes_sp = negative.loc[
        negative.index.str.startswith("s__Alistipes_sp."),
    ].sort_values(["signed", "rank_score"], ascending=[True, False]).index.tolist()[:2]
    phoc = negative.loc[negative.index.str.startswith("s__Phocaeicola_")].sort_values(
        ["signed", "rank_score"], ascending=[True, False]
    ).index.tolist()[:1]
    slots: list[tuple[str | None, str, str]] = [
        (anchors["blautia_argi"], "Blautia argi", "exact source label"),
        (anchors["blautia_obeum"], "Blautia obeum", "exact source label"),
        (
            anchors["faecium"],
            "Bacteroides faecium",
            "source says Brachybacterium faecium, but the supplied 259-taxon universe contains Bacteroides faecium only",
        ),
        (anchors["blautia_hydrogenotrophica"], "Blautia hydrogenotrophica", "exact source taxon"),
        (anchors["bacteroides_humanifaecis"], "Bacteroides humanifaecis", "exact source taxon"),
        (anchors["bacteroides_eggerthii"], "Bacteroides eggerthii", "exact source taxon"),
        (faec_sp[0] if len(faec_sp) > 0 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (faec_sp[1] if len(faec_sp) > 1 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (phoc[0] if phoc else None, "Phocaeicola sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (faec_sp[2] if len(faec_sp) > 2 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (anchors["alistipes_dispar"], "Alistipes dispar", "source Alidtipes typo corrected"),
        (anchors["faecalibacterium_prausnitzii"], "Faecalibacterium prausnitzii", "exact source taxon"),
        (alistipes_sp[0] if len(alistipes_sp) > 0 else None, "Alistipes sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (alistipes_sp[1] if len(alistipes_sp) > 1 else None, "Alistipes sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (anchors["alistipes_shahii"], "Alistipes shahii", "source Alidtipes typo corrected"),
        (anchors["alistipes_megaguti"], "Alistipes megaguti", "source Alidtipes typo corrected"),
        (anchors["alistipes_communis"], "Alistipes communis", "source Alidtipes typo corrected"),
        (faec_sp[3] if len(faec_sp) > 3 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (anchors["faecalibacterium_taiwanense"], "Faecalibacterium taiwanense", "exact source taxon"),
        (faec_sp[4] if len(faec_sp) > 4 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
        (anchors["alistipes_indistinctus"], "Alistipes indistinctus", "source Alidtipes typo corrected"),
        (anchors["alistipes_ihumii"], "Alistipes ihumii", "source Alidtipes typo corrected"),
        (anchors["faecalibacterium_duncaniae"], "Faecalibacterium duncaniae", "exact source taxon"),
        (anchors["faecalibacterium_wellingii"], "Faecalibacterium wellingii", "exact source taxon"),
        (faec_sp[5] if len(faec_sp) > 5 else None, "Faecalibacterium sp.", "ambiguous source slot resolved by within-genus contrast rank"),
    ]
    selected: list[dict[str, object]] = []
    used: set[str] = set()
    for slot, (taxon, display, note) in enumerate(slots, start=1):
        if taxon is None or taxon not in effects.index or taxon in used:
            continue
        record = effects.loc[taxon].to_dict()
        record.update(
            {
                "taxon": taxon,
                "source_visible_slot": slot,
                "display_label": display,
                "resolution_note": note,
            }
        )
        selected.append(record)
        used.add(taxon)
    # Defensive fill keeps the renderer usable if a future fixed universe
    # lacks an anchor, while making every substitute explicit in source data.
    for taxon, row in effects.sort_values("rank_score", ascending=False).iterrows():
        if len(selected) >= 25:
            break
        if taxon in used:
            continue
        record = row.to_dict()
        record.update(
            {
                "taxon": taxon,
                "source_visible_slot": len(selected) + 1,
                "display_label": _pretty_taxon(taxon),
                "resolution_note": "fallback strongest contrast because a source anchor was absent",
            }
        )
        selected.append(record)
        used.add(taxon)
    return pd.DataFrame(selected).sort_values("source_visible_slot").set_index("taxon")


def draw_dumbbell(ax: plt.Axes, ldr: GroupData, hdr: GroupData, tables: Path) -> pd.DataFrame:
    effects = _combined_effects(ldr, hdr).dropna(subset=["LDR_median", "HDR_median"])
    top = _source_calibrated_panel_c(effects)
    y = np.arange(len(top))[::-1]
    for yi, (_, row) in zip(y, top.iterrows()):
        ax.plot([row["LDR_median"], row["HDR_median"]], [yi, yi], color="#888888", lw=0.55, zorder=1)
    ax.scatter(top["LDR_median"], y, s=11, color=LDR_BLUE, edgecolor="white", linewidth=0.25, zorder=3, label="LDR")
    ax.scatter(top["HDR_median"], y, s=11, color=HDR_RED, edgecolor="white", linewidth=0.25, zorder=3, label="HDR")
    ax.axvline(0, color="#777777", lw=0.55, ls=(0, (2, 2)), zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(top["display_label"], fontstyle="italic", fontsize=5.25)
    ax.set_xlabel("Perturbation response")
    ax.tick_params(axis="y", pad=2)
    ax.legend(loc="lower right", handletextpad=0.35, borderaxespad=0.25)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.text(-0.58, 1.02, "C", transform=ax.transAxes, fontsize=11, va="bottom")
    export = top.reset_index(names="taxon")
    export.to_csv(tables / "Fig3C_top25_perturbation_dumbbell.csv", index=False)
    return export


def _genus(taxon: str) -> str:
    name = _pretty_taxon(taxon)
    return name.split()[0] if name else ""


def _module_id(value: object) -> int:
    match = re.fullmatch(r"(?:M)?(\d+)(?:\.0)?", str(value).strip())
    if match is None:
        raise ValueError(f"Unrecognized module label: {value!r}")
    return int(match.group(1))


def draw_effect_growth(ax: plt.Axes, ldr: GroupData, hdr: GroupData, tables: Path) -> pd.DataFrame:
    taxa = ldr.perturbation.index.intersection(hdr.perturbation.index, sort=False)
    records = []
    for group, data, group_color in [("LDR", ldr, LDR_BLUE), ("HDR", hdr, HDR_RED)]:
        effect_col = "median" if "median" in data.perturbation else "mean"
        frame = pd.DataFrame(
            {
                "taxon": taxa,
                "group": group,
                "perturbation_response": pd.to_numeric(data.perturbation.loc[taxa, effect_col], errors="coerce").to_numpy(),
                "growth_rate": pd.to_numeric(data.growth.reindex(taxa), errors="coerce").to_numpy(),
            }
        )
        frame["genus"] = frame["taxon"].map(_genus)
        records.append(frame)
        ax.scatter(
            frame["perturbation_response"],
            frame["growth_rate"],
            s=8,
            color=GREY,
            alpha=0.36,
            edgecolor="none",
            rasterized=True,
            zorder=1,
        )
    all_points = pd.concat(records, ignore_index=True)
    wide_x = all_points.pivot(index="taxon", columns="group", values="perturbation_response")
    wide_y = all_points.pivot(index="taxon", columns="group", values="growth_rate")
    for taxon in wide_x.dropna().index.intersection(wide_y.dropna().index):
        genus = _genus(taxon)
        if genus not in GENUS_COLORS:
            continue
        color = GENUS_COLORS[genus]
        x0, x1 = wide_x.loc[taxon, "LDR"], wide_x.loc[taxon, "HDR"]
        y0, y1 = wide_y.loc[taxon, "LDR"], wide_y.loc[taxon, "HDR"]
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="->", color=color, lw=0.65, alpha=0.72, shrinkA=2, shrinkB=2),
            zorder=2,
        )
        ax.scatter(x0, y0, s=13, facecolor=LDR_BLUE, edgecolor=color, linewidth=0.55, zorder=4)
        ax.scatter(x1, y1, s=13, facecolor=HDR_RED, edgecolor=color, linewidth=0.55, zorder=4)
    ax.axvline(0, color="#777777", lw=0.55, ls=(0, (2, 2)), zorder=0)
    ax.set_xlabel("Perturbation response")
    ax.set_ylabel("Growth rate")
    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=LDR_BLUE, markeredgecolor="white", markersize=4.5, label="LDR"),
        Line2D([], [], marker="o", color="none", markerfacecolor=HDR_RED, markeredgecolor="white", markersize=4.5, label="HDR"),
    ] + [Line2D([], [], color=c, lw=1.4, label=g) for g, c in GENUS_COLORS.items()]
    ax.legend(handles=handles, loc="upper left", ncol=1, handlelength=1.5, labelspacing=0.25)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.text(-0.11, 1.02, "D", transform=ax.transAxes, fontsize=11, va="bottom")
    all_points.to_csv(tables / "Fig3D_growth_vs_perturbation.csv", index=False)
    return all_points


def _module_palette(n: int) -> list[str]:
    bases = [
        "#E9D5D0", "#E8E4C8", "#DCE8C8", "#CDE5D9", "#CDE4E8", "#D2DDEC",
        "#DDD4EA", "#E8D4E3", "#EFDAC8", "#D8E6D2", "#D5E5E0", "#DADBED",
    ]
    return [bases[i % len(bases)] for i in range(n)]


def _radial_order(data: GroupData) -> tuple[list[str], list[int]]:
    taxa_order: list[str] = []
    modules_order: list[int] = []
    # Preserve M0...Mn clockwise, matching the manuscript rather than sorting
    # sectors by size (which silently changes module identity between panels).
    for module in sorted(int(value) for value in data.modules.unique()):
        members = data.modules.index[data.modules == module]
        if len(members) > 2:
            corr = np.corrcoef(data.zscores.loc[members].fillna(0).to_numpy())
            dist = np.clip(1 - np.nan_to_num(corr, nan=0), 0, 2)
            np.fill_diagonal(dist, 0)
            link = hierarchy.linkage(squareform(dist, checks=False), method="average", optimal_ordering=True)
            members = members[hierarchy.leaves_list(link)]
        taxa_order.extend(members.tolist())
        modules_order.append(int(module))
    return taxa_order, modules_order


def draw_radial(ax: plt.Axes, data: GroupData, label: str, tables: Path) -> pd.DataFrame:
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()
    taxa_order, module_order = _radial_order(data)
    n_taxa = len(taxa_order)
    n_modules = len(module_order)
    gap = np.deg2rad(3.2)
    usable = 2 * np.pi - n_modules * gap
    width = usable / max(n_taxa, 1)
    module_colors = dict(zip(module_order, _module_palette(n_modules)))
    cmap = mpl.colormaps["RdBu_r"]
    norm = Normalize(-2.5, 2.5)
    angle = 0.0
    taxon_angles: dict[str, float] = {}
    module_bounds: dict[int, tuple[float, float]] = {}
    for module in module_order:
        members = [x for x in taxa_order if int(data.modules.loc[x]) == module]
        start = angle
        for taxon in members:
            center = angle + width / 2
            taxon_angles[taxon] = center
            vals = data.zscores.loc[taxon, TIME_ORDER].to_numpy(float)
            for j, value in enumerate(vals):
                bottom = 0.34 + j * (0.64 / len(TIME_ORDER))
                ax.bar(
                    center,
                    0.64 / len(TIME_ORDER) + 0.001,
                    width=width * 1.015,
                    bottom=bottom,
                    color=cmap(norm(value)),
                    edgecolor="none",
                    align="center",
                    rasterized=True,
                )
            angle += width
        end = angle
        module_bounds[module] = (start, end)
        mid = (start + end) / 2
        ax.bar(
            mid,
            0.135,
            width=max(end - start, width),
            bottom=0.19,
            color=module_colors[module],
            edgecolor="white",
            linewidth=0.45,
        )
        rotation = 90 - np.degrees(mid)
        if rotation < -90:
            rotation += 180
        if rotation > 90:
            rotation -= 180
        ax.text(mid, 0.255, f"M{module}", ha="center", va="center", fontsize=4.6, rotation=rotation, rotation_mode="anchor")
        angle += gap
    # Delicate stage rings improve readability without changing the source grammar.
    theta = np.linspace(0, 2 * np.pi, 720)
    for j in [0, 4, 10, 16]:
        radius = 0.34 + j * (0.64 / len(TIME_ORDER))
        ax.plot(theta, np.full_like(theta, radius), color="white", lw=0.33, alpha=0.9, zorder=3)

    # Main panel network: an undirected force-directed module-level Spearman
    # graph, as in the source center.  Node area is module membership and node
    # fill is mean perturbation response (blue=inhibited, red=tolerant).
    graph = nx.Graph()
    graph.add_nodes_from(module_order)
    for _, edge in data.network.iterrows():
        try:
            source, target = _module_id(edge["source"]), _module_id(edge["target"])
        except (TypeError, ValueError):
            continue
        if source not in module_order or target not in module_order:
            continue
        rho = float(edge["rho"])
        graph.add_edge(source, target, rho=rho)
    size_map = data.modules.value_counts()
    effect_col = "median" if "median" in data.perturbation else "mean"
    module_effect = pd.to_numeric(data.perturbation[effect_col], errors="coerce").groupby(data.modules).mean()
    effect_norm = mpl.colors.TwoSlopeNorm(vmin=-10.0, vcenter=0.0, vmax=5.0)
    center_ax = ax.inset_axes([0.405, 0.405, 0.19, 0.19], transform=ax.transAxes)
    center_ax.set_axis_off()
    if graph.number_of_edges():
        pos = nx.spring_layout(graph, seed=20261009, weight=None, iterations=250, k=0.78)
    else:
        pos = nx.circular_layout(graph)
    edge_colors = ["#E28E2C" if graph[u][v]["rho"] >= 0 else "#7B3294" for u, v in graph.edges]
    edge_widths = [0.35 + 0.9 * abs(float(graph[u][v]["rho"])) for u, v in graph.edges]
    nx.draw_networkx_edges(graph, pos, ax=center_ax, edge_color=edge_colors, width=edge_widths, alpha=0.68)
    node_sizes = [13 + 34 * float(size_map.loc[module]) / max(float(size_map.max()), 1) for module in graph.nodes]
    node_colors = [mpl.colormaps["RdBu_r"](effect_norm(float(module_effect.get(module, 0)))) for module in graph.nodes]
    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=center_ax,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="#4F4F4F",
        linewidths=0.45,
    )
    # Keep the bottom-row labels inside each expanded polar axes.  The former
    # negative x-offset clipped E at the physical page edge and shifted F into
    # the inter-panel gutter after the radial panels were enlarged.
    ax.text(0.0, 1.025, label, transform=ax.transAxes, fontsize=11, va="bottom", ha="left")
    ax.text(
        0.02,
        0.96,
        f"{data.group} · {n_modules} modules",
        transform=ax.transAxes,
        fontsize=6.2,
        fontweight="bold",
        color=HDR_RED if data.group == "HDR" else LDR_BLUE,
        va="top",
    )
    ax.text(
        0.98,
        0.965,
        "Rings, inner → outer: A1–A4 | B1–B6 | C1–C6\n"
        "Spokes: species · inner band: module (M) · centre: module network",
        transform=ax.transAxes,
        fontsize=4.7,
        color="#222222",
        ha="right",
        va="top",
        linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "#D0D0D0", "linewidth": 0.35, "alpha": 0.92, "pad": 1.5},
        zorder=20,
    )
    order_table = pd.DataFrame(
        {
            "taxon": taxa_order,
            "module": [int(data.modules.loc[x]) for x in taxa_order],
            "theta_radians": [taxon_angles[x] for x in taxa_order],
        }
    )
    order_table.to_csv(tables / f"Fig3{label}_{data.group}_radial_order.csv", index=False)
    data.zscores.loc[taxa_order, TIME_ORDER].to_csv(tables / f"Fig3{label}_{data.group}_stage_zscores.csv")
    data.network.to_csv(tables / f"Fig3{label}_{data.group}_module_spearman_edges_drawn.csv", index=False)
    return order_table


def render(
    project_root: Path,
    output: Path,
    development_fallback: bool = False,
) -> dict[str, object]:
    publication_style()
    figures = output / "figures"
    tables = output / "tables"
    inputs = output / "inputs"
    qa = output / "qa"
    for directory in (figures, tables, qa):
        directory.mkdir(parents=True, exist_ok=True)
    ldr = load_group("LDR", tables, inputs, project_root, development_fallback)
    hdr = load_group("HDR", tables, inputs, project_root, development_fallback)

    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    outer = fig.add_gridspec(
        3,
        2,
        left=0.065,
        right=0.985,
        bottom=0.035,
        top=0.985,
        width_ratios=[1, 1],
        height_ratios=[0.285, 0.285, 0.43],
        wspace=0.18,
        hspace=0.12,
    )
    # Caption/raster audit establishes A=HDR (red) and B=LDR (blue).
    draw_cluster_panel(fig, outer[0, 0], hdr, HDR_RED, "A")
    draw_cluster_panel(fig, outer[0, 1], ldr, LDR_BLUE, "B")
    ax_c = fig.add_subplot(outer[1, 0])
    ax_d = fig.add_subplot(outer[1, 1])
    # Panel C needs a manuscript-style label gutter.  Shifting only this axes
    # preserves the equal top-row heatmap geometry and the wide panel D field.
    cpos = ax_c.get_position()
    ax_c.set_position([cpos.x0 + 0.125, cpos.y0, cpos.width - 0.125, cpos.height])
    draw_dumbbell(ax_c, ldr, hdr, tables)
    draw_effect_growth(ax_d, ldr, hdr, tables)
    # Source panel content proves E=HDR (17 modules), F=LDR (12 modules),
    # correcting the transposed group names in the supplied caption.
    ax_e = fig.add_subplot(outer[2, 0], projection="polar")
    ax_f = fig.add_subplot(outer[2, 1], projection="polar")
    # The source uses near-full-width radial panels.  Expand only the bottom
    # row into its generous side gutters while keeping the upper grid aligned.
    ax_e.set_position([0.012, -0.015, 0.480, 0.430])
    ax_f.set_position([0.508, -0.015, 0.480, 0.430])
    draw_radial(ax_e, hdr, "E", tables)
    draw_radial(ax_f, ldr, "F", tables)

    # Explicit shared keys make the radial encoding independently readable.
    cax = fig.add_axes([0.315, 0.018, 0.175, 0.010])
    sm = mpl.cm.ScalarMappable(norm=Normalize(-2.5, 2.5), cmap="RdBu_r")
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[-2.5, 0, 2.5])
    cb.set_label("Outer heatmap rings: within-species abundance z-score", fontsize=5.4, labelpad=1)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=4.9, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)

    node_cax = fig.add_axes([0.560, 0.018, 0.175, 0.010])
    node_sm = mpl.cm.ScalarMappable(
        norm=mpl.colors.TwoSlopeNorm(vmin=-10.0, vcenter=0.0, vmax=5.0),
        cmap="RdBu_r",
    )
    node_cb = fig.colorbar(node_sm, cax=node_cax, orientation="horizontal", ticks=[-10, 0, 5])
    node_cb.set_label("Centre-node fill: mean perturbation response", fontsize=5.4, labelpad=1)
    node_cb.ax.xaxis.set_label_position("top")
    node_cb.ax.tick_params(labelsize=4.9, length=1.5, pad=1)
    node_cb.outline.set_linewidth(0.4)

    legend_handles = [
        Line2D([], [], color="#E28E2C", lw=1.2, label="positive Spearman ρ"),
        Line2D([], [], color="#7B3294", lw=1.2, label="negative Spearman ρ"),
        Line2D([], [], marker="o", color="#555555", markerfacecolor="white", lw=0,
               markersize=3.2, label="small module"),
        Line2D([], [], marker="o", color="#555555", markerfacecolor="white", lw=0,
               markersize=6.0, label="large module"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower left",
        bbox_to_anchor=(0.018, 0.003),
        ncol=2,
        fontsize=4.9,
        handlelength=1.5,
        columnspacing=0.8,
        handletextpad=0.4,
        borderaxespad=0,
        title="Centre network: edge colour = association sign; node area = module size",
        title_fontsize=5.1,
        frameon=False,
    )

    stem = figures / "Fig3_restored"
    fig.savefig(stem.with_suffix(".svg"), dpi=VECTOR_RASTER_DPI, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), dpi=VECTOR_RASTER_DPI, facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=180, facecolor="white")
    rgba_tiff = stem.with_name(stem.name + ".rgba.tmp.tif")
    fig.savefig(
        rgba_tiff,
        dpi=TIFF_DPI,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    with Image.open(rgba_tiff) as image:
        image.convert("RGB").save(
            stem.with_suffix(".tif"),
            compression="tiff_lzw",
            dpi=(TIFF_DPI, TIFF_DPI),
        )
    rgba_tiff.unlink()
    plt.close(fig)

    outputs = [stem.with_suffix(x) for x in [".svg", ".pdf", ".png", ".tif"]]
    manifest: dict[str, object] = {
        "figure": "Fig3",
        "renderer": str(Path(__file__).resolve()),
        "backend": "Python/matplotlib",
        "source_mode": {"HDR": hdr.source_mode, "LDR": ldr.source_mode},
        "development_only": not all(x.source_mode.startswith("formal") for x in (hdr, ldr)),
        "final_dimensions_inches": list(FIGSIZE),
        "tiff_dpi": TIFF_DPI,
        "vector_raster_layer_dpi": VECTOR_RASTER_DPI,
        "expected_tiff_pixels": [round(FIGSIZE[0] * TIFF_DPI), round(FIGSIZE[1] * TIFF_DPI)],
        "panel_mapping": {
            "A": "HDR co-clustering matrix",
            "B": "LDR co-clustering matrix",
            "C": "top 25 between-group perturbation-response contrasts",
            "D": "growth versus perturbation response; arrows LDR to HDR",
            "E": "HDR radial longitudinal heatmap and module Spearman network",
            "F": "LDR radial longitudinal heatmap and module Spearman network",
        },
        "caption_correction": "E=HDR and F=LDR, established from the 17-versus-12 module content in the source raster.",
        "network_semantics": {
            "main_radial": "undirected module-level Spearman association across 16 canonical stages",
            "separate_audit": "directed MDSINE2 interactions with Bayes factor >10; not substituted into the radial center",
        },
        "outputs": [
            {"path": str(p), "bytes": p.stat().st_size, "sha256": _sha256(p)} for p in outputs
        ],
    }
    (output / "render_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "fig3_mdsine2",
    )
    parser.add_argument(
        "--development-fallback",
        action="store_true",
        help="Use historical short-chain exports only to test layout; marks manifest development_only.",
    )
    args = parser.parse_args()
    manifest = render(args.project_root.resolve(), args.output.resolve(), args.development_fallback)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
