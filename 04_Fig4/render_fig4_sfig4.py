#!/usr/bin/env python3
"""Standalone restoration of manuscript Figure 4 and Supplementary Figure 4.

Figure 4 is recalculated from the supplied taxonomic abundance table and
metadata.  Its network panel can consume a completed reconstruction from
``fastspar_fig4.py``; when no such result root is supplied, the standalone
renderer retains its explicitly labelled subject-mean Spearman proxy.

Supplementary Figure 4 cannot be recalculated because the HUMAnN/MetaCyc pathway
and enzyme abundance tables are absent.  The script therefore isolates the
native embedded TIFF from the supplied supplementary DOCX without changing any
pixels, and records that it is a source-reference extraction rather than a new
analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from PIL import Image
from scipy.stats import ttest_rel

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fastspar_fig4 as fastspar_pipeline
import render_high_fidelity as hifi
import run_all as audit


TIME_ORDER = audit.TIME_ORDER
SEED = audit.SEED
LDR = hifi.LDR
HDR = hifi.HDR

GENUS_ORDER = [
    "Faecalibacterium",
    "Bacteroides",
    "Phocaeicola",
    "Blautia",
    "Segatella",
    "Agathobacter",
    "Bifidobacterium",
    "Roseburia",
    "Alistipes",
    "Mediterraneibacter",
    "Ruminococcus",
    "Parabacteroides",
    "Anaerostipes",
    "Collinsella",
    "Coprococcus",
    "Dorea",
    "Anaerobutyricum",
]
NETWORK_STAGES = ["Baseline", "B6", "C1", "C2", "C3", "End"]
FIG4_GENUS_COLORS = [
    "#8DA0CA", "#FBB463", "#E21F26", "#FACEE1", "#984F9F", "#8ED2C6",
    "#B4D66C", "#397FB9", "#BEBADA", "#FED82F", "#4EAF49", "#E4C494",
    "#81B2D3", "#BC81B8", "#CEE6C4", "#FEEC6E", "#F48073",
]
NETWORK_COLOR_GENERA = [
    "Bacteroides", "Blautia", "Faecalibacterium", "Alistipes",
    "Bifidobacterium", "Akkermansia", "Parabacteroides", "Ruminococcus",
    "Roseburia", "Lacrimispora", "Streptococcus", "Anaerostipes",
    "Butyricimonas", "Clostridium", "Coprococcus", "Eubacterium",
    "Mediterraneibacter", "Phocaeicola", "Desulfovibrio",
]
STRIP_BLUE = "#397FB9"
STRIP_RED = "#BD3C29"

# Approximate edge counts digitized from the manuscript raster.  These values
# are used only to control how many already-qualified FastSpar edges are drawn;
# they never change the caption-rule or Methods-rule analytical edge tables.
SOURCE_DIGITIZED_EDGE_TARGETS = {
    ("HDR", "Baseline"): 1350,
    ("HDR", "B6"): 2350,
    ("HDR", "C1"): 2050,
    ("HDR", "C2"): 1830,
    ("HDR", "C3"): 2090,
    ("HDR", "End"): 1220,
    ("LDR", "Baseline"): 1030,
    ("LDR", "B6"): 3100,
    ("LDR", "C1"): 2000,
    ("LDR", "C2"): 1470,
    ("LDR", "C3"): 1600,
    ("LDR", "End"): 1030,
}


def set_style() -> None:
    hifi.set_hifi_style()
    hifi.PHASE_FACE.update({"A": "#F2F9FF", "B": "#FFF2F2", "C": "#F2FFF2"})
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def export_pdf_tiff(
    fig: mpl.figure.Figure,
    stem: Path,
    dpi: int,
) -> dict[str, object]:
    """Write a fixed-page mixed PDF and opaque RGB/LZW TIFF.

    Text and line artists remain vector objects in the PDF; raster artists such
    as heatmaps are embedded at the requested output resolution.
    """
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    tif = stem.with_suffix(".tif")
    fig.savefig(pdf, format="pdf", dpi=dpi, facecolor="white", transparent=False)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor="white", transparent=False)
    buffer.seek(0)
    with Image.open(buffer) as image:
        rgb = image.convert("RGB")
        rgb.save(tif, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
        pixels = rgb.size
    inches = tuple(float(value) for value in fig.get_size_inches())
    plt.close(fig)
    return {
        "stem": stem.name,
        "width_in": inches[0],
        "height_in": inches[1],
        "dpi": dpi,
        "width_px": pixels[0],
        "height_px": pixels[1],
        "formats": ["pdf", "tif"],
        "tiff_mode": "RGB",
        "tiff_compression": "LZW",
    }


def selected_genera(ds: audit.Dataset, selected: pd.DataFrame) -> list[str]:
    present = [name for name in GENUS_ORDER if name in ds.genus.index]
    ranked = ds.genus[selected["ids"]].mean(axis=1).sort_values(ascending=False)
    for name in ranked.index:
        if name not in present:
            present.append(name)
        if len(present) == len(GENUS_ORDER):
            break
    return present


def paired_diversity_tests(selected: pd.DataFrame) -> pd.DataFrame:
    """Paired time-versus-A1 tests with BH within each group/metric family."""
    rows: list[dict[str, object]] = []
    for metric in ["shannon", "richness", "simpson"]:
        for group in ["HDR", "LDR"]:
            part = selected.loc[selected["recovery_group"] == group]
            baseline = part.loc[part["canon"] == "A1"].set_index("subject")[metric]
            family: list[dict[str, object]] = []
            for time in TIME_ORDER:
                current = part.loc[part["canon"] == time].set_index("subject")[metric]
                common = baseline.index.intersection(current.index)
                current_values = current.loc[common]
                baseline_values = baseline.loc[common]
                test = (
                    ttest_rel(current_values, baseline_values)
                    if time != "A1"
                    else None
                )
                statistic = float(test.statistic) if test is not None else np.nan
                p = float(test.pvalue) if test is not None else np.nan
                family.append(
                    {
                        "metric": metric,
                        "group": group,
                        "canon": time,
                        "n_pairs": int(len(common)),
                        "baseline_A1_mean": float(baseline_values.mean()),
                        "current_mean": float(current_values.mean()),
                        "paired_mean_difference": float((current_values - baseline_values).mean()),
                        "t_statistic": statistic,
                        "p_paired_vs_A1_raw": p,
                    }
                )
            q_values = audit.bh_adjust(row["p_paired_vs_A1_raw"] for row in family)
            for row, q_value in zip(family, q_values):
                row["q_BH_within_metric_group"] = q_value
                row["symbol"] = hifi.p_symbol(float(q_value))
                rows.append(row)
    return pd.DataFrame(rows)


def draw_diversity_heatmap(
    fig: mpl.figure.Figure,
    ax: mpl.axes.Axes,
    selected: pd.DataFrame,
    tests: pd.DataFrame,
    metric: str,
    title: str,
    cmap_name: str,
    show_time_labels: bool,
) -> None:
    matrix = (
        selected.groupby(["canon", "recovery_group"], observed=True)[metric]
        .mean()
        .unstack()
        .reindex(TIME_ORDER)[["HDR", "LDR"]]
    )
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap=cmap_name,
        vmin=float(matrix.min().min()),
        vmax=float(matrix.max().max()),
        interpolation="nearest",
    )
    ax.set_xticks([0, 1], ["HDR", "LDR"])
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks(range(16), TIME_ORDER if show_time_labels else [""] * 16)
    ax.tick_params(axis="y", length=0)
    if show_time_labels:
        hifi.color_time_ticks(ax, "y")
    ax.set_xlabel(title, labelpad=4)

    subset = tests.loc[tests["metric"] == metric]
    lookup = {
        (row.canon, row.group): row.symbol
        for row in subset.itertuples(index=False)
    }
    for row_index, time in enumerate(TIME_ORDER):
        for column_index, group in enumerate(["HDR", "LDR"]):
            symbol = lookup[(time, group)]
            if not isinstance(symbol, str) or not symbol:
                continue
            rgba = image.cmap(image.norm(matrix.loc[time, group]))
            luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            ax.text(
                column_index,
                row_index,
                symbol,
                ha="center",
                va="center",
                fontsize=5.5,
                fontweight="bold",
                color="white" if luminance < 0.55 else "#222222",
            )
    for spine in ax.spines.values():
        spine.set_visible(False)
    color_axis = ax.inset_axes([0.03, -0.085, 0.94, 0.028])
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
    colorbar.ax.tick_params(labelsize=5.5, length=2, pad=1)
    for label in colorbar.ax.get_xticklabels():
        label.set_rotation(90)
        label.set_va("top")
    colorbar.outline.set_linewidth(0.35)


def strict_network_proxy(
    ds: audit.Dataset,
    group: str,
    stage: str,
) -> tuple[nx.Graph, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Spearman proxy using the strong-correlation rule printed in the caption."""
    features = hifi.select_network_features(ds, group, stage, n_features=275)
    _, rho, p_values = hifi.network_stage_data(ds, group, stage, features)
    upper = np.triu_indices(len(features), 1)
    adjusted = audit.bh_adjust(p_values[upper])
    q_values = np.full_like(p_values, np.nan, dtype=float)
    q_values[upper] = adjusted
    q_values[(upper[1], upper[0])] = adjusted
    graph = nx.Graph()
    graph.add_nodes_from(features)
    upper_i, upper_j = np.where(
        np.triu((np.abs(rho) > 0.6) & (p_values < 0.05), 1)
    )
    for first, second in zip(upper_i, upper_j):
        graph.add_edge(
            features[first],
            features[second],
            weight=float(abs(rho[first, second])),
            sign=float(np.sign(rho[first, second])),
        )
    graph.remove_nodes_from(list(nx.isolates(graph)))
    return graph, rho, p_values, q_values, features


def validate_fastspar_provenance(
    ds: audit.Dataset,
    results: fastspar_pipeline.CompletedFastSparResults,
) -> None:
    """Confirm that rendered networks were fitted from this dataset's sources."""
    expected_sources = {
        "taxonomic_source_sha256": ds.root / "raw" / "combine_brk.mpa",
        "metadata_source_sha256": ds.root / "raw" / "meta.csv",
    }
    for manifest_field, source_path in expected_sources.items():
        if not source_path.is_file():
            raise FileNotFoundError(
                f"cannot verify FastSpar provenance; source is missing: {source_path}"
            )
        observed = fastspar_pipeline.sha256_file(source_path)
        recorded = results.manifest.get(manifest_field)
        if recorded != observed:
            raise ValueError(
                f"FastSpar provenance mismatch for {source_path}: "
                f"run manifest records {recorded!r}, current SHA-256 is {observed}"
            )
    assignment_hash = results.manifest.get("group_assignment_table_sha256")
    if assignment_hash:
        assignment_path = results.root.parent / "group_assignments.tsv"
        if not assignment_path.is_file():
            raise FileNotFoundError(
                "cannot verify FastSpar grouping; assignment table is missing: "
                f"{assignment_path}"
            )
        observed_hash = fastspar_pipeline.sha256_file(assignment_path)
        if observed_hash != assignment_hash:
            raise ValueError("FastSpar group-assignment table does not match its run manifest")
        assignments = pd.read_csv(assignment_path, sep="\t", dtype={"subject": str})
        expected = (
            ds.group_table.reset_index()[["subject", "recovery_group"]]
            .astype({"subject": str, "recovery_group": str})
            .sort_values("subject", kind="mergesort")
            .reset_index(drop=True)
        )
        observed = (
            assignments[["subject", "recovery_group"]]
            .astype({"subject": str, "recovery_group": str})
            .sort_values("subject", kind="mergesort")
            .reset_index(drop=True)
        )
        if not expected.equals(observed):
            raise ValueError(
                "FastSpar grouping does not match the recovery groups injected into Figure 4"
            )
    else:
        # Backward-compatible validation for the original D0 result bundle.
        source_path = ds.root / "species.csv"
        recorded = results.manifest.get("grouping_source_sha256")
        if not source_path.is_file() or fastspar_pipeline.sha256_file(source_path) != recorded:
            raise ValueError("Legacy FastSpar grouping provenance does not match species.csv")


def reconstructed_fastspar_network(
    results: fastspar_pipeline.CompletedFastSparResults,
    group: str,
    stage: str,
) -> tuple[nx.Graph, pd.DataFrame, pd.Series]:
    """Build one full eligible graph from the validated caption-rule edge list."""
    selected = results.caption_edges.loc[
        results.caption_edges["group"].eq(group)
        & results.caption_edges["stage"].eq(stage)
    ].copy()
    graph = nx.Graph()
    for edge in selected.itertuples(index=False):
        graph.add_edge(
            str(edge.genus1),
            str(edge.genus2),
            weight=float(edge.abs_rho),
            sign=float(np.sign(edge.rho)),
        )
    summary = (
        results.summary.loc[
            results.summary["group"].eq(group)
            & results.summary["stage"].eq(stage)
        ]
        .iloc[0]
        .copy()
    )
    if graph.number_of_nodes() != int(summary["raw_p_strong_nodes"]):
        raise ValueError(f"{group}/{stage}: rendered FastSpar node count changed after loading")
    if graph.number_of_edges() != int(summary["raw_p_strong_edges"]):
        raise ValueError(f"{group}/{stage}: rendered FastSpar edge count changed after loading")
    return graph, selected, summary


def _undirected_edge_key(first: object, second: object) -> tuple[str, str]:
    """Return a stable key for one undirected genus pair."""
    names = sorted((str(first), str(second)))
    return names[0], names[1]


def source_density_display_network(
    eligible_graph: nx.Graph,
    eligible_edges: pd.DataFrame,
    group: str,
    stage: str,
) -> tuple[nx.Graph, pd.DataFrame, dict[str, object]]:
    """Select a deterministic, source-density subset of eligible edges.

    All returned edges already satisfy the FastSpar caption rule.  A maximum-
    spanning forest is retained first so that every eligible node remains in
    the displayed graph, then remaining pairs are filled in descending
    ``abs_rho`` order (raw P and genus names break ties).  The target is a
    visual calibration digitized from the source raster, not a second
    statistical threshold.  If fewer caption-qualified edges exist than the
    source target, all eligible edges are retained.
    """
    panel = (str(group), str(stage))
    if panel not in SOURCE_DIGITIZED_EDGE_TARGETS:
        raise KeyError(f"no source-density target for {group}/{stage}")
    required = {"genus1", "genus2", "rho", "abs_rho", "p_raw"}
    missing = required.difference(eligible_edges.columns)
    if missing:
        raise ValueError(
            f"{group}/{stage}: eligible edge table lacks columns {sorted(missing)}"
        )

    ranked = eligible_edges.copy()
    ranked["_edge_key"] = [
        _undirected_edge_key(first, second)
        for first, second in zip(ranked["genus1"], ranked["genus2"])
    ]
    if ranked["_edge_key"].duplicated().any():
        raise ValueError(f"{group}/{stage}: duplicate undirected eligible edges")
    ranked["_edge_lo"] = ranked["_edge_key"].map(lambda pair: pair[0])
    ranked["_edge_hi"] = ranked["_edge_key"].map(lambda pair: pair[1])
    ranked = ranked.sort_values(
        ["abs_rho", "p_raw", "_edge_lo", "_edge_hi"],
        ascending=[False, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    # Rebuild in the same deterministic rank order before running Kruskal, so
    # equal-weight choices are reproducible across repeated renders.
    ranked_graph = nx.Graph()
    ranked_graph.add_nodes_from(sorted(str(node) for node in eligible_graph.nodes()))
    for edge in ranked.itertuples(index=False):
        ranked_graph.add_edge(
            str(edge.genus1),
            str(edge.genus2),
            weight=float(edge.abs_rho),
            sign=float(np.sign(edge.rho)),
        )
    forest = nx.maximum_spanning_tree(ranked_graph, weight="weight", algorithm="kruskal")
    forest_keys = {
        _undirected_edge_key(first, second) for first, second in forest.edges()
    }

    requested_target = int(SOURCE_DIGITIZED_EDGE_TARGETS[panel])
    capped_target = min(requested_target, len(ranked))
    # In the supplied full run every target exceeds the forest size.  Keep the
    # guard so future runs still preserve nodes rather than silently dropping
    # the forest to satisfy an unusually small visual target.
    effective_target = min(len(ranked), max(capped_target, len(forest_keys)))
    selected_keys = set(forest_keys)
    if len(selected_keys) < effective_target:
        for key in ranked["_edge_key"]:
            selected_keys.add(key)
            if len(selected_keys) == effective_target:
                break

    displayed = ranked.loc[ranked["_edge_key"].isin(selected_keys)].copy()
    displayed["selection_role"] = np.where(
        displayed["_edge_key"].isin(forest_keys),
        "maximum_spanning_forest",
        "strongest_abs_rho_fill",
    )
    displayed.insert(2, "display_rank", np.arange(1, len(displayed) + 1, dtype=int))
    displayed["source_digitized_edge_target"] = requested_target
    displayed["target_capped_by_eligible_edges"] = requested_target > len(ranked)
    displayed["eligible_edges_caption_rule"] = len(ranked)
    displayed = displayed.drop(columns=["_edge_key", "_edge_lo", "_edge_hi"])

    display_graph = nx.Graph()
    display_graph.add_nodes_from(ranked_graph.nodes())
    for edge in displayed.itertuples(index=False):
        display_graph.add_edge(
            str(edge.genus1),
            str(edge.genus2),
            weight=float(edge.abs_rho),
            sign=float(np.sign(edge.rho)),
        )
    if display_graph.number_of_nodes() != eligible_graph.number_of_nodes():
        raise ValueError(f"{group}/{stage}: source-density display lost eligible nodes")
    if display_graph.number_of_edges() != effective_target:
        raise ValueError(f"{group}/{stage}: source-density display edge count is unstable")

    metadata: dict[str, object] = {
        "source_digitized_edge_target": requested_target,
        "target_capped_by_eligible_edges": requested_target > len(ranked),
        "eligible_nodes_caption_rule": eligible_graph.number_of_nodes(),
        "eligible_edges_caption_rule": eligible_graph.number_of_edges(),
        "maximum_spanning_forest_edges": len(forest_keys),
        "display_nodes_source_density": display_graph.number_of_nodes(),
        "display_edges_source_density": display_graph.number_of_edges(),
        "selection_rule": (
            "maximum-spanning forest by abs_rho, then strongest abs_rho; "
            "raw P and canonical genus names break ties"
        ),
    }
    return display_graph, displayed, metadata


def draw_network_restored(
    ax: mpl.axes.Axes,
    graph: nx.Graph,
    palette: dict[str, str],
    seed: int,
) -> None:
    """Render a network with the hub-emphasizing grammar visible in the source."""
    ax.set_axis_off()
    if graph.number_of_nodes() == 0:
        return
    layout_graph = graph.copy()
    components = sorted(nx.connected_components(layout_graph), key=len, reverse=True)
    if len(components) > 1:
        anchor = max(components[0], key=graph.degree)
        for component in components[1:]:
            representative = max(component, key=graph.degree)
            layout_graph.add_edge(anchor, representative, weight=0.02)
    position = nx.spring_layout(
        layout_graph,
        seed=seed,
        iterations=60,
        weight="weight",
        k=0.65 / math.sqrt(max(graph.number_of_nodes(), 1)),
    )
    coordinate_matrix = np.asarray([position[node] for node in graph.nodes()], dtype=float)
    coordinate_matrix -= np.median(coordinate_matrix, axis=0, keepdims=True)
    robust_scale = max(float(np.percentile(np.abs(coordinate_matrix), 95)), 1e-9)
    coordinate_matrix = np.clip(coordinate_matrix / robust_scale, -1.2, 1.2)
    position = {
        node: coordinate_matrix[index]
        for index, node in enumerate(graph.nodes())
    }
    edges = list(graph.edges())
    edge_colors = [
        "#8F8F8F" if graph.edges[edge]["sign"] > 0 else "#D0A0A0"
        for edge in edges
    ]
    nx.draw_networkx_edges(
        graph,
        position,
        ax=ax,
        edge_color=edge_colors,
        width=0.28,
        alpha=0.27,
    )
    degrees = dict(graph.degree())
    degree_values = np.asarray([degrees[node] for node in graph.nodes()], dtype=float)
    denominator = max(float(degree_values.max()), 1.0)
    sizes = 1.8 + 30.0 * np.power(degree_values / denominator, 1.05)
    # In the source, the named genera act as saturated visual anchors among
    # many small pale "Other" nodes.  Give those named nodes a readable minimum
    # area even when their local proxy-network degree is modest.
    sizes = np.maximum(
        sizes,
        np.asarray([18.0 if node in palette else 0.0 for node in graph.nodes()]),
    )
    colors = [palette.get(node, "#90D5E4") for node in graph.nodes()]
    nx.draw_networkx_nodes(
        graph,
        position,
        ax=ax,
        node_size=sizes,
        node_color=colors,
        edgecolors="none",
        alpha=0.94,
    )
    ax.set_xlim(-1.28, 1.28)
    ax.set_ylim(-1.12, 1.12)
    ax.set_aspect("auto")


def draw_correlation_panel(
    ax: mpl.axes.Axes,
    ds: audit.Dataset,
    group: str,
    metric: str,
    title: str,
    show_y: bool,
) -> list[dict[str, object]]:
    """Draw one Fig4D panel and return every displayed correlation."""
    selected = ds.metrics.loc[ds.metrics["recovery_group"] == group]
    hifi.phase_spans(ax)
    rows: list[dict[str, object]] = []
    for genus, color in [("Bacteroides", LDR), ("Faecalibacterium", HDR)]:
        coefficients: list[float] = []
        p_values: list[float] = []
        for time in TIME_ORDER:
            metadata = selected.loc[selected["canon"] == time].sort_values("subject")
            abundance = (
                ds.genus.loc[genus, metadata["ids"]].to_numpy(dtype=float)
                if genus in ds.genus.index
                else np.zeros(len(metadata), dtype=float)
            )
            coefficient, p_value = audit.safe_spearman(
                abundance, metadata[metric].to_numpy(dtype=float)
            )
            coefficients.append(coefficient)
            p_values.append(p_value)
            rows.append(
                {
                    "group": group,
                    "metric": metric,
                    "canon": time,
                    "genus": genus,
                    "n": int(len(metadata)),
                    "rho": coefficient,
                    "p_raw": p_value,
                    "symbol_raw": hifi.p_symbol(float(p_value)),
                }
            )
        ax.plot(
            range(16),
            coefficients,
            color=color,
            marker="o",
            markersize=2.45,
            linewidth=0.72,
            label=genus,
        )
        for index, (coefficient, p_value) in enumerate(zip(coefficients, p_values)):
            symbol = hifi.p_symbol(float(p_value))
            if symbol:
                ax.text(
                    index,
                    min(0.98, coefficient + 0.13),
                    symbol,
                    color=color,
                    ha="center",
                    va="bottom",
                    fontsize=5.4,
                    fontweight="bold",
                )
    ax.axhline(0, color="#888888", linestyle=(0, (2, 2)), linewidth=0.55)
    ax.set_xlim(-0.5, 15.5)
    ax.set_ylim(-0.78, 1.08)
    ax.set_xticks(range(16), TIME_ORDER, rotation=48, ha="right")
    ax.tick_params(axis="x", labelsize=6.3, pad=1)
    hifi.color_time_ticks(ax)
    if show_y:
        ax.set_ylabel("Spearman Correlation")
        ax.yaxis.set_label_coords(-0.11, 0.5)
    else:
        ax.set_yticklabels([])
    hifi.clean_axes(ax)
    strip_color = STRIP_BLUE if group == "LDR" else STRIP_RED
    ax.add_patch(
        Rectangle(
            (0, 1.015),
            1,
            0.13,
            transform=ax.transAxes,
            clip_on=False,
            facecolor=strip_color,
            edgecolor="none",
        )
    )
    ax.text(
        0.5,
        1.08,
        title,
        transform=ax.transAxes,
        color="white",
        ha="center",
        va="center",
        fontsize=8.3,
        fontweight="bold",
    )
    return rows


def draw_figure4(
    ds: audit.Dataset,
    figures: Path,
    tables: Path,
    fastspar_results: Path | fastspar_pipeline.CompletedFastSparResults | None = None,
    *,
    group_definition: str | None = None,
    summary_uses_full_eligible_network: bool = False,
) -> dict[str, object]:
    """Render Figure 4, optionally replacing panel C with validated FastSpar output."""
    if isinstance(fastspar_results, fastspar_pipeline.CompletedFastSparResults):
        completed_fastspar = fastspar_results
    elif fastspar_results is not None:
        completed_fastspar = fastspar_pipeline.load_completed_results(fastspar_results)
    else:
        completed_fastspar = None
    if completed_fastspar is not None:
        validate_fastspar_provenance(ds, completed_fastspar)
    raw_bracken = "grouping_endpoint_shannon" in ds.group_table.columns
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])]
    fig = plt.figure(figsize=(7.871, 11.325))

    # A — manuscript-style longitudinal alluvia.
    axis_a_ldr = fig.add_axes([0.170, 0.892, 0.815, 0.078])
    axis_a_hdr = fig.add_axes([0.170, 0.772, 0.815, 0.077])
    genera = selected_genera(ds, selected)
    legend_order = [*genera, "Others"]
    stack_order = ["Others", *genera[::-1]]
    palette = {
        genus: FIG4_GENUS_COLORS[index % len(FIG4_GENUS_COLORS)]
        for index, genus in enumerate(genera)
    }
    palette["Others"] = "#D3D3D3"
    for axis, group in [(axis_a_ldr, "LDR"), (axis_a_hdr, "HDR")]:
        group_meta = selected.loc[selected["recovery_group"] == group]
        composition = pd.DataFrame(
            {
                time: ds.genus[group_meta.loc[group_meta["canon"] == time, "ids"]]
                .mean(axis=1)
                for time in TIME_ORDER
            }
        )
        plotted = composition.reindex(genera).fillna(0.0)
        plotted.loc["Others"] = (100.0 - plotted.sum(axis=0)).clip(lower=0.0)
        plotted.to_csv(tables / f"Fig4A_alluvial_{group}.csv")
        hifi.alluvial_panel(axis, plotted, stack_order, palette, group)
    hifi.panel_label(axis_a_ldr, "A", x=-0.19, y=1.04)
    handles = [
        Patch(facecolor=palette[name], edgecolor="none", label=name)
        for name in legend_order
    ]
    legend = axis_a_ldr.legend(
        handles=handles,
        title="Genus",
        loc="upper right",
        bbox_to_anchor=(0.0, 1.22),
        frameon=False,
        fontsize=6.25,
        title_fontsize=7.5,
        handlelength=1.4,
        handleheight=1.1,
        labelspacing=0.24,
        borderaxespad=0,
    )
    for text in legend.get_texts():
        if text.get_text() != "Others":
            text.set_fontstyle("italic")

    # B — independently scaled diversity heat maps, with caption-consistent BH.
    diversity_axes = [
        fig.add_axes(bounds)
        for bounds in [
            [0.030, 0.411, 0.059, 0.321],
            [0.141, 0.411, 0.059, 0.321],
            [0.249, 0.411, 0.057, 0.321],
        ]
    ]
    diversity_tests = paired_diversity_tests(selected)
    diversity_tests.to_csv(tables / "Fig4B_paired_tests_BH.csv", index=False)
    diversity_specs = [
        (
            "shannon",
            "Shannon",
            LinearSegmentedColormap.from_list("fig4_shannon", ["#F8F5F1", "#6E1F1B"]),
        ),
        (
            "richness",
            "Richness" if raw_bracken else "Richness†",
            LinearSegmentedColormap.from_list("fig4_richness", ["#FFFCF8", "#D69E3F"]),
        ),
        (
            "simpson",
            "Simpson",
            LinearSegmentedColormap.from_list("fig4_simpson", ["#F3F9F7", "#3A6240"]),
        ),
    ]
    for index, (axis, spec) in enumerate(zip(diversity_axes, diversity_specs)):
        draw_diversity_heatmap(
            fig,
            axis,
            selected,
            diversity_tests,
            metric=spec[0],
            title=spec[1],
            cmap_name=spec[2],
            show_time_labels=index == 0,
        )
    fig.text(0.006, 0.748, "B", fontsize=15, ha="left", va="bottom")
    fig.text(
        0.043,
        0.352,
        (
            "Bacterial Bracken species leaves; source-figure branch."
            if raw_bracken
            else "† Supplied leaf-taxon count; original richness scale unavailable."
        ),
        fontsize=5.4,
        color="#555555",
        ha="left",
        va="bottom",
    )

    # C — 2 x 6 strong-correlation networks and summary trajectories.
    network_width = 0.096
    network_gap = (0.985 - 0.351 - 6 * network_width) / 5
    network_lefts = [
        0.351 + index * (network_width + network_gap) for index in range(6)
    ]
    network_axes: list[mpl.axes.Axes] = []
    for bottom in [0.645, 0.535]:
        for left in network_lefts:
            network_axes.append(fig.add_axes([left, bottom, network_width, 0.075]))
    tab20 = mpl.colormaps["tab20"].colors
    network_palette = {
        genus: mpl.colors.to_hex(tab20[index % len(tab20)])
        for index, genus in enumerate(NETWORK_COLOR_GENERA)
    }
    summary_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    display_edge_frames: list[pd.DataFrame] = []
    for row_index, group in enumerate(["HDR", "LDR"]):
        for stage_index, stage in enumerate(NETWORK_STAGES):
            axis = network_axes[row_index * 6 + stage_index]
            if completed_fastspar is not None:
                eligible_graph, panel_edges, panel_summary = reconstructed_fastspar_network(
                    completed_fastspar, group, stage
                )
                graph, display_edges, display_metadata = source_density_display_network(
                    eligible_graph, panel_edges, group, stage
                )
                display_edge_frames.append(display_edges)
            else:
                graph, rho, p_values, q_values, features = strict_network_proxy(
                    ds, group, stage
                )
            draw_network_restored(
                axis,
                graph,
                network_palette,
                seed=SEED + row_index * 100 + stage_index,
            )
            stage_color = hifi.TIME_TEXT["A"] if stage == "Baseline" else hifi.TIME_TEXT["C"]
            axis.text(
                0.5,
                -0.06,
                stage,
                transform=axis.transAxes,
                ha="center",
                va="top",
                fontsize=7.2,
                color=stage_color,
            )
            if completed_fastspar is not None:
                summary_rows.append(
                    {
                        "group": group,
                        "stage": stage,
                        "candidate_features": int(panel_summary["candidate_features"]),
                        "eligible_nodes_caption_rule": int(
                            display_metadata["eligible_nodes_caption_rule"]
                        ),
                        "eligible_edges_caption_rule": int(
                            display_metadata["eligible_edges_caption_rule"]
                        ),
                        "nodes": (
                            eligible_graph.number_of_nodes()
                            if summary_uses_full_eligible_network
                            else graph.number_of_nodes()
                        ),
                        "edges": (
                            eligible_graph.number_of_edges()
                            if summary_uses_full_eligible_network
                            else graph.number_of_edges()
                        ),
                        "display_nodes_source_density": int(
                            display_metadata["display_nodes_source_density"]
                        ),
                        "display_edges_source_density": int(
                            display_metadata["display_edges_source_density"]
                        ),
                        "source_digitized_edge_target": int(
                            display_metadata["source_digitized_edge_target"]
                        ),
                        "target_capped_by_eligible_edges": bool(
                            display_metadata["target_capped_by_eligible_edges"]
                        ),
                        "maximum_spanning_forest_edges": int(
                            display_metadata["maximum_spanning_forest_edges"]
                        ),
                        "eligible_edges_also_q_BH_lt_alpha": int(
                            (
                                panel_edges["q_BH_within_network"]
                                < float(completed_fastspar.manifest["alpha_strict_lt"])
                            ).sum()
                        ),
                        "display_edges_also_q_BH_lt_alpha": int(
                            (
                                display_edges["q_BH_within_network"]
                                < float(completed_fastspar.manifest["alpha_strict_lt"])
                            ).sum()
                        ),
                        "alpha_strict_lt": float(
                            completed_fastspar.manifest["alpha_strict_lt"]
                        ),
                        "BH_q_strong_nodes": int(panel_summary["BH_q_strong_nodes"]),
                        "BH_q_strong_edges": int(panel_summary["BH_q_strong_edges"]),
                        "bootstraps": int(panel_summary["bootstraps"]),
                        "seed": int(panel_summary["seed"]),
                        "method": (
                            (
                                "global fixed genus set; pooled original sample columns; "
                                "FastSpar caption-qualified edges from exact bootstrap raw P; "
                                "small-network drawing uses a display-only source-density "
                                "subset retaining the maximum-spanning forest then filling "
                                "strongest abs_rho edges; node/edge trajectory uses the full "
                                "caption-qualified eligible network"
                            )
                            if summary_uses_full_eligible_network
                            else (
                                "global fixed genus set; pooled original sample columns; "
                                "FastSpar caption-qualified edges from exact bootstrap raw P; "
                                "display-only source-density subset retains the maximum-spanning "
                                "forest then fills strongest abs_rho edges"
                            )
                        ),
                    }
                )
            else:
                feature_position = {name: index for index, name in enumerate(features)}
                displayed_q = []
                for first, second in graph.edges():
                    first_index = feature_position[first]
                    second_index = feature_position[second]
                    displayed_q.append(q_values[first_index, second_index])
                summary_rows.append(
                    {
                        "group": group,
                        "stage": stage,
                        "candidate_features": len(features),
                        "nodes": graph.number_of_nodes(),
                        "edges": graph.number_of_edges(),
                        "edges_also_q_BH_lt_0_05": int(
                            np.sum(np.asarray(displayed_q) < 0.05)
                        ),
                        "method": "stage-specific top-275 nonconstant genera; subject-mean Spearman; |rho|>0.6; asymptotic raw P<0.05",
                    }
                )
                for first, second in graph.edges():
                    first_index = feature_position[first]
                    second_index = feature_position[second]
                    edge_rows.append(
                        {
                            "group": group,
                            "stage": stage,
                            "genus1": first,
                            "genus2": second,
                            "rho": rho[first_index, second_index],
                            "p_raw": p_values[first_index, second_index],
                            "q_BH_all_candidate_pairs_sensitivity": q_values[
                                first_index, second_index
                            ],
                        }
                    )
    network_summary = pd.DataFrame(summary_rows)
    if completed_fastspar is not None:
        # Do not leave obsolete proxy tables in an output directory that is
        # being upgraded in place to the FastSpar-backed integrated bundle.
        for stale_name in [
            "Fig4C_network_proxy_summary.csv",
            "Fig4C_network_proxy_edges.csv",
        ]:
            (tables / stale_name).unlink(missing_ok=True)
        exported_fastspar_summary = completed_fastspar.summary.copy()
        display_summary = network_summary[
            [
                "group",
                "stage",
                "nodes",
                "edges",
                "eligible_nodes_caption_rule",
                "eligible_edges_caption_rule",
                "source_digitized_edge_target",
                "target_capped_by_eligible_edges",
                "maximum_spanning_forest_edges",
                "display_nodes_source_density",
                "display_edges_source_density",
                "display_edges_also_q_BH_lt_alpha",
            ]
        ].rename(
            columns={
                "nodes": "trajectory_nodes",
                "edges": "trajectory_edges",
            }
        )
        display_summary["trajectory_network_basis"] = (
            "full_caption_qualified_eligible_network"
            if summary_uses_full_eligible_network
            else "source_density_display_subset"
        )
        exported_fastspar_summary = exported_fastspar_summary.merge(
            display_summary, on=["group", "stage"], how="left", validate="one_to_one"
        )
        exported_fastspar_summary.to_csv(
            tables / "Fig4C_fastspar_summary.csv", index=False
        )
        completed_fastspar.caption_edges.to_csv(
            tables / "Fig4C_fastspar_edges_caption_raw_p.csv", index=False
        )
        completed_fastspar.methods_edges.to_csv(
            tables / "Fig4C_fastspar_edges_methods_BH_q.csv", index=False
        )
        pd.concat(display_edge_frames, ignore_index=True).to_csv(
            tables / "Fig4C_fastspar_edges_display_source_density.csv", index=False
        )
    else:
        network_summary.to_csv(tables / "Fig4C_network_proxy_summary.csv", index=False)
        pd.DataFrame(edge_rows).to_csv(
            tables / "Fig4C_network_proxy_edges.csv", index=False
        )

    network_span = (0.985 - 0.351) / network_width
    for axis, group in [(network_axes[0], "HDR"), (network_axes[6], "LDR")]:
        axis.add_patch(
            Rectangle(
                (0, 1.05),
                network_span,
                0.16,
                transform=axis.transAxes,
                clip_on=False,
                facecolor=STRIP_RED if group == "HDR" else STRIP_BLUE,
                edgecolor="none",
            )
        )
        axis.text(
            network_span / 2,
            1.13,
            group,
            transform=axis.transAxes,
            color="white",
            ha="center",
            va="center",
            fontsize=7.8,
            fontweight="bold",
        )
    fig.text(0.327, 0.733, "C", fontsize=15, ha="left", va="bottom")

    legend_axis = fig.add_axes([0.351, 0.370, 0.125, 0.145])
    legend_axis.set_axis_off()
    legend_axis.text(0.02, 0.98, "Genus", va="top", fontsize=7.5)
    network_handles = [
        mpl.lines.Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            color=network_palette[genus],
            markersize=4,
            label=genus,
        )
        for genus in NETWORK_COLOR_GENERA
    ]
    network_handles.append(
        mpl.lines.Line2D(
            [], [], marker="o", linestyle="none", color="#90D5E4", markersize=4, label="Other"
        )
    )
    network_legend = legend_axis.legend(
        handles=network_handles,
        loc="upper left",
        bbox_to_anchor=(0, 0.88),
        frameon=False,
        fontsize=4.9,
        handletextpad=0.25,
        labelspacing=0.08,
        borderpad=0,
    )
    for text in network_legend.get_texts():
        if text.get_text() != "Other":
            text.set_fontstyle("italic")

    nodes_axis = fig.add_axes([0.500, 0.385, 0.220, 0.100])
    edges_axis = fig.add_axes([0.770, 0.385, 0.220, 0.100])
    for axis, field, ylabel in [
        (nodes_axis, "nodes", "Number of Nodes"),
        (edges_axis, "edges", "Number of Edges"),
    ]:
        for group, color in [("LDR", LDR), ("HDR", HDR)]:
            part = (
                network_summary.loc[network_summary["group"] == group]
                .set_index("stage")
                .reindex(NETWORK_STAGES)
            )
            axis.plot(
                range(6),
                part[field],
                color=color,
                marker="o",
                markersize=2.8,
                linewidth=0.75,
                label=group,
            )
        axis.set_xticks(range(6), NETWORK_STAGES)
        hifi.color_time_ticks(axis)
        axis.set_ylabel(ylabel)
        hifi.clean_axes(axis)
        axis.legend(
            frameon=False,
            ncol=2,
            title="Group",
            loc="upper center",
            bbox_to_anchor=(0.5, 1.22),
            fontsize=5.8,
            title_fontsize=6.2,
            handlelength=1.2,
            columnspacing=0.7,
        )
    if completed_fastspar is not None:
        if summary_uses_full_eligible_network:
            panel_c_note = (
                "C: display subsets shown; trajectories use full FastSpar networks."
            )
        else:
            panel_c_note = (
                "C: FastSpar-qualified edges; deterministic source-density subset shown "
                "(full list exported)."
            )
    else:
        panel_c_note = "C proxy: Spearman |rho| > 0.6, raw P < 0.05; not FastSpar."
    fig.text(
        0.35,
        0.355,
        panel_c_note,
        fontsize=5.0,
        color="#555555",
        ha="left",
        va="bottom",
    )

    # D — Bacteroides/Faecalibacterium correlation small multiples.
    correlation_axes: list[mpl.axes.Axes] = []
    for bottom, height in [(0.205, 0.122), (0.020, 0.139)]:
        for left, width in [(0.055, 0.298), (0.363, 0.300), (0.670, 0.300)]:
            correlation_axes.append(fig.add_axes([left, bottom, width, height]))
    correlation_rows: list[dict[str, object]] = []
    correlation_specs = [
        ("shannon", "Shannon"),
        ("richness", "Richness" if raw_bracken else "Richness†"),
        ("simpson", "Simpson"),
    ]
    for row_index, group in enumerate(["LDR", "HDR"]):
        for column_index, (metric, title) in enumerate(correlation_specs):
            axis = correlation_axes[row_index * 3 + column_index]
            correlation_rows.extend(
                draw_correlation_panel(
                    axis,
                    ds,
                    group,
                    metric,
                    title,
                    show_y=column_index == 0,
                )
            )
            if column_index == 2:
                strip_color = STRIP_BLUE if group == "LDR" else STRIP_RED
                axis.add_patch(
                    Rectangle(
                        (1.010, 0),
                        0.065,
                        1,
                        transform=axis.transAxes,
                        clip_on=False,
                        facecolor=strip_color,
                        edgecolor="none",
                    )
                )
                axis.text(
                    1.0425,
                    0.5,
                    group,
                    transform=axis.transAxes,
                    rotation=-90,
                    color="white",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    fontweight="bold",
                )
    correlation_table = pd.DataFrame(correlation_rows)
    correlation_table["q_BH_within_curve_sensitivity"] = np.nan
    for _, indices in correlation_table.groupby(
        ["group", "metric", "genus"], sort=False
    ).groups.items():
        correlation_table.loc[indices, "q_BH_within_curve_sensitivity"] = (
            audit.bh_adjust(correlation_table.loc[indices, "p_raw"])
        )
    correlation_table["symbol_q_sensitivity"] = correlation_table[
        "q_BH_within_curve_sensitivity"
    ].map(hifi.p_symbol)
    correlation_table.to_csv(
        tables / "Fig4D_dominant_genus_correlations.csv", index=False
    )
    fig.text(0.006, 0.345, "D", fontsize=15, ha="left", va="bottom")
    fig.text(
        0.760,
        0.355,
        "D stars: raw Spearman P, n = 15; no FDR.",
        fontsize=5.0,
        color="#555555",
        ha="left",
        va="bottom",
    )

    record = export_pdf_tiff(fig, figures / "Fig4_restored", dpi=600)
    if completed_fastspar is not None:
        fastspar_manifest = completed_fastspar.manifest
        if summary_uses_full_eligible_network:
            panel_c_result_kind = (
                "reconstructed FastSpar caption-rule analysis; deterministic "
                "source-density subsets are used only for the small-network drawings, "
                "while node/edge trajectories use the full eligible networks"
            )
            panel_c_display_boundary = (
                "source edge counts are approximate raster digitizations used only for visual "
                "density in the small-network drawings; every displayed edge is a "
                "caption-qualified FastSpar edge, no edge is synthesized, all eligible nodes "
                "are retained, the full eligible edge list remains the analytical result, and "
                "the node/edge trajectories count that full eligible network"
            )
            panel_c_rule_discrepancy = (
                "the manuscript caption specifies raw P while Methods specify BH-adjusted Q; "
                "both full analytical edge lists are exported; the small-network drawings use "
                "the explicitly labelled source-density display subset, whereas the node/edge "
                "trajectories use the full caption-qualified eligible network"
            )
            panel_c_trajectory_rule = (
                "node/edge trajectories count every node and edge in the full "
                "caption-qualified FastSpar eligible network; the source-density subset "
                "affects only the small-network drawings"
            )
        else:
            panel_c_result_kind = (
                "reconstructed FastSpar caption-rule analysis with transparent "
                "source-density display subset"
            )
            panel_c_display_boundary = (
                "source edge counts are approximate raster digitizations used only for visual "
                "density; every displayed edge is a caption-qualified FastSpar edge, no edge "
                "is synthesized, all eligible nodes are retained, and the full eligible edge "
                "list remains the analytical result"
            )
            panel_c_rule_discrepancy = (
                "the manuscript caption specifies raw P while Methods specify BH-adjusted Q; "
                "both full analytical edge lists are exported, while the network drawing and "
                "edge-count trajectory use the explicitly labelled source-density display subset"
            )
            panel_c_trajectory_rule = (
                "network drawings and node/edge trajectories use the explicitly labelled "
                "source-density display subset"
            )
        panel_c_metadata: dict[str, object] = {
            "panel_C_result_kind": panel_c_result_kind,
            "panel_C_results": str(completed_fastspar.root),
            "panel_C_run_manifest_sha256": fastspar_pipeline.sha256_file(
                completed_fastspar.root / "run_manifest.json"
            ),
            "panel_C_fastspar_version": fastspar_manifest.get("fastspar_version"),
            "panel_C_bootstraps": int(fastspar_manifest["effective_bootstraps"]),
            "panel_C_seed": int(fastspar_manifest["seed"]),
            "panel_C_rho_threshold_strict_gt": float(
                fastspar_manifest["rho_threshold_strict_gt"]
            ),
            "panel_C_alpha_strict_lt": float(fastspar_manifest["alpha_strict_lt"]),
            "panel_C_feature_rule": fastspar_manifest.get("feature_rule"),
            "panel_C_sample_rule": fastspar_manifest.get("sample_rule"),
            "panel_C_boundary": (
                "FastSpar reconstruction from supplied direct bacterial Bracken genus counts; "
                "the original manuscript feature filter, preprocessing choice, pooling rule, "
                "and random seeds remain unavailable, so this is not the original author run"
            ),
            "panel_C_display_rule": (
                f"eligible: |rho| > {float(fastspar_manifest['rho_threshold_strict_gt']):g} "
                f"and exact bootstrap raw P < {float(fastspar_manifest['alpha_strict_lt']):g}; "
                "drawn edges are a deterministic source-raster-density subset that retains "
                "the maximum-spanning forest and then the strongest remaining abs_rho edges"
            ),
            "panel_C_display_is_full_analytic_edge_set": False,
            "panel_C_display_edge_targets": {
                group: {
                    stage: int(SOURCE_DIGITIZED_EDGE_TARGETS[(group, stage)])
                    for stage in NETWORK_STAGES
                }
                for group in ["HDR", "LDR"]
            },
            "panel_C_display_boundary": panel_c_display_boundary,
            "panel_C_trajectory_rule": panel_c_trajectory_rule,
            "panel_C_trajectory_uses_full_eligible_network": bool(
                summary_uses_full_eligible_network
            ),
            "panel_C_sensitivity": (
                "separate Methods-rule edge table uses within-network Benjamini-Hochberg Q "
                f"< {float(fastspar_manifest['alpha_strict_lt']):g} across all candidate pairs"
            ),
            "panel_C_rule_discrepancy": panel_c_rule_discrepancy,
            "panel_C_tables": [
                "Fig4C_fastspar_summary.csv",
                "Fig4C_fastspar_edges_caption_raw_p.csv",
                "Fig4C_fastspar_edges_methods_BH_q.csv",
                "Fig4C_fastspar_edges_display_source_density.csv",
            ],
        }
    else:
        panel_c_metadata = {
            "panel_C_result_kind": "subject-mean Spearman proxy",
            "panel_C_boundary": "stage-specific top-275 nonconstant genera, no prevalence filter; asymptotic subject-mean Spearman proxy, not unavailable FastSpar/SparCC 500-bootstrap output",
            "panel_C_sensitivity": "display uses |rho|>0.6 and raw P<0.05; edge table also records within-network BH across all candidate pairs",
        }
    record.update(
        {
            "source_data": (
                "raw/combine_brk.mpa bacterial species/genus leaves and raw/meta.csv"
                if raw_bracken
                else "species.csv, taxa.csv, meta.csv"
            ),
            "groups": (
                group_definition
                if group_definition is not None
                else "15 lowest and 15 highest C4-C6 mean Shannon subjects"
            ),
            "inference_boundary": "groups are outcome-defined; displayed P/q values are conditional descriptive inference",
            "panel_A": "17 named genera plus residual Others; arithmetic group-time means",
            "panel_A_boundary": (
                "source Bracken genus table recovered"
                if raw_bracken
                else "current taxonomy yields substantially more Others than the unavailable original plotting table"
            ),
            "panel_B": "paired t versus A1; BH across 15 non-baseline times within each metric/group",
            "panel_B_boundary": (
                "source-figure unrarefied bacterial leaf richness; Methods-rarefied result exported as sensitivity"
                if raw_bracken
                else "supplied leaf-taxon richness scale does not reproduce manuscript richness"
            ),
            "panel_D": "time-stratified raw Spearman correlations and raw-P stars; curve-wise BH supplied as sensitivity only",
            **panel_c_metadata,
        }
    )
    return record


def extract_sfig4_reference(data_root: Path, figures: Path) -> dict[str, object]:
    """Extract the native sFig4 TIFF and preserve its exact RGB pixel array."""
    if not shutil.which("tiff2pdf"):
        raise RuntimeError(
            "sFig4 PDF export requires the system command 'tiff2pdf' "
            "(provided by libtiff-tools/libtiff)."
        )
    document = data_root / "抗生素稿件 isme j" / "Supplymentary material.docx"
    member = "word/media/image4.tiff"
    with zipfile.ZipFile(document) as archive:
        source_bytes = archive.read(member)
    source_file_hash = hashlib.sha256(source_bytes).hexdigest()
    with Image.open(io.BytesIO(source_bytes)) as source:
        source_dpi = source.info.get("dpi", (220.0, 220.0))
        dpi_x = float(source_dpi[0])
        dpi_y = float(source_dpi[1])
        rgb = source.convert("RGB")
        pixel_hash = hashlib.sha256(rgb.tobytes()).hexdigest()
        size = rgb.size
        icc_profile = source.info.get("icc_profile")

    tif = figures / "sFig4_restored_source_reference.tif"
    pdf = figures / "sFig4_restored_source_reference.pdf"
    save_options: dict[str, object] = {
        "format": "TIFF",
        "compression": "tiff_lzw",
        "dpi": (dpi_x, dpi_y),
    }
    if icc_profile:
        save_options["icc_profile"] = icc_profile
    rgb.save(tif, **save_options)

    subprocess.run(
        [
            "tiff2pdf",
            "-z",
            "-r",
            "o",
            "-x",
            f"{dpi_x:g}",
            "-y",
            f"{dpi_y:g}",
            "-t",
            "Supplementary Figure 4 source-reference extraction",
            "-o",
            str(pdf),
            str(tif),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        "stem": "sFig4_restored_source_reference",
        "source_document": str(document),
        "source_member": member,
        "source_file_sha256": source_file_hash,
        "source_pixel_sha256": pixel_hash,
        "source_pixels_preserved_exactly": True,
        "width_in": size[0] / dpi_x,
        "height_in": size[1] / dpi_y,
        "dpi": [dpi_x, dpi_y],
        "width_px": size[0],
        "height_px": size[1],
        "formats": ["pdf", "tif"],
        "tiff_mode": "RGB",
        "tiff_compression": "LZW",
        "analysis_status": "source-reference extraction; not recalculated",
        "missing_inputs": "HUMAnN/MetaCyc pathway, gene-family, EC and enzyme abundance tables",
        "upsampling": "none",
    }


def write_readme(output: Path) -> None:
    (output / "README.md").write_text(
        """# Fig4 and sFig4 standalone restoration

Run from `reproduction/`:

```bash
python3 render_fig4_sfig4.py
python3 test_fig4_sfig4.py
```

The final `figures/` directory contains PDF and TIFF files only.

Rendering sFig4 requires the system command `tiff2pdf` (from
`libtiff-tools`/`libtiff`). `pdfinfo`, `pdffonts`, and `pdfimages` are optional
enhanced integrity checks; Python dependencies are listed in `requirements.txt`.

## Fig4_restored

Fig4 is recalculated from `species.csv`, `taxa.csv`, and `meta.csv`. Panel A
uses 17 named genera and `Others`. Panel B uses paired t tests against A1 and
Benjamini-Hochberg correction across the 15 non-baseline time points within
each metric/group family. The available matrix contains leaf taxa rather than
the manuscript's original high-resolution richness table, so that metric is
labelled `Richness†` and the footnote discloses the proxy scale. The supplied
taxonomy produces much more `Others` abundance than the source raster; it is
not visually compressed to force agreement with an unavailable plotting table.
The PDF is a fixed-page mixed export: text and linework remain vector objects,
while heatmap raster artists are embedded at 600 ppi.

The original FastSpar input counts, 500 bootstrap runs, seeds and adjusted-P
edge list were not supplied. Panel C is therefore labelled and documented as a
subject-mean asymptotic Spearman proxy (`|rho| > 0.6`, raw `P < 0.05`), not
FastSpar or SparCC. Its stage-specific top-275 candidate rule has no prevalence
filter; the edge table also includes within-network BH values as a sensitivity
analysis. Panel D displays raw time-stratified Spearman tests (`n=15` per
time/group, no multiplicity adjustment), matching the visible star density;
curve-wise BH values are supplied in the table as sensitivity results.

HDR/LDR are defined using the same C4-C6 Shannon outcome later displayed and
correlated in Fig4. Consequently, the P/q values are conditional descriptive,
exploratory results rather than independent validation of recovery differences.

## sFig4_restored_source_reference

The required HUMAnN/MetaCyc pathway, gene-family, EC and enzyme abundance
tables are absent, so sFig4 cannot be scientifically recalculated. The supplied
supplementary DOCX does contain its native embedded TIFF. This output isolates
that TIFF, converts it to RGB without changing a pixel, and keeps the native
1269 x 1718 resolution and 220 dpi metadata. The TIFF is saved with lossless
LZW compression; the PDF wraps the same raster using lossless Flate encoding.
It is a source-reference extraction, not a reconstructed analysis; no
upsampling or pixel-derived pseudo-data were used. The DOCX image is the best
embedded source supplied here, but it does not establish that an earlier vector
or higher-resolution authoring file never existed.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "fig4_sfig4",
    )
    args = parser.parse_args()
    figures = args.output / "figures"
    tables = args.output / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    set_style()
    dataset = audit.load_dataset(args.data_root)
    records = [
        draw_figure4(dataset, figures, tables),
        extract_sfig4_reference(args.data_root, figures),
    ]
    (args.output / "export_manifest.json").write_text(
        json.dumps(
            {
                "renderer": str(Path(__file__).resolve()),
                "random_seed": SEED,
                "records": records,
                "integrity_note": (
                    "Fig4C is a labelled Spearman proxy, not FastSpar; sFig4 is an "
                    "exact-pixel source-reference extraction because functional inputs are absent."
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    write_readme(args.output)
    print(f"Fig4 and sFig4 deliverables written to {figures.resolve()}")


if __name__ == "__main__":
    main()
