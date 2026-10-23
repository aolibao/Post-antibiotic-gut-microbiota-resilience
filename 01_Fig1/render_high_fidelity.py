#!/usr/bin/env python3
"""Render manuscript-faithful figures without changing the audited analysis outputs.

The statistical calculations are taken from :mod:`run_all`.  This module is a
presentation layer: it reconstructs the panel geometry, visual hierarchy and plot
grammar of the supplied manuscript figures.  Panels whose original method cannot
be recovered from the supplied tables remain explicitly labelled as proxies.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Circle, Ellipse, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
from PIL import Image, ImageDraw, ImageFont
from scipy.interpolate import PchipInterpolator
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import chi2, mannwhitneyu, rankdata, spearmanr, t as student_t, ttest_ind, ttest_rel, wilcoxon

import run_all as audit
import render_selected_panels as selected_panels


SEED = audit.SEED
TIME_ORDER = audit.TIME_ORDER
BASELINE = audit.BASELINE
EXPOSURE = audit.EXPOSURE
RECOVERY = audit.RECOVERY
ENDPOINT = audit.ENDPOINT

# Reverse calibration of sFig2A against the supplied raster.  The manuscript
# does not report a detection threshold for this panel; therefore the exact
# displayed threshold is kept visible in the panel and in sensitivity tables.
SFIG2A_GENUS_PRESENCE_THRESHOLD_PERCENT = 0.0001

BLUE = "#0072B5"
RED = "#BC3C29"
GREEN = "#238443"
LDR = "#0072B5"
HDR = "#BC3C29"
GRID = "#D9D9D9"
PHASE_FACE = {"A": "#F0F8FF", "B": "#FFF6F0", "C": "#F0FFFC"}
TIME_TEXT = {"A": "#2C69A0", "B": "#A52A2A", "C": "#237A45"}
BOX_COLORS = [
    "#9DC1D6", "#78A8C5", "#568BB1", "#3773A5",
    "#EA947C", "#DE735F", "#CA4F48", "#A83435", "#861E23",
    "#9FCD9C", "#83BD82", "#68AD6A", "#4D9954", "#348441", "#1F6E34", "#0D582B",
]
GENUS_COLORS = [
    "#8DA0CB", "#FDB863", "#E31A1C", "#E7B6D2", "#88419D", "#80CDC1",
    "#A6D854", "#2C7FB8", "#B3ABD2", "#FFD11A", "#3AA657", "#DFC27D",
    "#80B1D3", "#BC80BD", "#CCEBC5", "#FFED6F", "#F8766D", "#BDBDBD",
    "#1F9E89", "#F28E2B", "#59A14F", "#AF7AA1", "#EDC948", "#76B7B2",
]


def set_hifi_style() -> None:
    """Typography and line weights inferred from the supplied raster figures."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.4,
            "axes.titlesize": 9.2,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def phase_of(time: str) -> str:
    return "A" if time.startswith("A") else "B" if time in EXPOSURE else "C"


def compact_species_label(name: str) -> str:
    """Return a manuscript-style display label while preserving full CSV names."""
    parts = name.split()
    if len(parts) < 2:
        return name.replace(" ", "_")
    if parts[0].startswith("GGB"):
        return "_".join(parts)
    if parts[0] == "Candidatus" and len(parts) >= 3:
        return f"Ca._{parts[1]}_{parts[-1]}"
    if "bacterium" in parts or "unclassified" in parts:
        family = parts[0].replace("aceae", "")[:10] + "."
        descriptor = "uncl." if "unclassified" in parts else "bact."
        suffix = "-".join(parts[2:]) if len(parts) > 2 else ""
        return f"{family}_{descriptor}" + (f"_{suffix}" if suffix else "")
    if parts[1] == "sp":
        return f"{parts[0][0]}._sp._{'-'.join(parts[2:])}"
    return f"{parts[0][0]}._{'_'.join(parts[1:])}"


def color_time_ticks(ax: mpl.axes.Axes, axis: str = "x") -> None:
    labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for label in labels:
        text = label.get_text()
        if text:
            label.set_color(TIME_TEXT.get(phase_of(text), "black"))


def panel_label(ax: mpl.axes.Axes, text: str, x: float = -0.12, y: float = 1.08, size: float = 15) -> None:
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="normal", va="bottom", ha="left")


def phase_spans(ax: mpl.axes.Axes) -> None:
    ax.axvspan(-0.5, 3.5, color=PHASE_FACE["A"], zorder=-20)
    ax.axvspan(3.5, 8.5, color=PHASE_FACE["B"], zorder=-20)
    ax.axvspan(8.5, 15.5, color=PHASE_FACE["C"], zorder=-20)


def clean_axes(ax: mpl.axes.Axes, grid: bool = False) -> None:
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.8)
    ax.tick_params(color="#333333", labelcolor="#222222")
    if grid:
        ax.grid(axis="y", color=GRID, linestyle="--", linewidth=0.65, alpha=0.85)
        ax.set_axisbelow(True)


def add_antibiotic_bracket(ax: mpl.axes.Axes, y: float = -0.16) -> None:
    trans = ax.get_xaxis_transform()
    ax.plot([4, 8], [y, y], color="#9B2D2D", lw=0.8, transform=trans, clip_on=False)
    ax.plot([4, 4], [y, y + 0.02], color="#9B2D2D", lw=0.8, transform=trans, clip_on=False)
    ax.plot([8, 8], [y, y + 0.02], color="#9B2D2D", lw=0.8, transform=trans, clip_on=False)
    ax.text(6, y - 0.025, "Antibiotic exposure", color="#9B2D2D", ha="center", va="top", transform=trans, fontsize=8.2)


def export(fig: mpl.figure.Figure, stem: Path, dpi: int = 600) -> dict[str, object]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    fig.savefig(png, dpi=dpi, facecolor="white", transparent=False)
    fig.savefig(pdf, facecolor="white", transparent=False)
    fig.savefig(svg, facecolor="white", transparent=False)
    inches = tuple(float(x) for x in fig.get_size_inches())
    plt.close(fig)
    with Image.open(png) as opened:
        rgba = opened.convert("RGBA")
        bg = Image.new("RGBA", rgba.size, (255, 255, 255, 255))
        bg.alpha_composite(rgba)
        rgb = bg.convert("RGB")
        rgb.save(png, dpi=(dpi, dpi), optimize=True)
        pixels = rgb.size
    return {"stem": stem.name, "width_in": inches[0], "height_in": inches[1], "dpi": dpi, "width_px": pixels[0], "height_px": pixels[1]}


def phase_boxplot(
    ax: mpl.axes.Axes,
    arrays: list[np.ndarray],
    ylabel: str,
    pvalues: list[float] | None = None,
    ylim: tuple[float, float] | None = None,
    pvalue_color: str = "#222222",
    box_width: float = 0.78,
) -> None:
    bp = ax.boxplot(
        arrays,
        positions=np.arange(16),
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#333333", "linewidth": 0.9},
        whiskerprops={"color": "#444444", "linewidth": 0.8},
        capprops={"color": "#444444", "linewidth": 0.8},
        boxprops={"edgecolor": "#444444", "linewidth": 0.75},
    )
    for patch, color in zip(bp["boxes"], BOX_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.92)
    ax.set_xticks(range(16), TIME_ORDER)
    color_time_ticks(ax)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    clean_axes(ax, grid=True)
    if pvalues is not None:
        lo, hi = ax.get_ylim()
        y = hi - 0.035 * (hi - lo)
        for i, p in enumerate(pvalues):
            if i == 0 or not np.isfinite(p):
                continue
            label = "p<0.001" if p < 0.001 else f"p={p:.3f}"
            # A small horizontal offset keeps the vertical annotation from being
            # bisected by the boxplot whisker.  The translucent white backing is
            # deliberately tight so it improves legibility without reading as a
            # separate annotation tile.
            ax.text(
                i + 0.14,
                y,
                label,
                rotation=90,
                ha="center",
                va="top",
                fontsize=6.2,
                color=pvalue_color,
                zorder=8,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.08},
            )


def smooth_mean_line(ax: mpl.axes.Axes, means: np.ndarray, color: str = "#C51B2A") -> None:
    x = np.arange(len(means), dtype=float)
    valid = np.isfinite(means)
    dense = np.linspace(x[valid].min(), x[valid].max(), 400)
    curve = PchipInterpolator(x[valid], means[valid])(dense)
    ax.plot(dense, curve, color=color, lw=1.35, ls=(0, (3, 2)), zorder=5)


def streamgraph(
    ax: mpl.axes.Axes,
    table: pd.DataFrame,
    colors: list[str],
    labels: list[str],
) -> None:
    """Draw exact time-point strata joined by manuscript-style alluvia."""
    x = np.arange(table.shape[1], dtype=float)
    values = np.clip(table.to_numpy(dtype=float), 0, None)
    values *= 100.0 / np.maximum(values.sum(axis=0, keepdims=True), 1e-12)
    cumulative = np.vstack([np.zeros(values.shape[1]), np.cumsum(values, axis=0)])
    stratum_width = 1.0 / 3.0

    # Quintic smoothstep connectors keep the first and second derivatives
    # horizontal at both strata, as in the supplied Figure 1D.
    for interval in range(len(x) - 1):
        xx = np.linspace(
            interval + stratum_width / 2,
            interval + 1 - stratum_width / 2,
            80,
        )
        t = np.linspace(0.0, 1.0, 80)
        ease = 6 * t**5 - 15 * t**4 + 10 * t**3
        boundaries = (
            cumulative[:, interval, None] * (1.0 - ease)[None, :]
            + cumulative[:, interval + 1, None] * ease[None, :]
        )
        for row, color in enumerate(colors):
            ax.fill_between(
                xx,
                boundaries[row],
                boundaries[row + 1],
                color=color,
                alpha=0.60,
                linewidth=0,
                zorder=1,
            )

    for column, xpos in enumerate(x):
        for row, (color, label) in enumerate(zip(colors, labels)):
            ax.add_patch(
                Rectangle(
                    (xpos - stratum_width / 2, cumulative[row, column]),
                    stratum_width,
                    values[row, column],
                    facecolor=color,
                    edgecolor="none",
                    alpha=1.0,
                    zorder=2,
                    label=label if column == 0 else "_nolegend_",
                )
            )
    ax.set_xlim(-0.60, 15.60)
    ax.set_ylim(-5, 105)
    ax.set_xticks(range(16), TIME_ORDER, rotation=48, ha="right")
    ax.set_yticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Time Point")
    ax.set_ylabel("Relative Abundance")
    clean_axes(ax)


def workflow_human_icon(
    ax: mpl.axes.Axes,
    x: float,
    y: float,
    scale: float,
    color: str,
) -> None:
    """Small neutral vector human pictogram for the workflow panel."""
    ax.add_patch(Circle((x, y + 1.55 * scale), .27 * scale, fc=color, ec="white", lw=.35, zorder=6))
    ax.plot([x, x], [y + 1.25 * scale, y + .40 * scale], color=color, lw=3.2 * scale / .02,
            solid_capstyle="round", zorder=5)
    ax.plot([x, x - .58 * scale], [y + 1.05 * scale, y + .52 * scale], color=color, lw=1.25,
            solid_capstyle="round", zorder=5)
    ax.plot([x, x + .58 * scale], [y + 1.05 * scale, y + .52 * scale], color=color, lw=1.25,
            solid_capstyle="round", zorder=5)
    ax.plot([x, x - .42 * scale], [y + .44 * scale, y - .20 * scale], color=color, lw=1.45,
            solid_capstyle="round", zorder=5)
    ax.plot([x, x + .42 * scale], [y + .44 * scale, y - .20 * scale], color=color, lw=1.45,
            solid_capstyle="round", zorder=5)


def workflow_mouse_icon(
    ax: mpl.axes.Axes,
    x: float,
    y: float,
    scale: float,
) -> None:
    """Compact mouse silhouette redrawn with vector primitives."""
    body = "#5B5552"
    ax.add_patch(Ellipse((x, y), 1.45 * scale, .72 * scale, angle=-10, fc=body, ec="none", zorder=5))
    ax.add_patch(Circle((x + .72 * scale, y + .12 * scale), .30 * scale, fc=body, ec="none", zorder=6))
    ax.add_patch(Circle((x + .62 * scale, y + .39 * scale), .15 * scale, fc="#827873", ec="none", zorder=6))
    ax.add_patch(Circle((x + .92 * scale, y + .17 * scale), .035 * scale, fc="#111111", ec="none", zorder=7))
    verts = [
        (x - .67 * scale, y + .04 * scale),
        (x - 1.15 * scale, y + .25 * scale),
        (x - 1.35 * scale, y + .62 * scale),
        (x - 1.58 * scale, y + .50 * scale),
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    ax.add_patch(PathPatch(MplPath(verts, codes), fill=False, ec="#88736C", lw=.85, zorder=4))
    # A tiny injection cue, matching the intervention iconography in the source.
    ax.plot([x - .18 * scale, x + .16 * scale], [y + .65 * scale, y + .34 * scale],
            color="#6B737C", lw=1.2, zorder=8)
    ax.plot([x - .44 * scale, x - .18 * scale], [y + .90 * scale, y + .65 * scale],
            color="#98C9D5", lw=2.4, solid_capstyle="butt", zorder=8)


def workflow_molecule(
    ax: mpl.axes.Axes,
    x: float,
    y: float,
    scale: float,
    label: str,
    fused: bool = False,
) -> None:
    """Minimal skeletal-formula motif used only as a workflow icon."""
    angles = np.linspace(0, 2 * np.pi, 7) + np.pi / 6
    first = np.c_[x + scale * np.cos(angles), y + scale * np.sin(angles)]
    ax.plot(first[:, 0], first[:, 1], color="#343434", lw=.65, zorder=4)
    ax.plot([x + .87 * scale, x + 1.38 * scale], [y + .50 * scale, y + .82 * scale], color="#343434", lw=.65)
    ax.text(x + 1.52 * scale, y + .90 * scale, "OH", color="#C73A3A", fontsize=3.5, ha="left", va="center")
    if fused:
        second_x = x - 1.35 * scale
        second = np.c_[second_x + .82 * scale * np.cos(angles), y + .82 * scale * np.sin(angles)]
        ax.plot(second[:, 0], second[:, 1], color="#343434", lw=.65, zorder=4)
        ax.plot([x - .45 * scale, second_x + .70 * scale], [y, y], color="#343434", lw=.65)
        ax.text(second_x - .12 * scale, y + .25 * scale, "NH", fontsize=3.4,
                color="#2F5F9A", ha="center", va="center")
    else:
        ax.text(x - .82 * scale, y + .50 * scale, "N", fontsize=3.5,
                color="#2F5F9A", ha="center", va="center")
    ax.text(x, y - 1.62 * scale, label, fontsize=5.0, ha="center", va="top", color="#444444")


def workflow_qpcr_icon(ax: mpl.axes.Axes, x: float, y: float, scale: float) -> None:
    ax.add_patch(Rectangle((x - .60 * scale, y - .42 * scale), 1.20 * scale, .78 * scale,
                           fc="#4D5664", ec="#2F3540", lw=.45, zorder=4))
    ax.add_patch(Polygon([(x - .53 * scale, y + .36 * scale), (x + .37 * scale, y + .36 * scale),
                          (x + .58 * scale, y + .58 * scale), (x - .30 * scale, y + .58 * scale)],
                         closed=True, fc="#6F7B8D", ec="#2F3540", lw=.4, zorder=5))
    ax.add_patch(Rectangle((x - .38 * scale, y - .18 * scale), .42 * scale, .32 * scale,
                           fc="#8CB9D7", ec="white", lw=.35, zorder=6))
    ax.add_patch(Circle((x + .36 * scale, y - .16 * scale), .055 * scale, fc="#72C78C", ec="none", zorder=6))


def workflow_sequencer_icon(ax: mpl.axes.Axes, x: float, y: float, scale: float) -> None:
    for offset, color in [(0, "#5B6570"), (.12, "#747F8A"), (.24, "#919AA2")]:
        ax.add_patch(Rectangle((x - .58 * scale + offset * scale, y - .42 * scale + offset * scale),
                               1.05 * scale, .72 * scale, fc=color, ec="white", lw=.35, zorder=4))
    for index, base in enumerate("ACGT"):
        ax.text(x + (.02 + index * .14) * scale, y - (.13 + index * .05) * scale, base,
                color=["#39A96B", "#2878B5", "#D94A45", "#E2A72E"][index],
                fontsize=3.4, rotation=25, ha="center", va="center", zorder=7)


def draw_workflow_panel(ax: mpl.axes.Axes) -> None:
    """Vector redraw of the study workflow; no analytical result is invented."""
    ax.set_axis_off()
    sections = [(0.00, 0.29, "Stage 1: Cohort grouping"), (0.32, 0.69, "Stage 2: Bioinformatic analysis"), (0.72, 1.00, "Stage 3: In vivo experiment")]
    for left, right, title in sections:
        ax.add_patch(Rectangle((left, 0.88), right - left, 0.08, facecolor="#9FD0D2", edgecolor="none"))
        ax.add_patch(Rectangle((left, 0.02), right - left, 0.84, facecolor="white", edgecolor="#21A9E0", lw=1.0, linestyle=(0, (2, 2))))
        ax.text((left + right) / 2, 0.92, title, ha="center", va="center", fontsize=10.5)
    for x in (0.305, 0.705):
        ax.annotate("", xy=(x + 0.025, 0.92), xytext=(x - 0.025, 0.92), arrowprops={"arrowstyle": "simple", "fc": "#D2D2D2", "ec": "white", "lw": 1.2, "mutation_scale": 20})

    # Cohort and ranking path.
    ax.text(0.145, 0.82, "Metagenomic data", ha="center", fontsize=8.2, fontweight="bold")
    ax.text(0.145, 0.785, "Yaffe, E. et al.", ha="center", fontsize=4.8,
            color="#527E86", fontstyle="italic")
    workflow_human_icon(ax, .145, .675, .030, "#D9C9BE")
    ax.text(0.145, 0.635, "CIP exposure cohort, n=65", ha="center", fontsize=7.9)
    ax.plot([0.09, 0.20], [0.585, 0.585], color="#333333", lw=1.0)
    for x, lab in zip([0.09, 0.145, 0.20], ["Pre-", "During-", "Post-"]):
        ax.plot([x, x], [0.568, 0.602], color="#333333", lw=0.8)
        ax.text(x, 0.542, lab, ha="center", fontsize=6.6)
    ax.text(0.075, 0.475, "Sample correlation\nanalysis", ha="center", fontsize=6.2)
    ax.text(0.215, 0.475, "Diversity ranking\nanalysis", ha="center", fontsize=6.2)
    correlation_icon = np.array(
        [[1.0, .55, -.20, -.45], [.55, 1.0, .15, -.35], [-.20, .15, 1.0, .70], [-.45, -.35, .70, 1.0]]
    )
    ax.imshow(correlation_icon, cmap="RdBu_r", vmin=-1, vmax=1,
              extent=(0.040, 0.110, 0.350, 0.440), aspect="auto")
    for pos in np.linspace(.040, .110, 5):
        ax.plot([pos, pos], [.350, .440], color="white", lw=.22, alpha=.8)
        ax.plot([.040, .110], [.350 + (pos - .040) / .070 * .090] * 2,
                color="white", lw=.22, alpha=.8)
    ranked_x = np.linspace(.178, .247, 18)
    ranked_height = np.linspace(.020, .100, 18) + .004 * np.sin(np.arange(18) * 1.3)
    ax.bar(ranked_x, ranked_height, bottom=.350, width=.0038,
           color="#69B88B", edgecolor="white", linewidth=.12)
    ax.annotate("", xy=(0.145, 0.292), xytext=(0.145, 0.338), arrowprops={"arrowstyle": "->", "lw": .9})
    ax.plot([0.075, 0.215], [0.292, 0.292], color="#333333", lw=.9)
    ax.plot([0.075, 0.075], [0.292, 0.254], color="#333333", lw=.9)
    ax.plot([0.215, 0.215], [0.292, 0.254], color="#333333", lw=.9)
    workflow_human_icon(ax, .075, .205, .016, "#2B83BA")
    workflow_human_icon(ax, .215, .205, .016, "#D65A73")
    ax.text(0.075, 0.115, "P1 = Low-diversity\nrecovery cohort I\n(n=15)", ha="center", va="center", fontsize=6.2, fontweight="bold")
    ax.text(0.215, 0.115, "P4 = High-diversity\nrecovery cohort II\n(n=15)", ha="center", va="center", fontsize=6.2, fontweight="bold")

    # Bioinformatics modules.
    modules = [
        (0.40, 0.71, "Microbial composition\nanalysis", "bars"),
        (0.56, 0.71, "Microbial diversity\nanalysis", "scatter"),
        (0.40, 0.42, "Differential abundance\nanalysis", "heat"),
        (0.56, 0.42, "Interaction network\nanalysis", "network"),
        (0.40, 0.14, "Functional\nanalysis", "volcano"),
        (0.56, 0.14, "Correlation\nanalysis", "corr"),
    ]
    for x, y, label, kind in modules:
        ax.text(x, y + 0.08, label, ha="center", fontsize=7.25)
        if kind == "bars":
            for k in range(4):
                left = x - 0.05
                yy = y - 0.015 * k
                ax.add_patch(Rectangle((left, yy), 0.022, 0.012, color=GENUS_COLORS[k]))
                ax.add_patch(Rectangle((left + 0.022, yy), 0.030, 0.012, color=GENUS_COLORS[k + 5]))
                ax.add_patch(Rectangle((left + 0.052, yy), 0.048, 0.012, color=GENUS_COLORS[k + 9]))
        elif kind == "scatter":
            rng = np.random.default_rng(10)
            for col, cx in zip(["#53B6AD", "#E56B5D", "#8C6BB1"], [x - 0.035, x + 0.035, x]):
                ax.scatter(cx + rng.normal(0, 0.012, 35), y - 0.025 + rng.normal(0, 0.025, 35), s=2, color=col, alpha=0.7)
        elif kind == "heat":
            rng = np.random.default_rng(11)
            ax.imshow(rng.normal(size=(6, 6)), cmap="RdYlBu_r", extent=(x - 0.05, x + 0.05, y - 0.08, y + 0.015), aspect="auto")
        elif kind == "network":
            pts = np.array([[-.04, .00], [-.01, .04], [.03, .03], [.05, -.01], [.01, -.06], [-.045, -.045]]) + [x, y - .03]
            for i, j in [(0, 1), (0, 3), (1, 2), (1, 4), (2, 3), (3, 4), (4, 5), (0, 5)]:
                ax.plot(pts[[i, j], 0], pts[[i, j], 1], color="#333", lw=.8)
            ax.scatter(pts[:, 0], pts[:, 1], s=25, c=[BLUE, "#F2C14E", BLUE, "#B75D69", BLUE, "#B75D69"], edgecolors="#333", lw=.5)
        elif kind == "volcano":
            rng = np.random.default_rng(12)
            xx = rng.normal(0, .032, 120); yy = np.abs(xx) * 1.7 + rng.random(120) * .035
            ax.scatter(x + xx, y - .065 + yy, s=2, c=np.where(xx > .025, "#F28E2B", np.where(xx < -.025, "#4E9BD4", "#AAAAAA")), alpha=.8)
        else:
            ptsx = np.linspace(x - .05, x + .05, 8); ptsy = y - .07 + .7 * (ptsx - x + .05) + np.array([0,.01,-.005,.008,-.002,.012,-.006,.006])
            ax.scatter(ptsx, ptsy, s=9, color="#7A4E9D"); ax.plot([x-.05,x+.05],[y-.07,y], color="#B99CD1", lw=.8)

    # Small skeletal-formula motifs and causal arrows reproduce the source
    # workflow's molecular-intervention visual density without using raster art.
    workflow_molecule(ax, .650, .700, .012, "NAM")
    workflow_molecule(ax, .650, .420, .011, "Trp", fused=True)
    ax.annotate("", xy=(.650, .490), xytext=(.650, .630),
                arrowprops={"arrowstyle": "->", "lw": .7, "color": "#555555"})
    ax.annotate("", xy=(.585, .180), xytext=(.632, .360),
                arrowprops={"arrowstyle": "->", "lw": .7, "color": "#555555"})

    # In-vivo branch.
    workflow_mouse_icon(ax, .790, .765, .030)
    workflow_mouse_icon(ax, .920, .765, .030)
    ax.text(0.79, 0.695, "Saline", ha="center", fontsize=7.5)
    ax.text(0.92, 0.695, "Trp or NAM", ha="center", fontsize=7.5)
    ax.plot([0.79, 0.92], [0.650, 0.650], color="#333", lw=.9)
    ax.plot([0.79, 0.79], [0.650, 0.678], color="#333", lw=.9)
    ax.plot([0.92, 0.92], [0.650, 0.678], color="#333", lw=.9)
    ax.annotate("", xy=(0.855, 0.535), xytext=(0.855, 0.650), arrowprops={"arrowstyle": "->", "lw": .9})
    for dx, dy, angle in [(-.018, .008, 20), (0, 0, -10), (.018, .010, 15), (.007, -.013, -25)]:
        ax.add_patch(Ellipse((.855 + dx, .500 + dy), .018, .010, angle=angle,
                             fc="#72513D", ec="#4F372A", lw=.25, zorder=5))
    ax.text(0.855, 0.455, "Fecal", ha="center", fontsize=7.6)
    ax.plot([0.79, 0.92], [0.405, 0.405], color="#333", lw=.9)
    ax.annotate("", xy=(0.79, 0.320), xytext=(0.79, 0.405), arrowprops={"arrowstyle": "->", "lw": .9})
    ax.annotate("", xy=(0.92, 0.320), xytext=(0.92, 0.405), arrowprops={"arrowstyle": "->", "lw": .9})
    workflow_qpcr_icon(ax, .790, .267, .050)
    workflow_sequencer_icon(ax, .920, .267, .052)
    ax.text(0.79, 0.175, "RT-qPCR", ha="center", fontsize=7.6)
    ax.text(0.92, 0.175, "16S rDNA\nsequencing", ha="center", fontsize=7.6)
    # Keep the panel label inside the fixed canvas; y=.99 clipped the glyph cap.
    panel_label(ax, "A", x=-0.015, y=0.95, size=16)


def figure1_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    stats = audit.paired_time_stats(ds, "shannon").set_index("canon").reindex(TIME_ORDER)
    jac_stats = []
    controls = ds.metrics.loc[ds.metrics["abx"] == "N"]
    for time in TIME_ORDER:
        if time == "A1":
            jac_stats.append(np.nan)
            continue
        now = exposed.loc[exposed["canon"] == time, "jaccard_to_A1"]
        reference = controls.loc[controls["canon"] == time, "jaccard_to_A1"]
        jac_stats.append(float(ttest_ind(now, reference, equal_var=False).pvalue))
    pd.DataFrame(
        {
            "canon": TIME_ORDER,
            "p_visual_annotation": jac_stats,
            "contrast": "CIP-exposed vs unexposed control at the same canonical time",
            "test": "Welch t-test; exploratory visual-layer annotation because manuscript contrast is unspecified",
        }
    ).to_csv(tables / "figure1C_visual_annotation_tests.csv", index=False)

    fig = plt.figure(figsize=(7.778, 8.600))
    # Panel boxes measured from the supplied 7778 x 8600 source TIFF.  Fixed
    # axes preserve the original physical panel sizes on the 7.778 x 8.600 in
    # canvas; in particular the former GridSpec made panel D too small.
    ax_a = fig.add_axes([0.035, 0.660, 0.950, 0.325])
    ax_b = fig.add_axes([0.04551, 0.35395, 0.37233, 0.25279])
    # Preserve the measured top edge while reserving enough bottom space for
    # the antibiotic bracket and its caption below the tick labels.
    ax_c = fig.add_axes([0.04474, 0.05500, 0.37310, 0.23547])
    ax_d = fig.add_axes([0.49087, 0.19372, 0.49344, 0.41349])
    draw_workflow_panel(ax_a)

    shannon_arrays = [exposed.loc[exposed["canon"] == t, "shannon"].to_numpy() for t in TIME_ORDER]
    phase_boxplot(ax_b, shannon_arrays, "Shannon-index", stats["q_BH_vs_A1"].to_list(), ylim=(1.85, 5.05))
    ax_b.set_ylabel("Shannon-index", fontsize=8.5)
    ax_b.yaxis.set_label_coords(-.075, .5)
    smooth_mean_line(ax_b, np.array([np.mean(x) for x in shannon_arrays]))
    add_antibiotic_bracket(ax_b)
    panel_label(ax_b, "B", x=-0.06, y=1.02)

    jac_arrays = [exposed.loc[exposed["canon"] == t, "jaccard_to_A1"].to_numpy() for t in TIME_ORDER]
    phase_boxplot(ax_c, jac_arrays, "Paired Jaccard Dist to A1", jac_stats, ylim=(-0.04, 0.91))
    ax_c.set_ylabel("Paired Jaccard Dist to A1", fontsize=8.5)
    ax_c.yaxis.set_label_coords(-.075, .5)
    smooth_mean_line(ax_c, np.array([np.mean(x) for x in jac_arrays]))
    add_antibiotic_bracket(ax_c, y=-0.150)
    panel_label(ax_c, "C", x=-0.06, y=1.02)

    ids_by_time = exposed.groupby("canon", observed=True)["ids"].apply(list)
    mean_species = pd.DataFrame({t: ds.species_named[ids_by_time[t]].mean(axis=1) for t in TIME_ORDER})
    # The manuscript legend is recovered exactly when feature selection uses
    # all 1,040 samples, while the plotted longitudinal means use the 960 CIP
    # samples.  The legacy SGB branch retains its original exposed-only rule.
    top = (
        ds.species_named.mean(axis=1).nlargest(30).index
        if "grouping_endpoint_shannon" in ds.group_table.columns
        else mean_species.mean(axis=1).nlargest(30).index
    )
    stream = mean_species.loc[top].copy()
    stream.loc["Others"] = (100.0 - stream.sum(axis=0)).clip(lower=0)
    # Put Others at the bottom, matching the supplied plot, then reverse the remaining rank.
    order = ["Others"] + list(stream.drop(index="Others").index[::-1])
    stream = stream.loc[order]
    stream.to_csv(tables / "figure1_streamgraph_high_fidelity.csv")
    cols = ["#D3D3D3", *selected_panels.FIG1D_COLORS[::-1]]
    streamgraph(ax_d, stream, cols, order)
    panel_label(ax_d, "D", x=-0.08, y=1.02)
    handles, labels = ax_d.get_legend_handles_labels()
    # Three explicit columns reproduce the manuscript's 10 / 10 / 11 split.
    ranked_handles = handles[::-1]
    ranked_labels = labels[::-1]
    for left, start, stop in [(.4909, 0, 10), (.6662, 10, 20), (.8049, 20, 31)]:
        legend = fig.legend(
            ranked_handles[start:stop], ranked_labels[start:stop],
            loc="upper left", bbox_to_anchor=(left, .1250), ncol=1,
            frameon=False, fontsize=5.8, handlelength=.86, handleheight=.86,
            handletextpad=.38, labelspacing=.10, borderaxespad=0,
        )
        for txt in legend.get_texts():
            if txt.get_text() != "Others":
                txt.set_fontstyle("italic")
    return export(fig, out / "Figure1_high_fidelity")


def p_symbol(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def trajectory_panel(
    ax: mpl.axes.Axes,
    ds: audit.Dataset,
    metric: str,
    ylabel: str,
    show_legend: bool = True,
) -> list[dict[str, object]]:
    phase_spans(ax)
    rows: list[dict[str, object]] = []
    xx = np.arange(16)
    for group, color in [("LDR", LDR), ("HDR", HDR)]:
        part = ds.metrics.loc[ds.metrics["recovery_group"] == group]
        summary = part.groupby("canon", observed=True)[metric].agg(["mean", "sem"]).reindex(TIME_ORDER)
        ax.fill_between(xx, summary["mean"] - summary["sem"], summary["mean"] + summary["sem"], color=color, alpha=0.20, lw=0)
        ax.plot(xx, summary["mean"], color=color, lw=1.1, marker="o", markersize=2.9, markeredgecolor="white", markeredgewidth=0.35, label=group)
    ymax = -np.inf
    ymin = np.inf
    for i, time in enumerate(TIME_ORDER):
        low = ds.metrics.loc[(ds.metrics["recovery_group"] == "LDR") & (ds.metrics["canon"] == time), metric]
        high = ds.metrics.loc[(ds.metrics["recovery_group"] == "HDR") & (ds.metrics["canon"] == time), metric]
        test = ttest_ind(high, low, equal_var=False)
        rows.append({"metric": metric, "canon": time, "p_welch_raw": float(test.pvalue), "symbol": p_symbol(float(test.pvalue))})
        ymax = max(ymax, float(low.mean() + low.sem()), float(high.mean() + high.sem()))
        ymin = min(ymin, float(low.mean() - low.sem()), float(high.mean() - high.sem()))
    span = max(ymax - ymin, 1e-6)
    ax.set_ylim(ymin - 0.12 * span, ymax + 0.26 * span)
    star_y = ymax + 0.11 * span
    for i, row in enumerate(rows):
        if row["symbol"]:
            ax.text(i, star_y, row["symbol"], ha="center", va="center", fontsize=7.0, fontweight="bold")
    ax.set_xlim(-0.5, 15.5)
    ax.set_xticks(range(16), TIME_ORDER)
    color_time_ticks(ax)
    ax.set_ylabel(ylabel)
    clean_axes(ax)
    if show_legend:
        ax.legend(frameon=False, loc="lower right", handlelength=1.6, handletextpad=0.35, borderpad=0.1)
    return rows


def figure2_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    raw_bracken = "grouping_endpoint_shannon" in ds.group_table.columns
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"]
    pivot = exposed.pivot(index="subject", columns="canon", values="shannon").reindex(columns=TIME_ORDER)
    corr = pivot.corr(method="spearman")
    fig = plt.figure(figsize=(8.137, 9.447))
    # Fixed panel boxes measured from the supplied 2441 x 2834 source TIFF.
    ax_a = fig.add_axes([.056125, .695836, .324867, .279817])
    # Keep a two-pixel safety inset from the fixed page edge at 600 dpi.
    ax_b = fig.add_axes([.495289, .696895, .502711, .275935])
    ax_c = fig.add_axes([.049980, .435427, .286768, .214185])
    ax_d = fig.add_axes([.378943, .437191, .287177, .214185])
    ax_e = fig.add_axes([.702581, .434368, .287177, .214185])
    ax_f = fig.add_axes([.057354, .108680, .918066, .288638])

    # A — correlation heat map.
    corr_cmap = LinearSegmentedColormap.from_list("manuscript_corr", ["#0B4F8A", "#67A9CF", "#F7F7F7", "#F4A582", "#87001E"])
    im = ax_a.imshow(corr.to_numpy(), vmin=0.2, vmax=1.0, cmap=corr_cmap, interpolation="nearest", aspect="equal")
    ax_a.set_xticks(range(16), TIME_ORDER, rotation=0)
    ax_a.set_yticks(range(16), TIME_ORDER)
    color_time_ticks(ax_a, "x"); color_time_ticks(ax_a, "y")
    ax_a.set_xticks(np.arange(-0.5, 16, 1), minor=True); ax_a.set_yticks(np.arange(-0.5, 16, 1), minor=True)
    ax_a.grid(which="minor", color="white", lw=0.8); ax_a.tick_params(which="minor", length=0)
    ax_a.add_patch(Rectangle((12.5, 12.5), 3, 3, fill=False, ec="#111111", lw=1.2, ls=(0, (4, 2))))
    clean_axes(ax_a)
    cax_a = fig.add_axes([.3855, .6960, .0168, .2790])
    cbar = fig.colorbar(im, cax=cax_a)
    # Keep this title close to its colour bar so it does not collide with the
    # neighbouring panel-B y label on the manuscript's fixed canvas.
    cbar.set_label("Spearman r", rotation=90, labelpad=-4, fontsize=7.4)
    cbar.outline.set_linewidth(0.5)
    panel_label(ax_a, "A", x=-0.13, y=0.95)

    # B — ranked grouping panel.  The default branch preserves the supplied
    # endpoint-Shannon reconstruction.  Reviewer-response reruns can inject an
    # alternative scalar definition through the four ``alternative_*`` columns
    # without changing any other Figure 2 panel or its fixed geometry.
    alternative_grouping = "alternative_metric_value" in ds.group_table.columns
    if alternative_grouping:
        direction = str(ds.group_table["alternative_high_recovery_direction"].iloc[0])
        if direction not in {"larger", "smaller"}:
            raise ValueError(f"Unknown recovery direction for Figure 2B: {direction!r}")
        # Keep the manuscript's visual grammar (LDR at left, HDR at right).
        # Thus smaller-is-better metrics are sorted in descending order.
        rank_table = ds.group_table.sort_values(
            ["alternative_metric_value", "subject_sort_key"],
            ascending=[direction == "larger", True],
            kind="mergesort",
        )
        x = np.arange(len(rank_table))
        values_b = rank_table["alternative_metric_value"].to_numpy(dtype=float)
        fills = [
            "#86C5DA" if group == "LDR" else "#E58A7B" if group == "HDR" else "#D7D7D7"
            for group in rank_table["recovery_group"]
        ]
        edges = [
            LDR if group == "LDR" else HDR if group == "HDR" else "#8A8A8A"
            for group in rank_table["recovery_group"]
        ]
        ax_b.bar(
            x,
            values_b,
            width=0.76,
            color=fills,
            edgecolor=edges,
            linewidth=0.45,
            zorder=3,
        )
        finite = values_b[np.isfinite(values_b)]
        if not len(finite):
            raise ValueError("Alternative Figure 2B grouping metric contains no finite values")
        ymin, ymax = float(finite.min()), float(finite.max())
        span = max(ymax - ymin, max(abs(ymax), 1.0) * 0.08)
        y0 = min(0.0, ymin - 0.06 * span)
        y1 = ymax + 0.22 * span
        ax_b.set_ylim(y0, y1)
        box_bottom = y0 + 0.02 * (y1 - y0)
        box_height = 0.86 * (y1 - y0)
        ax_b.add_patch(Rectangle((-0.75, box_bottom), 15.5, box_height, fill=False, ec="#2D9CDB", lw=1.2, ls=(0, (2, 1))))
        ax_b.add_patch(Rectangle((44.25, box_bottom), 15.5, box_height, fill=False, ec="#BB2E3C", lw=1.2, ls=(0, (2, 1))))
        label_y = y1 - 0.04 * (y1 - y0)
        ax_b.text(7, label_y, "LDR", color="#2D9CDB", ha="center", va="top", fontsize=10)
        ax_b.text(29.5, label_y, "Middle", color="#777777", ha="center", va="top", fontsize=8.2)
        ax_b.text(52, label_y, "HDR", color="#F05A78", ha="center", va="top", fontsize=10)
        metric_label = str(ds.group_table["alternative_metric_label"].iloc[0])
        definition_id = str(ds.group_table["alternative_definition_id"].iloc[0])
        ax_b.set_xlabel(f"Subjects (ranked by {metric_label})")
        ax_b.set_ylabel(metric_label)
        ax_b.text(
            0.5,
            1.015,
            f"{definition_id}: blue/red bars are the selected n=15 groups",
            transform=ax_b.transAxes,
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#555555",
        )
        rank_table[
            [
                "alternative_definition_id",
                "alternative_metric_label",
                "alternative_high_recovery_direction",
                "alternative_metric_value",
                "recovery_group",
            ]
        ].to_csv(tables / "figure2_panelB_grouping.csv")
    else:
        rank_sort = "grouping_endpoint_shannon" if raw_bracken else "endpoint_shannon"
        rank_table = ds.group_table.sort_values(rank_sort)
        endpoint = exposed.loc[exposed["canon"].isin(ENDPOINT)].pivot(index="subject", columns="canon", values="shannon").loc[rank_table.index]
        x = np.arange(len(rank_table)); width = 0.18
        for j, (time, color) in enumerate(zip(ENDPOINT, ["#B4E1AF", "#67BC7D", "#338A56"])):
            ax_b.bar(x + (j - 1.5) * width, endpoint[time], width, color=color, linewidth=0, label=time, zorder=2)
        mean_edge = ["#4DBBD5" if g == "LDR" else "#E64B35" if g == "HDR" else "#3C5488" for g in rank_table["recovery_group"]]
        ax_b.bar(x + 1.5 * width, rank_table["endpoint_shannon"], width, facecolor="white", edgecolor=mean_edge, linewidth=0.65, label="Mean", zorder=3)
        y0 = max(1.9, float(endpoint.min().min()) - 0.08)
        y1 = float(endpoint.max().max()) + 0.23
        ax_b.set_ylim(y0, y1)
        ax_b.add_patch(Rectangle((-0.75, y0 + 0.01), 15.5, y1 - y0 - 0.24, fill=False, ec="#2D9CDB", lw=1.2, ls=(0, (2, 1))))
        ax_b.add_patch(Rectangle((44.25, y0 + 0.01), 15.5, y1 - y0 - 0.08, fill=False, ec="#BB2E3C", lw=1.2, ls=(0, (2, 1))))
        ax_b.text(7, y1 - 0.68, "LDR", color="#2D9CDB", ha="center", fontsize=10)
        ax_b.text(52, y1 - 0.12, "HDR", color="#F05A78", ha="center", fontsize=10)
        ax_b.set_xlabel("Subjects (Ranked by endpoint shannon)")
        ax_b.set_ylabel("Shannon Index")
        handles, labels = ax_b.get_legend_handles_labels()
        handles += [mpl.lines.Line2D([], [], marker="s", linestyle="none", color=LDR, markersize=6), mpl.lines.Line2D([], [], marker="s", linestyle="none", color=HDR, markersize=6)]
        labels += ["LDR", "HDR"]
        ax_b.legend(handles, labels, loc="upper left", ncol=2, frameon=False, columnspacing=.7, handletextpad=.3)
    ax_b.set_xlim(-1, 60); ax_b.set_xticks([])
    ax_b.grid(axis="y", color="#E4E4E4", ls="--", lw=.55); ax_b.set_axisbelow(True)
    clean_axes(ax_b)
    panel_label(ax_b, "B", x=-0.055, y=0.95)

    # C — baseline/endpoint association with actual supplied-table statistics.
    corr_rows = []
    for group, color, marker in [("LDR", LDR, "o"), ("HDR", "#D2644F", "o")]:
        part = rank_table.query("recovery_group == @group")
        ax_c.scatter(part["baseline_shannon"], part["endpoint_shannon"], s=25, color=color, edgecolor="white", linewidth=.35, zorder=3)
        coef = np.polyfit(part["baseline_shannon"], part["endpoint_shannon"], 1)
        xx = np.linspace(part["baseline_shannon"].min(), part["baseline_shannon"].max(), 100)
        ax_c.plot(xx, np.polyval(coef, xx), color=color, lw=1.0, ls=(0, (3, 2)))
        rho, p = audit.safe_spearman(part["baseline_shannon"].to_numpy(), part["endpoint_shannon"].to_numpy())
        corr_rows.append({"group": group, "rho": rho, "p": p, "n": len(part)})
    tables.joinpath("figure2_panelC_high_fidelity.csv").write_text(pd.DataFrame(corr_rows).to_csv(index=False), encoding="utf-8")
    ax_c.text(.05, .95, f"LDR: r={corr_rows[0]['rho']:.2f}, p={corr_rows[0]['p']:.4f}", transform=ax_c.transAxes, color=LDR, va="top", fontsize=8.5)
    ax_c.text(.05, .87, f"HDR: r={corr_rows[1]['rho']:.2f}, p={corr_rows[1]['p']:.4f}", transform=ax_c.transAxes, color=HDR, va="top", fontsize=8.5)
    ax_c.set_xlabel("Baseline Mean Shannon (A1–A4)")
    ax_c.set_ylabel("Endpoint Mean Shannon (C4–C6)")
    ax_c.yaxis.set_label_coords(-.11, .5)
    ax_c.legend(handles=[mpl.lines.Line2D([], [], marker="o", linestyle="none", color=HDR, label="HDR"), mpl.lines.Line2D([], [], marker="o", linestyle="none", color=LDR, label="LDR")], frameon=False, loc="lower right", handletextpad=.3, borderpad=.1)
    clean_axes(ax_c)
    panel_label(ax_c, "C", x=-0.15, y=1.02)

    # D/E — manuscript panel geometry, real values and actual Welch-test stars.
    richness_label = "Richness" if raw_bracken else "Observed leaf taxa†"
    stats_rows = trajectory_panel(ax_d, ds, "richness", richness_label)
    stats_rows += trajectory_panel(ax_e, ds, "shannon", "Shannon")
    ax_d.yaxis.set_label_coords(-.11, .5)
    pd.DataFrame(stats_rows).to_csv(tables / "figure2_trajectory_tests_high_fidelity.csv", index=False)
    panel_label(ax_d, "D", x=-0.16, y=1.02); panel_label(ax_e, "E", x=-0.16, y=1.02)
    if raw_bracken:
        fig.text(
            .52, .407,
            "Bacterial Bracken species leaves; source-figure branch (Methods-rarefied sensitivity exported separately)",
            fontsize=5.5, color="#555555", ha="center", va="center",
        )
    else:
        fig.text(
            .52, .407,
            "† Leaf-taxon count in supplied matrix; original richness scale unavailable",
            fontsize=5.8, color="#555555", ha="center", va="center",
        )

    # F — latest manuscript-calibrated heat map.  Values remain bounded
    # rank-biserial effects; the symmetric gamma lookup restores only the
    # source raster's colour occupancy, never the numerical scale.
    matrix, displayed = selected_panels.figure2f_tables(ds, tables)
    matrix.to_csv(tables / "figure2_effect_heatmap_high_fidelity.csv")
    pd.DataFrame(
        {
            "normalized_coordinate": np.linspace(
                -1.0, 1.0, len(selected_panels.FIG2F_RDBU_ANCHORS)
            ),
            "source_anchor_hex": selected_panels.FIG2F_RDBU_ANCHORS,
        }
    ).to_csv(tables / "Fig2F_source_palette_anchors.csv", index=False)
    values = matrix.to_numpy(dtype=float)
    is_lda = "lda_score_LDR_over_HDR" in displayed.columns
    heatmap_cmap = selected_panels.manuscript_fig2f_cmap(gamma=1.0 if is_lda else selected_panels.FIG2F_COLOR_GAMMA)
    heatmap_norm = (
        TwoSlopeNorm(vmin=-5.0, vcenter=0.0, vmax=5.0)
        if is_lda
        else TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    )
    fim = ax_f.pcolormesh(
        np.arange(values.shape[1] + 1),
        np.arange(values.shape[0] + 1),
        values,
        cmap=heatmap_cmap,
        norm=heatmap_norm,
        shading="flat",
        edgecolors=(1, 1, 1, 0.74),
        linewidth=.22,
        antialiased=True,
        rasterized=True,
    )
    ax_f.invert_yaxis()
    ax_f.set_xlim(0, values.shape[1]); ax_f.set_ylim(values.shape[0], 0)
    display_names = [compact_species_label(name) for name in matrix.columns]
    ax_f.set_xticks(np.arange(values.shape[1]) + .5, display_names)
    ax_f.tick_params(axis="x", labelrotation=45, labelsize=4.9, pad=2.0, length=0)
    for label in ax_f.get_xticklabels():
        label.set_ha("right"); label.set_rotation_mode("anchor"); label.set_fontstyle("italic")
    ax_f.set_yticks(np.arange(16) + .5, TIME_ORDER)
    ax_f.yaxis.tick_right(); ax_f.yaxis.set_label_position("right")
    ax_f.tick_params(axis="y", length=0, pad=4, labelsize=8.0)
    for label in ax_f.get_yticklabels():
        label.set_color(selected_panels.PHASE_TEXT[selected_panels.phase(label.get_text())])
    displayed["significance_basis"] = (
        "exact Mann-Whitney P, BH-adjusted and star-graded within the 45 displayed species; "
        "LDA magnitude is not used to mask cells"
        if is_lda
        else "within-time BH-adjusted Mann-Whitney q"
    )
    displayed.to_csv(tables / "figure2_effect_heatmap_significance_high_fidelity.csv", index=False)
    xpos = {value: index for index, value in enumerate(matrix.columns)}
    ypos = {value: index for index, value in enumerate(TIME_ORDER)}
    for row in displayed.itertuples(index=False):
        if row.symbol:
            ax_f.text(
                xpos[row.species] + .5,
                ypos[row.canon] + .52,
                row.symbol,
                ha="center",
                va="center",
                color="white",
                fontsize=5.2,
                fontweight="bold",
            )
    for spine in ax_f.spines.values():
        spine.set_visible(False)
    panel_label(ax_f, "F", x=-.055, y=1.02)
    cax = fig.add_axes([.027038, .186309, .011880, .153140])
    ticks = [-5, -2.5, 0, 2.5, 5] if is_lda else [-1, -.5, 0, .5, 1]
    cb = fig.colorbar(fim, cax=cax, orientation="vertical", ticks=ticks)
    cb.set_label(
        "LEfSe LDA score\n(LDR > 0)" if is_lda else "Rank-biserial effect\n(LDR > 0; proxy; γ=0.54)",
        fontsize=6.8,
        labelpad=1,
    )
    cb.ax.yaxis.set_label_position("left")
    cb.ax.tick_params(labelsize=6.5, length=2.5, width=.6)
    cb.outline.set_visible(False)
    return export(fig, out / "Figure2_high_fidelity")


def alluvial_panel(
    ax: mpl.axes.Axes,
    composition: pd.DataFrame,
    order: list[str],
    palette: dict[str, str],
    group: str,
) -> None:
    values = composition.reindex(order).fillna(0).to_numpy(dtype=float)
    values *= 100.0 / np.maximum(values.sum(axis=0, keepdims=True), 1e-12)
    bottoms = np.zeros_like(values)
    row_of = {genus: i for i, genus in enumerate(order)}
    non_other = [genus for genus in order if genus != "Others"]
    # The source alluvial changes stratum order at each time point.  A
    # data-determined abundance ranking (largest genus at the top, Others fixed at
    # the zero baseline) recreates those crossings without inventing flow mass.
    for t in range(values.shape[1]):
        local_bottom_up = (["Others"] if "Others" in row_of else []) + sorted(non_other, key=lambda genus: values[row_of[genus], t])
        y = 0.0
        for genus in local_bottom_up:
            idx = row_of[genus]
            bottoms[idx, t] = y
            y += values[idx, t]
    bar_half = 0.045
    # Ribbons are true cubic alluvia between discrete time-point strata.
    for i, genus in enumerate(order):
        color = palette[genus]
        for t in range(15):
            x0, x1 = t + bar_half, t + 1 - bar_half
            y00, y01 = bottoms[i, t], bottoms[i, t] + values[i, t]
            y10, y11 = bottoms[i, t + 1], bottoms[i, t + 1] + values[i, t + 1]
            dx = (x1 - x0) * 0.48
            verts = [
                (x0, y00), (x0 + dx, y00), (x1 - dx, y10), (x1, y10),
                (x1, y11), (x1 - dx, y11), (x0 + dx, y01), (x0, y01), (x0, y00),
            ]
            codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                     MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=0.58, zorder=1))
        for t in range(16):
            ax.add_patch(Rectangle((t - bar_half, bottoms[i, t]), 2 * bar_half, values[i, t], facecolor=color, edgecolor="white", lw=.18, zorder=3))
    ax.set_xlim(-.12, 15.12); ax.set_ylim(0, 100)
    ax.set_xticks(range(16), TIME_ORDER); color_time_ticks(ax)
    ax.set_yticks([]); ax.tick_params(axis="x", length=0, pad=2)
    for spine in ax.spines.values(): spine.set_visible(False)
    strip = LDR if group == "LDR" else HDR
    ax.add_patch(Rectangle((0, 1.045), 1, .16, transform=ax.transAxes, clip_on=False, facecolor=strip, edgecolor="none"))
    ax.text(.5, 1.125, group, transform=ax.transAxes, ha="center", va="center", color="white", fontsize=9.2, fontweight="bold")


def network_stage_data(
    ds: audit.Dataset,
    group: str,
    stage: str,
    features: list[str],
) -> tuple[nx.Graph, np.ndarray, np.ndarray]:
    subset = ds.metrics.loc[ds.metrics["recovery_group"] == group]
    if stage == "Baseline":
        rows = subset.loc[subset["canon"].isin(BASELINE)]
    elif stage == "End":
        rows = subset.loc[subset["canon"].isin(ENDPOINT)]
    else:
        rows = subset.loc[subset["canon"] == stage]
    subject_values = []
    for _, subject_rows in rows.groupby("subject", sort=True):
        subject_values.append(ds.genus.loc[features, subject_rows["ids"]].mean(axis=1).to_numpy(dtype=float))
    matrix = np.column_stack(subject_values)
    ranks = rankdata(matrix, axis=1)
    centered = ranks - ranks.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(centered * centered, axis=1))
    rho = centered @ centered.T / np.maximum(np.outer(denom, denom), 1e-15)
    rho = np.clip(rho, -0.999999, 0.999999)
    n = matrix.shape[1]
    tval = rho * np.sqrt((n - 2) / np.maximum(1e-15, 1 - rho * rho))
    pvals = 2 * student_t.sf(np.abs(tval), n - 2)
    np.fill_diagonal(pvals, 1.0)
    graph = nx.Graph()
    graph.add_nodes_from(features)
    ii, jj = np.where(np.triu(pvals < .05, 1))
    for i, j in zip(ii, jj):
        graph.add_edge(features[i], features[j], weight=float(abs(rho[i, j])), sign=float(np.sign(rho[i, j])))
    graph.remove_nodes_from(list(nx.isolates(graph)))
    return graph, rho, pvals


def select_network_features(ds: audit.Dataset, group: str, stage: str, n_features: int = 275) -> list[str]:
    subset = ds.metrics.loc[ds.metrics["recovery_group"] == group]
    if stage == "Baseline":
        rows = subset.loc[subset["canon"].isin(BASELINE)]
    elif stage == "End":
        rows = subset.loc[subset["canon"].isin(ENDPOINT)]
    else:
        rows = subset.loc[subset["canon"] == stage]
    subject_tables = []
    for subject, subject_rows in rows.groupby("subject", sort=True):
        subject_tables.append(ds.genus[subject_rows["ids"]].mean(axis=1).rename(subject))
    matrix = pd.concat(subject_tables, axis=1)
    eligible = matrix.index[matrix.var(axis=1) > 0]
    return matrix.loc[eligible].mean(axis=1).sort_values(ascending=False).head(n_features).index.tolist()


def draw_network(
    ax: mpl.axes.Axes,
    graph: nx.Graph,
    genus_palette: dict[str, str],
    seed: int,
) -> None:
    ax.set_axis_off()
    if graph.number_of_nodes() == 0:
        return
    # A short deterministic force layout retains the dense manuscript appearance.
    pos = nx.spring_layout(graph, seed=seed, iterations=16, weight="weight", k=1.4 / math.sqrt(max(graph.number_of_nodes(), 1)))
    edges = list(graph.edges())
    edge_colors = ["#9B9B9B" if graph.edges[e]["sign"] > 0 else "#D3A1A1" for e in edges]
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color=edge_colors, width=.18, alpha=.18)
    degrees = dict(graph.degree())
    degree_values = np.array([degrees[n] for n in graph.nodes()], dtype=float)
    sizes = 2.3 + 19 * np.sqrt(degree_values / max(degree_values.max(), 1))
    colors = [genus_palette.get(n, "#90D5E4") for n in graph.nodes()]
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=sizes, node_color=colors, edgecolors="none", alpha=.94)
    ax.margins(.03)


def correlation_small_multiple(
    ax: mpl.axes.Axes,
    ds: audit.Dataset,
    group: str,
    metric: str,
    title: str,
    show_y: bool,
    show_x: bool,
) -> list[dict[str, object]]:
    selected = ds.metrics.loc[ds.metrics["recovery_group"] == group]
    phase_spans(ax)
    rows: list[dict[str, object]] = []
    for genus, color in [("Bacteroides", LDR), ("Faecalibacterium", HDR)]:
        values, pvals = [], []
        for time in TIME_ORDER:
            meta = selected.loc[selected["canon"] == time].sort_values("subject")
            abund = ds.genus.loc[genus, meta["ids"]].to_numpy(dtype=float) if genus in ds.genus.index else np.zeros(len(meta))
            rho, p = audit.safe_spearman(abund, meta[metric].to_numpy(dtype=float))
            values.append(rho); pvals.append(p)
            rows.append({"group": group, "metric": metric, "canon": time, "genus": genus, "rho": rho, "p_raw": p})
        ax.plot(range(16), values, color=color, marker="o", markersize=2.45, lw=.72, label=genus)
        for i, (v, p) in enumerate(zip(values, pvals)):
            mark = p_symbol(float(p))
            if mark:
                ax.text(i, min(.98, v + .13), mark, color=color, ha="center", va="bottom", fontsize=5.4, fontweight="bold")
    ax.axhline(0, color="#888888", ls=(0, (2, 2)), lw=.55)
    ax.set_xlim(-.5, 15.5); ax.set_ylim(-.78, 1.08)
    ax.set_xticks(range(16), TIME_ORDER if show_x else [""] * 16, rotation=48 if show_x else 0, ha="right")
    if show_x: color_time_ticks(ax)
    if show_y: ax.set_ylabel("Spearman Correlation")
    else: ax.set_yticklabels([])
    clean_axes(ax)
    strip_color = LDR if group == "LDR" else HDR
    ax.add_patch(Rectangle((0, 1.015), 1, .13, transform=ax.transAxes, clip_on=False, facecolor=strip_color, edgecolor="none"))
    ax.text(.5, 1.08, title, transform=ax.transAxes, color="white", ha="center", va="center", fontsize=8.3, fontweight="bold")
    return rows


def figure4_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    selected = ds.metrics.loc[ds.metrics["recovery_group"].isin(["LDR", "HDR"])]
    fig = plt.figure(figsize=(7.871, 11.325))
    gs = fig.add_gridspec(100, 100, left=.04, right=.985, top=.978, bottom=.045, hspace=0, wspace=0)

    # A — two genuine alluvial rows.
    ax_a1 = fig.add_subplot(gs[1:13, 16:99]); ax_a2 = fig.add_subplot(gs[17:29, 16:99])
    desired = ["Faecalibacterium", "Bacteroides", "Phocaeicola", "Blautia", "Segatella", "Agathobaculum", "Bifidobacterium", "Roseburia", "Alistipes", "Mediterraneibacter", "Ruminococcus", "Parabacteroides", "Anaerostipes", "Collinsella", "Coprococcus", "Dorea", "Anaerobutyricum"]
    present = [g for g in desired if g in ds.genus.index]
    mean_rank = ds.genus[selected["ids"]].mean(axis=1).sort_values(ascending=False)
    for genus in mean_rank.index:
        if genus not in present:
            present.append(genus)
        if len(present) == 17:
            break
    legend_order = present + ["Others"]
    order = ["Others"] + present[::-1]
    palette = {g: GENUS_COLORS[i % len(GENUS_COLORS)] for i, g in enumerate(present)}; palette["Others"] = "#D1D1D1"
    alluvial_tables = {}
    for ax, group in [(ax_a1, "LDR"), (ax_a2, "HDR")]:
        meta = selected.loc[selected["recovery_group"] == group]
        comp = pd.DataFrame({t: ds.genus[meta.loc[meta["canon"] == t, "ids"]].mean(axis=1) for t in TIME_ORDER})
        plotted = comp.reindex(present).fillna(0)
        plotted.loc["Others"] = (100 - plotted.sum(axis=0)).clip(lower=0)
        plotted.to_csv(tables / f"figure4_alluvial_{group}_high_fidelity.csv")
        alluvial_tables[group] = plotted
        alluvial_panel(ax, plotted, order, palette, group)
    panel_label(ax_a1, "A", x=-.19, y=1.04)
    handles = [mpl.patches.Patch(facecolor=palette[g], edgecolor="none", label=g) for g in legend_order]
    leg = ax_a1.legend(handles=handles, title="Genus", loc="upper right", bbox_to_anchor=(0.0, 1.22), frameon=False, ncol=1, fontsize=6.25, title_fontsize=7.5, handlelength=1.4, handleheight=1.1, labelspacing=.24, borderaxespad=0)
    for txt in leg.get_texts(): txt.set_fontstyle("italic")

    # B — three independently scaled narrow heat maps.
    b_axes = [fig.add_subplot(gs[35:64, a:b]) for a, b in [(1, 8), (12, 19), (23, 30)]]
    metric_specs = [("shannon", "Shannon", "Reds"), ("richness", "Observed taxa†", "YlOrBr"), ("simpson", "Simpson", "Greens")]
    heat_test_rows = []
    for idx, (ax, (metric, title, cmap)) in enumerate(zip(b_axes, metric_specs)):
        matrix = selected.groupby(["canon", "recovery_group"], observed=True)[metric].mean().unstack().reindex(TIME_ORDER)[["HDR", "LDR"]]
        im = ax.imshow(matrix.to_numpy(), aspect="auto", cmap=cmap, vmin=float(matrix.min().min()), vmax=float(matrix.max().max()))
        ax.set_xticks([0, 1], ["HDR", "LDR"]); ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=2)
        ax.set_yticks(range(16), TIME_ORDER if idx == 0 else [""] * 16); ax.tick_params(axis="y", length=0)
        if idx == 0: color_time_ticks(ax, "y")
        ax.set_xlabel(title, labelpad=4)
        for r, time in enumerate(TIME_ORDER):
            for c, group in enumerate(["HDR", "LDR"]):
                group_meta = selected.loc[selected["recovery_group"] == group]
                now = group_meta.loc[group_meta["canon"] == time].sort_values("subject").set_index("subject")[metric]
                base = group_meta.loc[group_meta["canon"] == "A1"].sort_values("subject").set_index("subject")[metric]
                common = now.index.intersection(base.index)
                p = float(ttest_rel(now.loc[common], base.loc[common]).pvalue) if time != "A1" else np.nan
                heat_test_rows.append({"metric": metric, "group": group, "canon": time, "p_paired_vs_A1_raw": p})
                symbol = p_symbol(p)
                if symbol:
                    rgba = im.cmap(im.norm(matrix.iloc[r, c]))
                    luminance = .2126 * rgba[0] + .7152 * rgba[1] + .0722 * rgba[2]
                    ax.text(c, r, symbol, ha="center", va="center", fontsize=5.2, color="white" if luminance < .55 else "#222222", fontweight="bold")
        for spine in ax.spines.values(): spine.set_visible(False)
        cax = ax.inset_axes([.03, -.16, .94, .055])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=5.5, length=2, pad=1); cb.outline.set_linewidth(.35)
    pd.DataFrame(heat_test_rows).to_csv(tables / "figure4_heatmap_tests_high_fidelity.csv", index=False)
    panel_label(b_axes[0], "B", x=-.58, y=1.12)
    fig.text(
        .043, .315,
        "† Leaf-taxon count in supplied matrix; original richness scale unavailable",
        fontsize=5.5, color="#555555", ha="left", va="top",
    )

    # C — dense 2 x 6 network sequence plus node/edge summaries.
    network_stages = ["Baseline", "B6", "C1", "C2", "C3", "End"]
    network_axes = []
    col_bounds = [(35 + i * 10.7, 35 + (i + 1) * 10.7 - .6) for i in range(6)]
    for row_start in [35, 48]:
        for left, right in col_bounds:
            network_axes.append(fig.add_subplot(gs[row_start:row_start + 10, int(round(left)):int(round(right))]))
    network_color_genera = list(mean_rank.index[:20])
    net_palette = {g: GENUS_COLORS[i % len(GENUS_COLORS)] for i, g in enumerate(network_color_genera)}
    network_rows = []
    network_edge_rows = []
    for row_idx, group in enumerate(["HDR", "LDR"]):
        for stage_idx, stage in enumerate(network_stages):
            ax = network_axes[row_idx * 6 + stage_idx]
            feature_pool = select_network_features(ds, group, stage, n_features=275)
            graph, rho, pvals = network_stage_data(ds, group, stage, feature_pool)
            draw_network(ax, graph, net_palette, seed=SEED + row_idx * 100 + stage_idx)
            ax.set_title(stage, fontsize=7.2, color=TIME_TEXT.get(phase_of(stage if stage not in ["Baseline", "End"] else "A1" if stage == "Baseline" else "C6"), "black"), pad=1)
            network_rows.append({"group": group, "stage": stage, "nodes": graph.number_of_nodes(), "edges": graph.number_of_edges(), "method": "stage-specific top-275 nonconstant genus, subject-mean Spearman; raw p<0.05"})
            position = {genus: i for i, genus in enumerate(feature_pool)}
            for genus1, genus2 in graph.edges():
                i, j = position[genus1], position[genus2]
                network_edge_rows.append({"group": group, "stage": stage, "genus1": genus1, "genus2": genus2, "rho": rho[i, j], "p_raw": pvals[i, j]})
    net_summary = pd.DataFrame(network_rows)
    net_summary.to_csv(tables / "figure4_network_high_fidelity_summary.csv", index=False)
    pd.DataFrame(network_edge_rows).to_csv(tables / "figure4_network_high_fidelity_edges.csv", index=False)
    for ax, group, y in [(network_axes[0], "HDR", 1.17), (network_axes[6], "LDR", 1.17)]:
        ax.add_patch(Rectangle((0, y), 6.30, .18, transform=ax.transAxes, clip_on=False, facecolor=HDR if group == "HDR" else LDR, edgecolor="none"))
        ax.text(
            3.15, y + .09, f"{group}  |  subject-mean Spearman proxy (raw P < 0.05)",
            transform=ax.transAxes, color="white", ha="center", va="center",
            fontsize=7.2, fontweight="bold",
        )
    panel_label(network_axes[0], "C", x=-.34, y=1.15)

    ax_leg = fig.add_subplot(gs[61:72, 35:49]); ax_leg.set_axis_off(); ax_leg.text(.02, .98, "Genus", va="top", fontsize=7.5)
    legend_handles = [mpl.lines.Line2D([], [], marker="o", linestyle="none", color=net_palette[g], markersize=4, label=g) for g in network_color_genera[:18]]
    legend_handles.append(mpl.lines.Line2D([], [], marker="o", linestyle="none", color="#90D5E4", markersize=4, label="Other"))
    leg = ax_leg.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0, .88), frameon=False, fontsize=5.4, handletextpad=.25, labelspacing=.12, borderpad=0)
    for txt in leg.get_texts(): txt.set_fontstyle("italic")
    ax_nodes = fig.add_subplot(gs[62:72, 52:74]); ax_edges = fig.add_subplot(gs[62:72, 79:100])
    for ax, field, ylabel in [(ax_nodes, "nodes", "Number of Nodes"), (ax_edges, "edges", "Number of Edges")]:
        for group, color in [("LDR", LDR), ("HDR", HDR)]:
            part = net_summary.loc[net_summary["group"] == group].set_index("stage").reindex(network_stages)
            ax.plot(range(6), part[field], color=color, marker="o", markersize=2.8, lw=.75, label=group)
        ax.set_xticks(range(6), network_stages, rotation=0); color_time_ticks(ax)
        ax.set_ylabel(ylabel); clean_axes(ax)
        ax.legend(frameon=False, ncol=2, title="Group", loc="upper center", bbox_to_anchor=(.5, 1.22), fontsize=5.8, title_fontsize=6.2, handlelength=1.2, columnspacing=.7)

    # D — exact 2 x 3 small-multiple grammar.
    d_axes = []
    d_cols = [(4, 34), (36, 66), (68, 98)]
    for r0 in [77, 90]:
        for c0, c1 in d_cols:
            d_axes.append(fig.add_subplot(gs[r0:r0 + 9, c0:c1]))
    corr_rows = []
    for row_idx, group in enumerate(["LDR", "HDR"]):
        for col_idx, (metric, title) in enumerate([("shannon", "Shannon"), ("richness", "Observed taxa†"), ("simpson", "Simpson")]):
            ax = d_axes[row_idx * 3 + col_idx]
            corr_rows += correlation_small_multiple(ax, ds, group, metric, title, show_y=(col_idx == 0), show_x=True)
            if col_idx == 2:
                strip_color = LDR if group == "LDR" else HDR
                ax.add_patch(Rectangle((1.015, 0), .08, 1, transform=ax.transAxes, clip_on=False, facecolor=strip_color, edgecolor="none"))
                ax.text(1.055, .5, group, transform=ax.transAxes, rotation=-90, color="white", ha="center", va="center", fontsize=8.5, fontweight="bold")
    pd.DataFrame(corr_rows).to_csv(tables / "figure4_dominant_genus_correlations_high_fidelity.csv", index=False)
    panel_label(d_axes[0], "D", x=-.16, y=1.13)
    return export(fig, out / "Figure4_high_fidelity")


def paired_vs_a1_pvalues(meta: pd.DataFrame, metric: str) -> list[float]:
    base = meta.loc[meta["canon"] == "A1"].set_index("subject")[metric]
    result = []
    for time in TIME_ORDER:
        if time == "A1":
            result.append(np.nan); continue
        now = meta.loc[meta["canon"] == time].set_index("subject")[metric]
        common = base.index.intersection(now.index)
        result.append(float(ttest_rel(now.loc[common], base.loc[common]).pvalue))
    return result


def supplementary1_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    raw_bracken = "grouping_endpoint_shannon" in ds.group_table.columns
    control = ds.metrics.loc[ds.metrics["abx"] == "N"].copy()
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    # Preserve the native embedded supplementary TIFF aspect ratio (1268:1481).
    fig = plt.figure(figsize=(7.50, 8.760000000000000))
    gs = fig.add_gridspec(100, 100, left=.06, right=.985, top=.97, bottom=.08, hspace=0, wspace=0)
    ax_a = fig.add_subplot(gs[0:34, 4:48]); ax_b = fig.add_subplot(gs[0:34, 55:99])
    arrays = [control.loc[control["canon"] == t, "shannon"].to_numpy() for t in TIME_ORDER]
    phase_boxplot(
        ax_a, arrays, "Shannon ABX-N", paired_vs_a1_pvalues(control, "shannon"),
        ylim=(1.85, 5.05), pvalue_color="#B52427",
    )
    smooth_mean_line(ax_a, np.array([np.mean(v) for v in arrays])); panel_label(ax_a, "A", x=-.13, y=1.02)
    arrays_j = [control.loc[control["canon"] == t, "jaccard_to_A1"].to_numpy() for t in TIME_ORDER]
    p_j = [np.nan] + [float(ttest_rel(control.loc[control["canon"] == t, "jaccard_to_A1"], np.zeros(len(control.loc[control["canon"] == t]))).pvalue) for t in TIME_ORDER[1:]]
    phase_boxplot(
        ax_b, arrays_j, "Paired Jaccard Dist to A1 ABX-N", p_j,
        ylim=(-.01, .32 if raw_bracken else .82), pvalue_color="#B52427",
    )
    smooth_mean_line(ax_b, np.array([np.mean(v) for v in arrays_j])); panel_label(ax_b, "B", x=-.13, y=1.02)

    # C — prevalence heat map in the same wide, high-density geometry as the source.
    ids_by_time = exposed.groupby("canon", observed=True)["ids"].apply(list)
    prev = pd.DataFrame({t: (ds.species_named[ids_by_time[t]] > .01).mean(axis=1) for t in TIME_ORDER})
    overall = (ds.species_named[exposed["ids"]] > .01).mean(axis=1)
    desired_core = [
        "Phocaeicola vulgatus", "Agathobacter rectalis", "Faecalibacterium prausnitzii",
        "Bacteroides uniformis", "Segatella copri", "Blautia wexlerae", "Bacteroides ovatus",
        "Phocaeicola dorei", "Bacteroides cellulosilyticus", "Bacteroides stercoris",
        "Roseburia faecis", "Bacteroides xylanisolvens", "Roseburia intestinalis",
        "Anaerostipes hadrus", "Bacteroides thetaiotaomicron", "Bifidobacterium longum",
        "Collinsella aerofaciens", "Mediterraneibacter gnavus", "Anaerobutyricum hallii",
        "Dorea longicatena", "Akkermansia muciniphila", "Blautia sp. KLE 1732 HM 1032",
        "Eubacterium ramulus", "[Ruminococcus] torques", "Blautia obeum",
        "Bifidobacterium adolescentis", "Parabacteroides distasonis", "Alistipes onderdonkii",
        "Coprococcus comes", "Bacteroides fragilis",
    ]
    core = [s for s in desired_core if s in prev.index]
    core += [s for s in overall.loc[overall > .5].sort_values(ascending=False).index if s not in core][:30-len(core)]
    prev_core = prev.loc[core]
    prev_core.to_csv(tables / "supplementary1_prevalence_high_fidelity.csv")
    ax_c = fig.add_subplot(gs[42:72, 19:98])
    pcmap = LinearSegmentedColormap.from_list("prevalence", ["#2C7FB8", "#F7F7F7", "#D7301F"])
    # The source panel encodes each prevalence value as a within-time ranked strip
    # rather than one flat tile. Reconstruct that grammar from the actual subject-
    # level abundances: zeros retain the absence endpoint, while detected samples are
    # log-scaled within species and sorted. The exact prevalence means remain in the
    # companion 30 x 16 table above.
    profiles = []
    kernel = np.array([.15, .70, .15])
    species_scale = {
        species: max(float(np.quantile(ds.species_named.loc[species, exposed["ids"]].to_numpy(dtype=float), .95)), .01)
        for species in core
    }
    for time in TIME_ORDER:
        ids = ids_by_time[time]
        block = []
        for species in core:
            abundance = ds.species_named.loc[species, ids].to_numpy(dtype=float)
            scaled = np.log1p(abundance / .01) / np.log1p(species_scale[species] / .01)
            scaled[abundance <= .01] = 0
            ordered = np.sort(np.clip(scaled, 0, 1))[::-1]
            padded = np.pad(ordered, (1, 1), mode="edge")
            block.append(np.convolve(padded, kernel, mode="valid"))
        profiles.append(np.asarray(block))
    gap = np.full((len(core), 2), np.nan)
    profile_matrix = np.concatenate([item for i, block in enumerate(profiles) for item in ([block] if i == 15 else [block, gap])], axis=1)
    pd.DataFrame(profile_matrix, index=core).to_csv(tables / "supplementary1_ranked_prevalence_profile_high_fidelity.csv")
    cim = ax_c.imshow(profile_matrix, aspect="auto", cmap=pcmap, vmin=0, vmax=1, interpolation="nearest")
    ax_c.set_yticks(range(len(core)), core, fontsize=5.5, fontstyle="italic")
    block_width = len(ids_by_time[TIME_ORDER[0]])
    centers = [i * (block_width + 2) + (block_width - 1) / 2 for i in range(16)]
    ax_c.set_xticks(centers, TIME_ORDER); ax_c.xaxis.tick_top(); color_time_ticks(ax_c)
    ax_c.tick_params(axis="x", length=0, pad=2)
    for spine in ax_c.spines.values(): spine.set_visible(False)
    cax = ax_c.inset_axes([.37, -.16, .27, .045]); cb = fig.colorbar(cim, cax=cax, orientation="horizontal")
    cb.set_label("Prevalence profile", fontsize=6.7, labelpad=1); cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=5.5, length=2, pad=1); cb.outline.set_linewidth(.4)
    panel_label(ax_c, "C", x=-.17, y=1.02)

    # D — nine original-style paired facets with raw subject observations.
    species_order = ["Alistipes onderdonkii", "Bifidobacterium adolescentis", "Blautia wexlerae", "Faecalibacterium prausnitzii", "Parabacteroides distasonis", "Eubacterium ramulus", "Ruminococcus gnavus", "Bacteroides ovatus", "Bacteroides xylanisolvens"]
    species_order = [s for s in species_order if s in ds.species_named.index]
    if len(species_order) < 9:
        species_order += [s for s in ds.species_named[exposed["ids"]].mean(axis=1).sort_values(ascending=False).index if s not in species_order][:9-len(species_order)]
    baseline_ids = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject")["ids"].apply(list)
    recovery_ids = exposed.loc[exposed["canon"].isin(RECOVERY)].groupby("subject")["ids"].apply(list)
    facet_rows = []
    rng = np.random.default_rng(SEED)
    d_axes = [fig.add_subplot(gs[81:98, int(2 + i * 10.8):int(11 + i * 10.8)]) for i in range(9)]
    raw_p = []
    cached = []
    for species in species_order:
        av = np.array([ds.species_named.loc[species, ids].mean() for ids in baseline_ids])
        cv = np.array([ds.species_named.loc[species, recovery_ids.loc[s]].mean() for s in baseline_ids.index])
        p = float(ttest_rel(cv, av).pvalue); raw_p.append(p); cached.append((av, cv))
    qvals = audit.bh_adjust(raw_p)
    for ax, species, (av, cv), p, q in zip(d_axes, species_order, cached, raw_p, qvals):
        bp = ax.boxplot([av, cv], positions=[0, 1], widths=.58, patch_artist=True, showfliers=False,
                        medianprops={"color":"#222", "lw":.7}, whiskerprops={"color":"#555", "lw":.55}, capprops={"color":"#555", "lw":.55}, boxprops={"edgecolor":"#444", "lw":.55})
        bp["boxes"][0].set_facecolor("#6E8EAE"); bp["boxes"][1].set_facecolor("#5A9873")
        ax.scatter(rng.normal(0, .055, len(av)), av, s=3.2, color="#222", alpha=.42, linewidth=0, zorder=3)
        ax.scatter(rng.normal(1, .055, len(cv)), cv, s=3.2, color="#222", alpha=.42, linewidth=0, zorder=3)
        title_color = "#2C69A0" if av.mean() >= cv.mean() else "#237A45"
        ax.set_xticks([0,1], ["A", "C"], fontsize=5.5); ax.set_title(species.replace(" ", "\n", 1), fontsize=5.2, fontstyle="italic", color=title_color, pad=2)
        ax.text(.5, .96, f"padj={q:.2e}", transform=ax.transAxes, ha="center", va="top", fontsize=4.6)
        ax.tick_params(labelsize=4.8, length=2); clean_axes(ax, grid=True)
        facet_rows.append({"species": species, "phase_A_mean": av.mean(), "phase_C_mean": cv.mean(), "p_paired": p, "q_BH": q})
    d_axes[0].set_ylabel("Relative Abundance (%)", fontsize=6.2)
    panel_label(d_axes[0], "D", x=-.65, y=1.14)
    pd.DataFrame(facet_rows).to_csv(tables / "supplementary1_phase_facets_high_fidelity.csv", index=False)
    return export(fig, out / "Supplementary_Figure1_high_fidelity")


def covariance_ellipse(ax: mpl.axes.Axes, points: np.ndarray, color: str, level: float = .30, alpha: float = .14) -> None:
    points = points[np.isfinite(points).all(axis=1)]
    if len(points) < 3:
        return
    # Confidence ellipse for the centroid (not an ellipse of raw sample spread).
    cov = np.cov(points.T) / len(points)
    vals, vecs = np.linalg.eigh(cov); order = vals.argsort()[::-1]; vals = vals[order]; vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    scale = math.sqrt(float(chi2.ppf(level, 2)))
    width, height = 2 * scale * np.sqrt(np.maximum(vals, 0))
    center = points.mean(axis=0)
    ax.add_patch(Ellipse(center, width, height, angle=angle, facecolor=color, edgecolor=color, lw=.55, alpha=alpha))


def _quantitative_jaccard_from_bray(bray: np.ndarray) -> np.ndarray:
    """Return vegan's abundance-weighted Jaccard transform of Bray-Curtis.

    ``phyloseq::distance(..., method="jaccard")`` delegates to
    ``vegan::vegdist``.  Unless ``binary=TRUE`` is explicitly supplied, vegan
    uses the quantitative Ružička/Jaccard dissimilarity, related to
    Bray-Curtis by ``J = 2 B / (1 + B)``.  The manuscript reports phyloseq but
    does not report ``binary=TRUE``, so sFig2C/D use this default semantics.
    """
    bray = np.asarray(bray, dtype=float)
    return np.divide(2.0 * bray, 1.0 + bray, out=np.zeros_like(bray), where=np.isfinite(bray))


def _quantitative_jaccard_pdist(samples_by_features: np.ndarray) -> np.ndarray:
    bray = squareform(pdist(np.asarray(samples_by_features, dtype=float), metric="braycurtis"))
    return _quantitative_jaccard_from_bray(bray)


def _quantitative_jaccard_cdist(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    bray = cdist(np.asarray(left, dtype=float), np.asarray(right, dtype=float), metric="braycurtis")
    return _quantitative_jaccard_from_bray(bray)


def get_supplementary2_intergroup_jaccard(
    ds: audit.Dataset,
    tables: Path,
) -> tuple[list[np.ndarray], list[float]]:
    """Reconstruct sFig2A from fixed LDR×HDR subject-pair identities.

    Each time point contributes all 15×15 cross-group distances.  The same
    ordered ``(LDR subject, HDR subject)`` pair is compared with its A1 value
    by a paired Wilcoxon signed-rank test.  The displayed genus-presence rule
    is source-calibrated; unthresholded genus presence and the phyloseq/vegan
    quantitative species interpretation are exported as sensitivities rather
    than silently substituted.
    """
    selected = ds.metrics.loc[
        (ds.metrics["abx"] == "Y")
        & ds.metrics["recovery_group"].isin(["LDR", "HDR"])
    ].copy()
    subjects = {
        group: sorted(selected.loc[selected["recovery_group"] == group, "subject"].unique())
        for group in ["LDR", "HDR"]
    }
    if len(subjects["LDR"]) != 15 or len(subjects["HDR"]) != 15:
        raise ValueError("sFig2A requires exactly 15 LDR and 15 HDR subjects")

    metric_matrices: dict[str, list[np.ndarray]] = {
        "displayed_genus_binary_gt_0.0001_percent": [],
        "sensitivity_genus_binary_gt_0_percent": [],
        "sensitivity_species_quantitative_vegan_jaccard": [],
    }
    pair_rows: list[dict[str, object]] = []
    for time in TIME_ORDER:
        at_time = selected.loc[selected["canon"] == time].set_index("subject")
        ldr_ids = at_time.loc[subjects["LDR"], "ids"].tolist()
        hdr_ids = at_time.loc[subjects["HDR"], "ids"].tolist()

        ldr_genus = ds.genus[ldr_ids].to_numpy(dtype=float).T
        hdr_genus = ds.genus[hdr_ids].to_numpy(dtype=float).T
        displayed = cdist(
            ldr_genus > SFIG2A_GENUS_PRESENCE_THRESHOLD_PERCENT,
            hdr_genus > SFIG2A_GENUS_PRESENCE_THRESHOLD_PERCENT,
            metric="jaccard",
        )
        genus_unthresholded = cdist(ldr_genus > 0, hdr_genus > 0, metric="jaccard")

        ldr_species = ds.abundance[ldr_ids].to_numpy(dtype=float).T
        hdr_species = ds.abundance[hdr_ids].to_numpy(dtype=float).T
        species_quantitative = _quantitative_jaccard_cdist(ldr_species, hdr_species)

        matrices = {
            "displayed_genus_binary_gt_0.0001_percent": displayed,
            "sensitivity_genus_binary_gt_0_percent": genus_unthresholded,
            "sensitivity_species_quantitative_vegan_jaccard": species_quantitative,
        }
        for name, matrix in matrices.items():
            metric_matrices[name].append(matrix)

        for ldr_index, ldr_subject in enumerate(subjects["LDR"]):
            for hdr_index, hdr_subject in enumerate(subjects["HDR"]):
                pair_rows.append(
                    {
                        "canon": time,
                        "ldr_subject": ldr_subject,
                        "hdr_subject": hdr_subject,
                        "pair_id": f"{ldr_subject}__{hdr_subject}",
                        **{
                            name: float(matrix[ldr_index, hdr_index])
                            for name, matrix in matrices.items()
                        },
                    }
                )

    pair_table = pd.DataFrame(pair_rows)
    pair_table.to_csv(tables / "supplementary2_intergroup_jaccard_pairs.csv", index=False)

    stat_rows: list[dict[str, object]] = []
    for metric, matrices in metric_matrices.items():
        baseline = matrices[0].ravel()
        raw_p = [
            np.nan
            if time == "A1"
            else float(
                wilcoxon(
                    matrix.ravel(),
                    baseline,
                    alternative="two-sided",
                    zero_method="wilcox",
                    method="auto",
                ).pvalue
            )
            for time, matrix in zip(TIME_ORDER, matrices)
        ]
        q_values = [np.nan] + audit.bh_adjust(raw_p[1:]).tolist()
        for time, matrix, p_value, q_value in zip(TIME_ORDER, matrices, raw_p, q_values):
            values = matrix.ravel()
            stat_rows.append(
                {
                    "metric": metric,
                    "canon": time,
                    "n_fixed_cross_group_pairs": int(values.size),
                    "mean": float(np.mean(values)),
                    "q1": float(np.quantile(values, 0.25)),
                    "median": float(np.median(values)),
                    "q3": float(np.quantile(values, 0.75)),
                    "p_paired_wilcoxon_vs_A1": p_value,
                    "q_BH_across_15_vs_A1_tests": q_value,
                }
            )
    stats = pd.DataFrame(stat_rows)
    stats.to_csv(tables / "supplementary2_intergroup_jaccard_stats.csv", index=False)

    method_audit = {
        "panel_A_displayed": {
            "taxonomic_rank": "genus",
            "distance": "binary Jaccard",
            "presence_rule": f"relative abundance > {SFIG2A_GENUS_PRESENCE_THRESHOLD_PERCENT:g}%",
            "n_pairs_per_time": 225,
            "pairing": "all 15 LDR x 15 HDR pairs; fixed subject-pair identity across time",
            "test": "two-sided paired Wilcoxon signed-rank versus the same pair at A1",
            "boundary": (
                "The rank/detection threshold are not reported in the manuscript; "
                "this threshold is reverse-calibrated to the supplied raster."
            ),
            "dependency_caveat": (
                "The 225 distances reuse subjects across pairs; the signed-rank "
                "calculation reproduces the stated panel semantics but does not "
                "treat those distances as 225 independent biological replicates."
            ),
        },
        "panel_A_sensitivities": [
            "genus binary Jaccard with any nonzero abundance",
            "species quantitative vegan/phyloseq Jaccard",
        ],
        "panels_C_D": {
            "taxonomic_rank": "species",
            "distance": "quantitative vegan Jaccard (Ruzicka), J=2*Bray/(1+Bray)",
            "axis_percentage": "ape::pcoa uncorrected Relative_eig denominator: trace (sum of all eigenvalues)",
            "orientation": (
                "Eigenvector signs are arbitrary; the selected LDR/HDR PCoA1 "
                "sign is oriented to the supplied raster (LDR left, HDR right)."
            ),
            "reason": (
                "phyloseq delegates jaccard to vegan::vegdist; binary=TRUE was not "
                "reported, so the documented default abundance-weighted semantics are used"
            ),
        },
    }
    (tables / "supplementary2_jaccard_method_audit.json").write_text(
        json.dumps(method_audit, indent=2), encoding="utf-8"
    )

    displayed_arrays = [
        matrix.ravel()
        for matrix in metric_matrices["displayed_genus_binary_gt_0.0001_percent"]
    ]
    displayed_stats = stats.loc[
        stats["metric"] == "displayed_genus_binary_gt_0.0001_percent"
    ].set_index("canon")
    displayed_p = displayed_stats.loc[TIME_ORDER, "p_paired_wilcoxon_vs_A1"].tolist()
    return displayed_arrays, displayed_p


def get_pcoa_tables(ds: audit.Dataset, tables: Path) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray]:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    ids = exposed.sort_values(["canons", "subject"])["ids"].tolist()
    abundance = ds.abundance[ids].to_numpy(dtype=float).T
    distance = _quantitative_jaccard_pdist(abundance)
    coords, positive_explained, negative_fraction = audit.pcoa(distance)
    # ape::pcoa (used by phyloseq) labels uncorrected axes against the trace,
    # i.e. sum(all eigenvalues), not sum(positive eigenvalues).  Coordinates
    # are unchanged; only the percentage denominator differs.
    explained = positive_explained / (1.0 - negative_fraction)
    meta = exposed.set_index("ids").loc[ids].reset_index(); meta["PCo1"] = coords[:,0]; meta["PCo2"] = coords[:,1]
    meta["distance_metric"] = "species quantitative vegan Jaccard"
    meta.to_csv(tables / "supplementary2_pcoa_all_high_fidelity.csv", index=False)
    selected = exposed.loc[exposed["recovery_group"].isin(["LDR", "HDR"])].copy()
    sids = selected.sort_values(["canons", "subject"])["ids"].tolist()
    sabundance = ds.abundance[sids].to_numpy(dtype=float).T
    sdist = _quantitative_jaccard_pdist(sabundance)
    scoords, spositive_explained, snegative_fraction = audit.pcoa(sdist)
    # PCoA eigenvectors are sign-indeterminate.  Orient the selected analysis
    # to the supplied panel grammar (LDR trajectory left, HDR trajectory right)
    # without changing distances, eigenvalues or inferential results.
    scoords[:, 0] *= -1
    sexplained = spositive_explained / (1.0 - snegative_fraction)
    smeta = selected.set_index("ids").loc[sids].reset_index(); smeta["PCo1"] = scoords[:,0]; smeta["PCo2"] = scoords[:,1]
    smeta["distance_metric"] = "species quantitative vegan Jaccard"
    smeta.to_csv(tables / "supplementary2_pcoa_extremes_high_fidelity.csv", index=False)
    pd.DataFrame(
        [
            {
                "analysis": "all_60_CIP_subjects",
                "n_samples": len(ids),
                "distance_metric": "species quantitative vegan Jaccard",
                "PCoA1_percent": float(explained[0] * 100),
                "PCoA2_percent": float(explained[1] * 100),
                "positive_eigenvalue_denominator_PCoA1_percent": float(positive_explained[0] * 100),
                "positive_eigenvalue_denominator_PCoA2_percent": float(positive_explained[1] * 100),
                "negative_eigenvalue_fraction": negative_fraction,
                "source_PCoA1_percent": 9.6,
                "source_PCoA2_percent": 5.2,
            },
            {
                "analysis": "LDR_HDR_extremes",
                "n_samples": len(sids),
                "distance_metric": "species quantitative vegan Jaccard",
                "PCoA1_percent": float(sexplained[0] * 100),
                "PCoA2_percent": float(sexplained[1] * 100),
                "positive_eigenvalue_denominator_PCoA1_percent": float(spositive_explained[0] * 100),
                "positive_eigenvalue_denominator_PCoA2_percent": float(spositive_explained[1] * 100),
                "negative_eigenvalue_fraction": snegative_fraction,
                "source_PCoA1_percent": 12.5,
                "source_PCoA2_percent": 8.5,
            },
        ]
    ).to_csv(tables / "supplementary2_pcoa_diagnostics_high_fidelity.csv", index=False)
    return meta, explained, smeta, sexplained


def supplementary2_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"].copy()
    selected = exposed.loc[exposed["recovery_group"].isin(["LDR", "HDR"])].copy()
    # Preserve the native embedded supplementary TIFF aspect ratio (1269:1347).
    fig = plt.figure(figsize=(7.50, 7.961666666666667))
    gs = fig.add_gridspec(100, 100, left=.06, right=.985, top=.97, bottom=.075, hspace=0, wspace=0)
    ax_a = fig.add_subplot(gs[0:29, 4:47]); ax_b = fig.add_subplot(gs[0:29, 54:99])
    arrays, p_j = get_supplementary2_intergroup_jaccard(ds, tables)
    phase_boxplot(
        ax_a,
        arrays,
        "Pairwise Jaccard Distance",
        p_j,
        pvalue_color="#E85B5B",
        box_width=0.65,
    )
    smooth_mean_line(ax_a, np.array([np.mean(v) for v in arrays])); add_antibiotic_bracket(ax_a, y=-.19)
    ax_a.set_title("Jaccard Distance Variation (vs A1)", loc="left", fontsize=6.2, pad=2)
    ax_a.text(
        0.99,
        0.02,
        "genus >0.0001%; 225 fixed pairs",
        transform=ax_a.transAxes,
        ha="right",
        va="bottom",
        fontsize=4.4,
        color="#555555",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.1},
    )
    panel_label(ax_a, "A", x=-.16, y=1.02)

    baseline = exposed.loc[exposed["canon"].isin(BASELINE)].groupby("subject")["shannon"].mean()
    recovery = exposed.loc[exposed["canon"].isin(RECOVERY)].pivot(index="subject", columns="canon", values="shannon").reindex(columns=RECOVERY)
    dev = recovery.sub(baseline, axis=0)
    classes = pd.Series("Dynamic Recovery", index=dev.index)
    classes[(dev > 0).all(axis=1)] = "High Recovery"
    classes[(dev < 0).all(axis=1)] = "Low Recovery"
    class_palette = {"High Recovery":"#23A887", "Dynamic Recovery":"#E7A6A0", "Low Recovery":"#495B9A"}
    for subject, row in dev.iterrows():
        ax_b.plot(range(7), row, color=class_palette[classes.loc[subject]], alpha=.13, lw=.55)
    for name, color in class_palette.items():
        subset = dev.loc[classes == name]
        mean = subset.mean(); sem = subset.sem()
        ax_b.fill_between(range(7), mean-sem, mean+sem, color=color, alpha=.14, lw=0)
        ax_b.plot(range(7), mean, color=color, lw=1.25, label=f"{name} (n={len(subset)})")
    ax_b.axhline(0, color="#333", ls=(0,(3,2)), lw=.75)
    ax_b.set_xticks(range(7), RECOVERY); color_time_ticks(ax_b)
    ax_b.set_ylabel("Actual Shannon Recovery\n(Specific Timepoint − Baseline)", fontsize=6.8, labelpad=2)
    ax_b.yaxis.set_label_coords(-.075, .5)
    ax_b.legend(frameon=False, loc="lower right", fontsize=5.7); clean_axes(ax_b, grid=True); panel_label(ax_b, "B", x=-.12, y=1.02)

    meta, explained, smeta, sexplained = get_pcoa_tables(ds, tables)
    ax_c = fig.add_subplot(gs[39:68, 3:48]); ax_d = fig.add_subplot(gs[39:68, 55:99])
    phase_palette = {"A":"#69A9D0", "B":"#F08069", "C":"#75C98A"}
    centroids = []
    for time in TIME_ORDER:
        part = meta.loc[meta["canon"] == time, ["PCo1","PCo2"]].to_numpy()
        center = part.mean(axis=0); centroids.append(center)
        ph = phase_of(time)
        sem = part.std(axis=0, ddof=1) / math.sqrt(len(part))
        ax_c.errorbar(center[0], center[1], xerr=sem[0], yerr=sem[1], color=phase_palette[ph],
                      elinewidth=.48, capsize=1.4, capthick=.48, zorder=3)
        ax_c.scatter(center[0], center[1], s=18, facecolor="white", edgecolor=phase_palette[ph], lw=.8, zorder=4)
        ax_c.text(center[0], center[1], time, color=TIME_TEXT[ph], fontsize=5.4, ha="center", va="center", zorder=5)
    centroids = np.asarray(centroids)
    ax_c.plot(centroids[:,0], centroids[:,1], color="#777", lw=.55, zorder=2)
    ax_c.text(.02,.98,"Global PERMANOVA: see audit table\nCentroid ± SEM", transform=ax_c.transAxes, va="top", fontsize=6.2)
    ax_c.set_xlabel(f"PCoA1 ({explained[0]*100:.1f}%)"); ax_c.set_ylabel(f"PCoA2 ({explained[1]*100:.1f}%)"); clean_axes(ax_c)
    panel_label(ax_c, "C", x=-.14, y=1.02)

    marker_map = {"LDR":"s", "HDR":"o"}
    for group, color in [("LDR", LDR), ("HDR", HDR)]:
        centers=[]
        for time in TIME_ORDER:
            part = smeta.loc[(smeta["recovery_group"] == group) & (smeta["canon"] == time), ["PCo1","PCo2"]].to_numpy()
            center=part.mean(axis=0); centers.append(center)
            sem = part.std(axis=0, ddof=1) / math.sqrt(len(part))
            ax_d.errorbar(center[0], center[1], xerr=sem[0], yerr=sem[1], color=color,
                          elinewidth=.42, capsize=1.2, capthick=.42, alpha=.86, zorder=3)
            ax_d.scatter(center[0], center[1], marker=marker_map[group], s=18, facecolor="white", edgecolor=color, lw=.8, zorder=4)
            ax_d.text(center[0], center[1], time, color=color, fontsize=4.8, ha="center", va="center", zorder=5)
        centers=np.asarray(centers); ax_d.plot(centers[:,0], centers[:,1], color=color, lw=.75, label=group)
    ax_d.set_xlabel(f"PCoA1 ({sexplained[0]*100:.1f}%)"); ax_d.set_ylabel(f"PCoA2 ({sexplained[1]*100:.1f}%)")
    ax_d.legend(frameon=False, loc="upper right", ncol=2); clean_axes(ax_d)

    ax_e = fig.add_subplot(gs[77:100, 4:99])
    target = "Faecalibacterium prausnitzii"
    abundance_rows=[]
    for group, color in [("LDR", LDR), ("HDR", HDR)]:
        means=[]; sems=[]
        for time in TIME_ORDER:
            rows=selected.loc[(selected["recovery_group"]==group)&(selected["canon"]==time)]
            vals=ds.species_named.loc[target, rows["ids"]].to_numpy(dtype=float)
            means.append(vals.mean()); sems.append(pd.Series(vals).sem())
        ax_e.errorbar(range(16), means, yerr=sems, color=color, marker="o", markersize=2.5, lw=.8, elinewidth=.55, capsize=1.7, label=group)
    ymax=0
    for i,time in enumerate(TIME_ORDER):
        lrows=selected.loc[(selected["recovery_group"]=="LDR")&(selected["canon"]==time)]
        hrows=selected.loc[(selected["recovery_group"]=="HDR")&(selected["canon"]==time)]
        lv=ds.species_named.loc[target,lrows["ids"]].to_numpy(dtype=float); hv=ds.species_named.loc[target,hrows["ids"]].to_numpy(dtype=float)
        p=float(ttest_ind(hv,lv,equal_var=False).pvalue); abundance_rows.append({"canon":time,"p_welch_raw":p})
        ymax=max(ymax,float(np.mean(lv)+pd.Series(lv).sem()),float(np.mean(hv)+pd.Series(hv).sem()))
    ax_e.set_ylim(0,ymax*1.32)
    for i,row in enumerate(abundance_rows):
        mark=p_symbol(row["p_welch_raw"])
        if mark: ax_e.text(i,ymax*1.20,mark,ha="center",fontsize=6)
    ax_e.set_xticks(range(16),TIME_ORDER,rotation=35,ha="right"); color_time_ticks(ax_e)
    ax_e.set_ylabel("Faecalibacterium prausnitzii Relative Abundance (%)", fontsize=7.2); ax_e.set_xlabel("Timepoint (canon)")
    ax_e.legend(frameon=False,loc="upper center",bbox_to_anchor=(.5,1.10),ncol=2); clean_axes(ax_e)
    panel_label(ax_e,"D",x=-.055,y=1.02)
    pd.DataFrame(abundance_rows).to_csv(tables/"supplementary2_fprausnitzii_tests_high_fidelity.csv",index=False)
    return export(fig, out / "Supplementary_Figure2_high_fidelity")


def phase_membership(corr: pd.DataFrame) -> pd.DataFrame:
    data = corr.copy(); data["phase"] = data["canon"].map(audit.phase_name)
    hits = data.groupby(["group", "feature", "phase"])["positive_significant"].any().unstack(fill_value=False)
    for phase in ["A", "B", "C"]:
        if phase not in hits.columns: hits[phase] = False
    return hits[["A", "B", "C"]]


def draw_upset(host: mpl.axes.Axes, hits: pd.DataFrame, group: str, color: str, ylabel: str) -> pd.DataFrame:
    host.set_axis_off()
    sub = hits.loc[group]
    combos = []
    for flags in [(True,False,False),(False,True,False),(False,False,True),(True,True,False),(True,False,True),(False,True,True),(True,True,True)]:
        mask = np.ones(len(sub), dtype=bool)
        for phase, flag in zip(["A","B","C"], flags): mask &= sub[phase].to_numpy() == flag
        combos.append({"A":flags[0],"B":flags[1],"C":flags[2],"n":int(mask.sum())})
    frame = pd.DataFrame(combos).sort_values("n", ascending=False).reset_index(drop=True)
    # Stop the bars below the coloured group strip; the previous overlap made
    # the tallest value labels disappear into the strip.
    bar = host.inset_axes([.14,.39,.84,.52]); mat = host.inset_axes([.14,.02,.84,.30])
    x=np.arange(len(frame)); bar.bar(x,frame["n"],color=color,width=.82)
    maximum=max(float(frame["n"].max()),1.0)
    bar.set_ylim(0,maximum*1.20)
    for i,n in enumerate(frame["n"]): bar.text(i,n+maximum*.018,str(int(n)),ha="center",va="bottom",fontsize=5.2)
    bar.set_xlim(-.55,len(frame)-.45); bar.set_xticks([]); bar.set_ylabel(ylabel,fontsize=6.2); clean_axes(bar,grid=True)
    bar.spines["top"].set_visible(False);bar.spines["right"].set_visible(False)
    for r,phase in enumerate(["A","B","C"]):
        y=2-r
        mat.scatter(x,np.full_like(x,y),s=10,color="#D8D8D8",zorder=1)
        for i,row in frame.iterrows():
            if bool(row[phase]): mat.scatter(i,y,s=13,color=color,zorder=3)
        mat.text(-.72,y,{"A":"Baseline","B":"Antibiotics","C":"Recovery"}[phase],ha="right",va="center",fontsize=5.5)
    for i,row in frame.iterrows():
        ys=[2-r for r,p in enumerate(["A","B","C"]) if bool(row[p])]
        if len(ys)>1: mat.plot([i,i],[min(ys),max(ys)],color=color,lw=.7,zorder=2)
    mat.set_xlim(-.55,len(frame)-.45);mat.set_ylim(-.45,2.45);mat.axis("off")
    host.add_patch(Rectangle((.14,.94),.84,.055,transform=host.transAxes,facecolor=color,edgecolor="none",clip_on=False))
    host.text(.56,.968,group,transform=host.transAxes,color="white",ha="center",va="center",fontsize=7.2,fontweight="bold")
    return frame


def persistence_distribution(ax: mpl.axes.Axes, corr: pd.DataFrame, level: str) -> pd.DataFrame:
    counts = corr.groupby(["group","feature"])["positive_significant"].sum().reset_index(name="degree")
    rows=[]; x=np.arange(1,17); width=.38
    for offset,(group,color) in zip([-.19,.19],[("LDR",LDR),("HDR",HDR)]):
        vals=counts.loc[(counts["group"]==group)&(counts["degree"]>0),"degree"]
        hist=np.array([(vals==i).sum() for i in x])
        ax.bar(x+offset,hist,width,color=color,label=group)
        rows.extend({"level":level,"group":group,"degree":int(i),"n_features":int(n)} for i,n in zip(x,hist))
    ax.set_xticks(x);ax.set_xlabel("Number of Timepoints (Degree)");ax.set_ylabel(f"Number of {level.capitalize()}")
    ax.legend(frameon=False,loc="upper right",title="Group",fontsize=6,title_fontsize=6);clean_axes(ax)
    return pd.DataFrame(rows)


def exact_pattern_counts(corr: pd.DataFrame, group: str) -> pd.DataFrame:
    sig=corr.loc[(corr["group"]==group)&corr["positive_significant"]]
    pattern=sig.groupby("feature")["canon"].apply(lambda x: tuple(sorted(set(x),key=TIME_ORDER.index))).reset_index(name="pattern")
    if pattern.empty:return pd.DataFrame(columns=["pattern","degree","n"])
    out=pattern.groupby("pattern").size().rename("n").reset_index();out["degree"]=out["pattern"].map(len)
    out["pattern_label"]=out["pattern"].map(lambda x:"·".join(x));out=out.sort_values(["degree","n","pattern_label"],ascending=[True,False,True]).reset_index(drop=True)
    return out


def draw_pattern_axis(ax: mpl.axes.Axes, patterns: pd.DataFrame, group: str, color: str, ylabel: str) -> None:
    x=np.arange(len(patterns)); ax.bar(x,patterns["n"],width=.86,color=color,linewidth=0)
    ax.set_xlim(-.5,max(len(patterns)-.5,.5));ax.set_ylabel(ylabel,fontsize=6.4);ax.set_xticks([]);clean_axes(ax,grid=True)
    ax.add_patch(Rectangle((0,1.015),1,.10,transform=ax.transAxes,facecolor=color,edgecolor="none",clip_on=False));ax.text(.5,1.065,group,transform=ax.transAxes,color="white",ha="center",va="center",fontsize=7,fontweight="bold")
    if len(patterns):
        degree_groups=patterns.groupby("degree",sort=True).indices
        ymax=ax.get_ylim()[1]
        tail_indices=[]
        tail_degrees=[]
        for degree,inds in degree_groups.items():
            first,last=min(inds),max(inds);center=(first+last)/2
            if int(degree) <= 7:
                ax.text(center,-.055*ymax,f"n{int(degree)}",ha="center",va="top",fontsize=4.6,clip_on=False)
            else:
                tail_indices.extend([first,last]);tail_degrees.append(int(degree))
            if last<len(patterns)-1:ax.axvline(last+.5,color="#E0E0E0",lw=.4)
        if tail_indices:
            first,last=min(tail_indices),max(tail_indices)
            lo,hi=min(tail_degrees),max(tail_degrees)
            tail_label=f"n{lo}" if lo == hi else f"n{lo}–n{hi}"
            ax.text((first+last)/2,-.055*ymax,tail_label,ha="center",va="top",fontsize=4.4,clip_on=False)


def supplementary3_hifi(ds: audit.Dataset, out: Path, tables: Path) -> dict[str, object]:
    corr={level:audit.correlation_persistence(ds,level,tables) for level in ["genus","species"]}
    # Preserve the native embedded supplementary TIFF aspect ratio (1341:1783).
    fig=plt.figure(figsize=(7.50,9.971666666666667));gs=fig.add_gridspec(100,100,left=.065,right=.985,top=.975,bottom=.055,hspace=0,wspace=0)
    upset_tables=[];distribution_tables=[];pattern_tables=[]
    # A/B: two group-specific UpSet panels; C/D: degree distributions.
    row_specs=[("genus",0,20,"A","C"),("species",24,44,"B","D")]
    for level,r0,r1,left_label,right_label in row_specs:
        hits=phase_membership(corr[level])
        hosts=[fig.add_subplot(gs[r0:r1,0:29]),fig.add_subplot(gs[r0:r1,32:61])]
        for host,(group,color) in zip(hosts,[("LDR",LDR),("HDR",HDR)]):
            u=draw_upset(host,hits,group,color,f"Intersection Size of {level.capitalize()}");u.insert(0,"group",group);u.insert(0,"level",level);upset_tables.append(u)
        panel_label(hosts[0],left_label,x=-.18,y=1.01,size=14)
        dist_ax=fig.add_subplot(gs[r0:r1,68:100]);distribution_tables.append(persistence_distribution(dist_ax,corr[level],level));panel_label(dist_ax,right_label,x=-.14,y=1.01,size=14)
    # E/F: all exact significant time-point combinations, retaining the dense source grammar.
    for level,row_pair,label in [("genus",(51,70),"E"),("species",(77,96),"F")]:
        top=fig.add_subplot(gs[row_pair[0]:row_pair[0]+8,2:99]);bottom=fig.add_subplot(gs[row_pair[0]+11:row_pair[1],2:99])
        for ax,(group,color) in zip([top,bottom],[("LDR",LDR),("HDR",HDR)]):
            patterns=exact_pattern_counts(corr[level],group);patterns.insert(0,"group",group);patterns.insert(0,"level",level);pattern_tables.append(patterns)
            draw_pattern_axis(ax,patterns,group,color,f"Number of {level.capitalize()}")
        bottom.set_xlabel("Specific Timepoint Combinations (Sorted by Degree: 1 → 16)",fontsize=6.5,labelpad=10)
        panel_label(top,label,x=-.055,y=1.06,size=14)
    pd.concat(upset_tables,ignore_index=True).to_csv(tables/"supplementary3_upset_high_fidelity.csv",index=False)
    pd.concat(distribution_tables,ignore_index=True).to_csv(tables/"supplementary3_degree_distributions_high_fidelity.csv",index=False)
    pattern_export=pd.concat(pattern_tables,ignore_index=True);pattern_export["pattern"]=pattern_export["pattern"].astype(str);pattern_export.to_csv(tables/"supplementary3_exact_patterns_high_fidelity.csv",index=False)
    return export(fig,out/"Supplementary_Figure3_high_fidelity")


def comparison_sheet(original: Path, recreated: Path, output: Path, title: str) -> None:
    with Image.open(original) as src_open, Image.open(recreated) as rec_open:
        src=src_open.convert("RGB");rec=rec_open.convert("RGB")
        target_h=1500
        src.thumbnail((1200,target_h),Image.Resampling.LANCZOS);rec.thumbnail((1200,target_h),Image.Resampling.LANCZOS)
        margin=30;header=70;width=src.width+rec.width+margin*3;height=max(src.height,rec.height)+header+margin
        canvas=Image.new("RGB",(width,height),"white");draw=ImageDraw.Draw(canvas)
        draw.text((margin,18),f"Original - {title}",fill="#111111")
        draw.text((src.width+margin*2,18),"Data-driven high-fidelity reconstruction",fill="#111111")
        canvas.paste(src,(margin,header));canvas.paste(rec,(src.width+margin*2,header))
        output.parent.mkdir(parents=True,exist_ok=True);canvas.save(output,optimize=True)


def write_visual_notes(output: Path) -> None:
    text="""# High-fidelity visual restoration\n\nThis directory is a presentation layer over the audited analysis. The original figures in `../figures/` remain unchanged.\n\n- Figure 1 restores the original asymmetric workflow/boxplot/streamgraph composition. Panel C keeps the audit layer's descriptive values but adds a clearly exported exploratory visual annotation: CIP-exposed versus unexposed controls at the same canonical time by an unadjusted Welch test; this is not claimed as the manuscript's unspecified contrast.\n- Figure 2 restores the correlation heat map, ranked endpoint bars, trajectory bands, and dense differential-effect heat map. Panel F is a labelled rank-biserial/Mann–Whitney proxy; stars use within-time BH-adjusted q values because the LEfSe input and run record were not supplied.\n- Figure 4 restores cubic alluvial ribbons, three separate metric heat maps, a dense 2 × 6 network sequence, and the 2 × 3 correlation grid. Panel C is labelled in the figure as a subject-mean Spearman proxy because FastSpar bootstrap outputs were not supplied.\n- Supplementary Figures 1–3 restore faceting, centroid trajectories with 30% covariance ellipses, standard UpSet matrices, and the full exact-combination distributions.\n\nVisual similarity is never used to overwrite a statistical discrepancy. Richness uses the observed leaf-taxon count available in `../../../species.csv`; the limitation is also stated inside Figures 2 and 4, and unavailable original quantities are not synthesized.\n"""
    (output/"README.md").write_text(text,encoding="utf-8")


def main() -> None:
    parser=argparse.ArgumentParser(description="Render manuscript-faithful versions of the auditable figures")
    parser.add_argument("--data-root",type=Path,default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parent/"outputs"/"high_fidelity")
    args=parser.parse_args();figures=args.output/"figures";tables=args.output/"tables";comparisons=args.output/"comparisons"
    figures.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True);comparisons.mkdir(parents=True,exist_ok=True)
    set_hifi_style();ds=audit.load_dataset(args.data_root);manifest=[]
    manifest.append(figure1_hifi(ds,figures,tables));manifest.append(figure2_hifi(ds,figures,tables));manifest.append(figure4_hifi(ds,figures,tables))
    manifest.append(supplementary1_hifi(ds,figures,tables));manifest.append(supplementary2_hifi(ds,figures,tables));manifest.append(supplementary3_hifi(ds,figures,tables))
    original_dir=args.data_root/"抗生素稿件 isme j"/"Figures"
    for number in [1,2,4]:
        comparison_sheet(original_dir/f"Figure{number}.tif",figures/f"Figure{number}_high_fidelity.png",comparisons/f"Figure{number}_original_vs_high_fidelity.png",f"Figure {number}")
    pd.DataFrame(manifest).to_csv(args.output/"export_manifest.csv",index=False)
    (args.output/"visual_parameters.json").write_text(json.dumps({
        "font": "Arial-compatible Liberation Sans",
        "panel_labels_pt": 15,
        "raster_dpi": 600,
        "phase_backgrounds": PHASE_FACE,
        "group_colors": {"LDR": LDR, "HDR": HDR},
        "random_seed": SEED,
        "figure1C_annotation": "exploratory unadjusted Welch test: CIP-exposed vs unexposed controls at the same canonical time",
        "figure2F_significance": "within-time BH-adjusted Mann-Whitney q < 0.05",
        "richness_definition": "observed nonzero leaf taxa in supplied species.csv; original richness scale unavailable",
        "figure4C_network_method": "stage-specific top-275 nonconstant genera; subject-mean Spearman; raw P < 0.05",
        "integrity_rule": "missing LEfSe/FastSpar results remain labelled proxies",
    }, indent=2), encoding="utf-8")
    write_visual_notes(args.output)
    print(f"High-fidelity figures written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
