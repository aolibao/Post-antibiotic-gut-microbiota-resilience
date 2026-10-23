#!/usr/bin/env python3
"""Standalone, manuscript-styled restoration of Figure 1D and Figure 2F.

The figures use only the supplied abundance table and metadata.  Figure 2F is
explicitly a Mann–Whitney/rank-biserial proxy because the manuscript LEfSe run,
feature table, and LDA output were not supplied.  Final figure formats are PDF
and RGB LZW-compressed TIFF only.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_all as audit


TIME_ORDER = audit.TIME_ORDER
SEED = audit.SEED
PHASE_TEXT = {"A": "#194B80", "B": "#861E23", "C": "#0D582B"}

# Exact low-to-high ColorBrewer RdBu anchors reverse-engineered from the
# manuscript Figure 2F color bar.  The source raster uses the same anchor
# family with approximately 100 interpolated bands.
FIG2F_RDBU_ANCHORS = [
    "#053061",
    "#2166AC",
    "#4393C3",
    "#92C5DE",
    "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7",
    "#F4A582",
    "#D6604D",
    "#B2182B",
    "#67001F",
]
FIG2F_COLOR_GAMMA = 0.54

# Legend-order palette inferred from the manuscript Figure 1D.  The color blocks
# progress blue → green → red → orange → purple → ochre; order is
# mapped to the supplied table's actual top-30 species, never to invented values.
FIG1D_COLORS = [
    "#A6CEE3", "#74AED1", "#438EC0", "#2D82AF", "#63A8A0",
    "#99CD91", "#98D277", "#6ABB54", "#3BA432", "#6E9E4C",
    "#B89B74", "#FA9594", "#F16667", "#E83739", "#E62F27",
    "#F06C45", "#F9A963", "#FDB055", "#FE982C", "#FE8103",
    "#ED8F47", "#D9A295", "#C3AAD2", "#A07FBC", "#7D54A5",
    "#825D99", "#B9A499", "#F0EB99", "#EAD27A", "#CD9551",
]


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def manuscript_fig2f_cmap(
    gamma: float = FIG2F_COLOR_GAMMA,
) -> ListedColormap:
    """Return the source-hued, saturation-calibrated Fig2F color ramp.

    The numerical normalization remains the bounded rank-biserial scale
    [-1, 0, 1].  Only the lookup ramp is warped symmetrically around the exact
    neutral midpoint.  Gamma was fitted to the 720 source-cell color positions
    because the unavailable LDA scores occupy the source color scale more
    strongly than the supplied rank-biserial proxy values.
    """
    if not 0 < gamma <= 1:
        raise ValueError("Fig2F color gamma must be in (0, 1].")
    source_lut = LinearSegmentedColormap.from_list(
        "fig2f_source_rd_bu_100",
        FIG2F_RDBU_ANCHORS,
        N=100,
    )
    positions = np.linspace(0.0, 1.0, 2049)
    centered = 2.0 * positions - 1.0
    source_positions = 0.5 * (
        1.0 + np.sign(centered) * np.abs(centered) ** gamma
    )
    cmap = ListedColormap(
        source_lut(source_positions),
        name=f"fig2f_source_calibrated_rd_bu_gamma_{gamma:g}",
    )
    return cmap.with_extremes(
        bad="#777777",
        under=FIG2F_RDBU_ANCHORS[0],
        over=FIG2F_RDBU_ANCHORS[-1],
    )


def phase(time: str) -> str:
    return time[0]


def p_symbol(p: float, q: float) -> str:
    """Gate annotations by BH q, then grade stars by the manuscript's raw-P rule."""
    if not np.isfinite(q) or q >= 0.05 or not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    return "*"


def export_pdf_tiff(fig: mpl.figure.Figure, stem: Path, dpi: int = 600) -> dict[str, object]:
    """Export a fixed-size vector PDF and opaque RGB/LZW TIFF."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    tif = stem.with_suffix(".tif")
    fig.savefig(pdf, format="pdf", dpi=dpi, facecolor="white", transparent=False)

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor="white", transparent=False)
    buffer.seek(0)
    with Image.open(buffer) as rendered:
        rgb = rendered.convert("RGB")
        rgb.save(tif, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
        pixels = rgb.size
    inches = tuple(float(v) for v in fig.get_size_inches())
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


def figure1d_table(ds: audit.Dataset) -> tuple[pd.DataFrame, list[str]]:
    exposed = ds.metrics.loc[ds.metrics["abx"] == "Y"]
    ids_by_time = exposed.groupby("canon", observed=True)["ids"].apply(list)
    means = pd.DataFrame(
        {time: ds.species_named[ids_by_time[time]].mean(axis=1) for time in TIME_ORDER}
    )
    legend_order = means.mean(axis=1).nlargest(30).index.tolist()
    selected = means.loc[legend_order].copy()
    selected.loc["Others"] = (100.0 - selected.sum(axis=0)).clip(lower=0)
    stack_order = ["Others", *legend_order[::-1]]
    return selected.loc[stack_order], legend_order


def draw_fig1d(ds: audit.Dataset, output: Path, tables: Path) -> dict[str, object]:
    table, legend_order = figure1d_table(ds)
    table.to_csv(tables / "Fig1D_top30_stream_table.csv")

    # Page and axes proportions were measured from the supplied 1000-dpi TIFF.
    fig = plt.figure(figsize=(4.70, 5.48))
    ax = fig.add_axes([0.140, 0.304, 0.815, 0.649])
    x = np.arange(len(TIME_ORDER), dtype=float)

    stack_order = table.index.tolist()
    color_map = {name: color for name, color in zip(legend_order, FIG1D_COLORS)}
    color_map["Others"] = "#D3D3D3"
    exact = table.to_numpy(dtype=float)
    exact *= 100.0 / np.maximum(exact.sum(axis=0, keepdims=True), 1e-12)
    cumulative = np.vstack([np.zeros(exact.shape[1]), np.cumsum(exact, axis=0)])
    stratum_width = 1.0 / 3.0

    # Manuscript grammar: quintic smoothstep connectors have horizontal first
    # and second derivatives at both strata.  They are deliberately translucent.
    for interval in range(len(TIME_ORDER) - 1):
        xx = np.linspace(interval + stratum_width / 2, interval + 1 - stratum_width / 2, 80)
        t = np.linspace(0.0, 1.0, 80)
        ease = 6 * t**5 - 15 * t**4 + 10 * t**3
        boundaries = (
            cumulative[:, interval, None] * (1.0 - ease)[None, :]
            + cumulative[:, interval + 1, None] * ease[None, :]
        )
        for row, name in enumerate(stack_order):
            ax.fill_between(
                xx,
                boundaries[row],
                boundaries[row + 1],
                color=color_map[name],
                alpha=0.60,
                linewidth=0,
                zorder=1,
            )

    for column, xpos in enumerate(x):
        for row, name in enumerate(stack_order):
            ax.add_patch(
                Rectangle(
                    (xpos - stratum_width / 2, cumulative[row, column]),
                    stratum_width,
                    exact[row, column],
                    facecolor=color_map[name],
                    edgecolor="none",
                    alpha=1.0,
                    zorder=2,
                )
            )

    ax.set_xlim(-0.60, 15.60)
    ax.set_ylim(-5, 105)
    ax.set_xticks(range(16), TIME_ORDER, rotation=48, ha="right")
    ax.set_yticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("Time Point")
    ax.set_ylabel("Relative Abundance")
    ax.tick_params(direction="out", length=3.2, color="#333333")
    for spine in ax.spines.values():
        spine.set_color("#2F2F2F")
        spine.set_linewidth(0.85)
    fig.text(0.014, 0.982, "D", fontsize=15, ha="left", va="top")

    handles = [Patch(facecolor=color_map[name], edgecolor="none", label=name) for name in legend_order]
    handles.append(Patch(facecolor=color_map["Others"], edgecolor="none", label="Others"))
    # Three explicit columns reproduce the manuscript's 10 / 10 / 11 split.
    for left, subset in zip([0.140, 0.430, 0.660], [handles[:10], handles[10:20], handles[20:]]):
        legend = fig.legend(
            handles=subset,
            loc="upper left",
            bbox_to_anchor=(left, 0.188),
            ncol=1,
            frameon=False,
            fontsize=5.8,
            handlelength=0.86,
            handleheight=0.86,
            handletextpad=0.38,
            labelspacing=0.10,
            borderaxespad=0,
        )
        for text in legend.get_texts():
            if text.get_text() != "Others":
                text.set_fontstyle("italic")

    record = export_pdf_tiff(fig, output / "Fig1D_streamgraph_restored", dpi=1000)
    record.update(
        {
            "source_data": "species.csv, taxa.csv, meta.csv",
            "selection": "top 30 named species by mean abundance in CIP-exposed subjects",
            "aggregation": "arithmetic mean relative abundance at each canonical time",
            "interpolation": "quintic smoothstep connectors between exact time-point strata",
            "missing_data": "no canonical group-time cell missing",
        }
    )
    return record


def figure2f_tables(
    ds: audit.Dataset, tables: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    full, _ = audit.differential_proxy(ds, tables)
    full = full.copy()
    # audit.differential_proxy is HDR minus LDR; the manuscript heat map places
    # LDR-enriched taxa at the positive/red end, so the displayed direction is reversed.
    full["effect_LDR_over_HDR"] = -full["rank_biserial_HDR_over_LDR"]
    ranking = (
        full.groupby("species")["effect_LDR_over_HDR"]
        .apply(lambda values: values.abs().mean())
        .nlargest(45)
    )
    selected = ranking.index
    signed_order = (
        full.loc[full["species"].isin(selected)]
        .groupby("species")["effect_LDR_over_HDR"]
        .mean()
        .sort_values(ascending=False)
        .index
    )
    matrix = (
        full.loc[full["species"].isin(signed_order)]
        .pivot(index="canon", columns="species", values="effect_LDR_over_HDR")
        .reindex(TIME_ORDER)
        .loc[:, signed_order]
    )
    annotations = full.loc[
        full["species"].isin(signed_order),
        ["canon", "species", "effect_LDR_over_HDR", "p", "q_BH"],
    ].copy()
    annotations["symbol"] = [p_symbol(p, q) for p, q in zip(annotations["p"], annotations["q_BH"])]
    annotations["displayed_significant"] = annotations["q_BH"] < 0.05
    annotations["method"] = "q<0.05 display gate; star level from raw Mann-Whitney P; BH within canonical time"
    matrix.to_csv(tables / "Fig2F_rank_biserial_matrix.csv")
    annotations.to_csv(tables / "Fig2F_significance_annotations.csv", index=False)
    return matrix, annotations


def draw_fig2f(ds: audit.Dataset, output: Path, tables: Path) -> dict[str, object]:
    matrix, annotations = figure2f_tables(ds, tables)
    pd.DataFrame(
        {
            "normalized_coordinate": np.linspace(-1.0, 1.0, len(FIG2F_RDBU_ANCHORS)),
            "source_anchor_hex": FIG2F_RDBU_ANCHORS,
        }
    ).to_csv(tables / "Fig2F_source_palette_anchors.csv", index=False)
    # Standalone proportions retain the supplied Figure 2F crop geometry.
    fig = plt.figure(figsize=(8.137, 4.800))
    # Leave a conservative right safety margin for phase labels in PDF viewers.
    ax = fig.add_axes([0.115, 0.340, 0.855, 0.565])

    cmap = manuscript_fig2f_cmap()
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    values = matrix.to_numpy(dtype=float)
    x_edges = np.arange(values.shape[1] + 1)
    y_edges = np.arange(values.shape[0] + 1)
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        values,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors=(1, 1, 1, 0.74),
        linewidth=0.22,
        antialiased=True,
    )
    ax.invert_yaxis()
    ax.set_xlim(0, values.shape[1])
    ax.set_ylim(values.shape[0], 0)
    display_names = [name.replace(" ", "_") for name in matrix.columns]
    ax.set_xticks(np.arange(values.shape[1]) + 0.5, display_names)
    ax.tick_params(axis="x", labelrotation=45, labelsize=5.8, pad=2.5, length=0)
    for label in ax.get_xticklabels():
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    ax.set_yticks(np.arange(16) + 0.5, TIME_ORDER)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", length=0, pad=4, labelsize=8.0)
    for label in ax.get_yticklabels():
        label.set_color(PHASE_TEXT[phase(label.get_text())])

    row_pos = {name: i for i, name in enumerate(TIME_ORDER)}
    col_pos = {name: i for i, name in enumerate(matrix.columns)}
    for row in annotations.itertuples(index=False):
        if not row.symbol:
            continue
        i = row_pos[row.canon]
        j = col_pos[row.species]
        ax.text(
            j + 0.5,
            i + 0.52,
            row.symbol,
            ha="center",
            va="center",
            fontsize=5.4,
            fontweight="bold",
            color="white",
        )

    for spine in ax.spines.values():
        spine.set_visible(False)
    cax = fig.add_axes([0.040, 0.485, 0.0105, 0.300])
    colorbar = fig.colorbar(mesh, cax=cax, orientation="vertical", ticks=[-1, -0.5, 0, 0.5, 1])
    colorbar.set_label("Rank-biserial effect", fontsize=7.2, labelpad=3)
    colorbar.ax.yaxis.set_label_position("left")
    colorbar.ax.tick_params(labelsize=6.8, length=2.5, width=0.6)
    colorbar.outline.set_visible(False)
    fig.text(0.006, 0.975, "F", fontsize=16, ha="left", va="top")
    fig.text(
        0.115,
        0.930,
        "Rank-biserial proxy (LDR > 0; Mann–Whitney + within-time BH; "
        "nonlinear colour ramp γ=0.54; not LEfSe LDA)",
        ha="left",
        va="bottom",
        fontsize=6.0,
        color="#444444",
    )

    record = export_pdf_tiff(fig, output / "Fig2F_effect_heatmap_restored", dpi=600)
    record.update(
        {
            "source_data": "species.csv, taxa.csv, meta.csv",
            "grouping": "15 subjects with lowest versus 15 with highest C4-C6 mean Shannon",
            "prefilter": "species prevalence >=20% and mean abundance >=0.01% across the 30 HDR/LDR subjects",
            "selection": "45 species with largest mean absolute rank-biserial effect",
            "column_order": "descending mean signed effect (LDR over HDR)",
            "test": "two-sided Mann-Whitney within time point",
            "multiplicity": "Benjamini-Hochberg within time point across all prefilter-eligible species",
            "effect_threshold": "none; no rank-biserial analogue of the manuscript |LDA|>2 threshold was imposed",
            "manuscript_discrepancy": "supplied raster has 45 columns although its caption says top 30; 45 retained for visual reconstruction",
            "palette": "source-sampled 100-level reverse ColorBrewer RdBu LUT (#053061 to #67001F; neutral anchor #F7F7F7)",
            "color_normalization": "TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1); raw plotted values unchanged",
            "color_ramp": "symmetric nonlinear lookup ramp calibrated to source-cell color occupancy",
            "color_ramp_gamma": FIG2F_COLOR_GAMMA,
            "source_palette_reference": "antibiotic manuscript Figures/Figure2.tif, panel F color bar and 720 cell fills",
            "integrity_boundary": "proxy; not LEfSe LDA",
        }
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "selected_panels",
    )
    args = parser.parse_args()
    figures = args.output / "figures"
    tables = args.output / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    set_style()
    ds = audit.load_dataset(args.data_root)
    records = [draw_fig1d(ds, figures, tables), draw_fig2f(ds, figures, tables)]
    (args.output / "export_manifest.json").write_text(
        json.dumps(
            {
                "renderer": str(Path(__file__).resolve()),
                "random_seed": SEED,
                "records": records,
                "integrity_note": "Figure 2F is a labelled rank-biserial/Mann-Whitney proxy, not reconstructed LEfSe LDA.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (args.output / "README.md").write_text(
        """# Standalone restored panels

Run from the `reproduction/` directory with:

```bash
python3 render_selected_panels.py
```

Final deliverables are under `figures/` and contain only PDF and RGB,
LZW-compressed TIFF files. Fig1D is exported at 1000 dpi and Fig2F at
600 dpi. The CSV files under `tables/` contain every plotted value and
annotation.

## Interpretation boundary

Fig1D recomputes the 30 most abundant named species in the supplied
CIP-exposed cohort and adds the residual as `Others`. The supplied
manuscript raster appears to have been produced from a different
taxonomy snapshot, so its species set cannot be matched exactly without
the original plotting table.

Fig2F reproduces the manuscript's 16-by-45 visual structure, but is
explicitly a rank-biserial/Mann-Whitney proxy. Positive red values mean
LDR > HDR. Stars are displayed only where within-time BH q < 0.05, then
graded by raw P (* < 0.05, ** < 0.01, *** < 0.001). It is not an LEfSe
LDA reconstruction because the original LEfSe inputs and outputs were
not supplied. The analysis first retains species with prevalence >=20%
and mean abundance >=0.01% across the 15 LDR and 15 HDR subjects, then
selects 45 by mean absolute effect. BH correction is performed among all
prefilter-eligible species within each time point. No effect threshold is
used as a substitute for the manuscript's |LDA| > 2 rule. The supplied
raster itself has 45 columns even though its caption says top 30; this
version retains 45 to reconstruct the visible panel structure. The color
family is rebuilt from the manuscript color bar's exact reverse-ColorBrewer
RdBu anchors (`#053061` to `#67001F`, neutral `#F7F7F7`) and its 100-level
interpolation. Because the source LDA values occupied the color scale more
strongly than the rank-biserial proxy, a disclosed symmetric nonlinear lookup
ramp (gamma = 0.54) restores the source saturation while the plotted values and
`TwoSlopeNorm(-1, 0, 1)` remain unchanged. Palette anchors are exported as
`tables/Fig2F_source_palette_anchors.csv`.
""",
        encoding="utf-8",
    )
    print(f"Standalone panels written to {figures.resolve()}")


if __name__ == "__main__":
    main()
