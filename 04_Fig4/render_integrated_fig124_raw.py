#!/usr/bin/env python3
"""Render the raw-calibrated figure set with completed FastSpar Fig4C results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import fastspar_fig4 as fastspar
import lefse_reconstruction as lefse
import load_raw_bracken as raw_loader
import render_fig4_sfig4 as fig4_renderer
import render_high_fidelity as hifi
import render_integrated_fig124 as integrated
import render_selected_panels as selected
import run_all as audit


def write_readme(output: Path, fig4_record: dict[str, object]) -> None:
    fastspar_bootstraps = int(fig4_record["panel_C_bootstraps"])
    fastspar_seed = int(fig4_record["panel_C_seed"])
    fastspar_rho = float(fig4_record["panel_C_rho_threshold_strict_gt"])
    fastspar_alpha = float(fig4_record["panel_C_alpha_strict_lt"])
    fastspar_results = str(fig4_record["panel_C_results"])
    fastspar_feature_rule = str(fig4_record["panel_C_feature_rule"])
    (output / "README.md").write_text(
        f"""# Raw-data-calibrated Fig1 / Fig2 / Fig4 restoration

This bundle is the source-faithful reconstruction of manuscript Fig1, Fig2,
Fig4 and their linked sFig1, sFig2 and sFig3. Final figures are supplied as
fixed-page PDF and opaque RGB/LZW TIFF.

## Input recovered from `raw/`

- `raw/combine_brk.mpa` is the original Kraken2/Bracken taxonomic table. Its
  all-cohort top-30 bacterial species reproduce the Fig1D legend item-for-item.
- The figure branch uses unrarefied bacterial species/genus leaves, normalized
  to proportions for Shannon and relative abundance. This reproduces Fig2D/E
  marker coordinates and gives the printed Fig2C statistics after rounding:
  LDR r=0.82, P=0.0002; HDR r=0.81, P=0.0003.
- HDR/LDR labels are the bottom/top 15 subjects in the supplied `species.csv`
  C4-C6 ranking. The grouping variable and plotted Bracken alpha values remain
  separate in the audit tables.

## Manuscript Methods discrepancy

The Methods state rarefaction to 2,693,000 reads per sample, but rarefied
values systematically miss the source Fig2D/E coordinates. The final figures
therefore use the demonstrably source-faithful unrarefied branch. A deterministic
Methods-faithful rarefied sensitivity table is retained under `tables/`; it is
not substituted into the source restoration.

## Figure-specific reconstruction

- Fig1D selects the top 30 species across all 1,040 samples, then plots the
  longitudinal mean in the 960 CIP samples, matching the original legend.
- Fig2F is no longer a rank-biserial stand-in. It runs the installed LEfSe
  1.1.2/MASS LDA bootstrap at each time point on all 268 eligible species.
  Every source cell retains its LDA magnitude, the sign follows the raw
  LDR-minus-HDR mean, and stars use exact Mann-Whitney P values BH-adjusted
  over the 45 displayed species. The raster has 45 columns despite its caption
  saying 30; their fixed source order and a pixel-digitized audit are exported.
- sFig2A displays all 225 fixed LDR×HDR genus-presence Jaccard pairs per
  time point (source-calibrated detection threshold >0.0001%); each pair is
  tested against its own A1 value by paired Wilcoxon. Unthresholded genus and
  quantitative-species alternatives are exported as sensitivities. sFig2C/D
  use the abundance-weighted vegan/phyloseq Jaccard default
  (`J = 2*Bray/(1+Bray)`), rather than binary presence Jaccard; axis labels use
  ape's uncorrected trace denominator, matching the phyloseq PCoA path.
- Fig4A/B/D and sFig1-sFig3 now use the recovered Bracken species/genus table.
  Fig4C uses the completed FastSpar reconstruction at `{fastspar_results}`:
  {fastspar_bootstraps} bootstraps, seed {fastspar_seed}, and the caption rule
  `|rho| > {fastspar_rho:g}` plus exact bootstrap raw `P < {fastspar_alpha:g}`.
  The separate Methods-rule edge table applies within-network BH
  `Q < {fastspar_alpha:g}` across all candidate pairs.
  Both are retained because the manuscript caption and Methods specify
  different significance rules. The complete caption/raw-P list is the primary
  analytical result and is exported unchanged. To restore the visibly lower
  network density in the manuscript raster, panel C draws a transparent,
  deterministic subset only: it first retains a maximum-spanning forest by
  absolute rho (preserving every eligible node), then fills the strongest
  remaining eligible edges to the approximate source-raster edge count. If the
  target exceeds the available eligible edges, all eligible edges are drawn.
  No edge is synthesized, and the edge-count trajectory reports this explicitly
  labelled display subset rather than silently redefining statistical inclusion.
  The candidate rule is: {fastspar_feature_rule}. This is a transparent rerun
  from the supplied Bracken data, not the unavailable original author run; no
  edge is copied from the raster.

## Reproduction

From `reproduction/`:

```bash
python3 fastspar_fig4.py
python3 render_integrated_fig124_raw.py --fastspar-results outputs/fastspar_fig4/results
python3 test_integrated_fig124_raw.py
```

The default FastSpar command writes `outputs/fastspar_fig4/results/` only after
running all 12 HDR/LDR-stage panels with 500 bootstraps each. The renderer
rejects benchmark, partial, hash-mismatched, or wrong-source result roots.

The `tables/` directory contains plotted matrices, inference tables, source
hashes, rarefaction sensitivity values and HDR/LDR assignments. Fig4C exports
`Fig4C_fastspar_summary.csv`, `Fig4C_fastspar_edges_caption_raw_p.csv`, and
`Fig4C_fastspar_edges_methods_BH_q.csv`. The exact drawn subset and its
maximum-spanning-forest/strongest-fill role are in
`Fig4C_fastspar_edges_display_source_density.csv`. `comparisons/` contains
Python-generated visual-QA sheets and is not an analytical output.

## Visualization workflow reference

Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). Scientific Agent
Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065.
https://doi.org/10.48550/arXiv.2609.00065
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
        default=Path(__file__).resolve().parent / "outputs" / "integrated_fig124_raw",
    )
    parser.add_argument(
        "--fastspar-results",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "fastspar_fig4" / "results",
        help=(
            "completed full-run result root written by fastspar_fig4.py; "
            "benchmark and partial result roots are rejected"
        ),
    )
    args = parser.parse_args()
    completed_fastspar = fastspar.load_completed_results(args.fastspar_results)
    figures = args.output / "figures"
    tables = args.output / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    hifi.export = integrated.integrated_export
    selected.figure2f_tables = lefse.figure2f_tables
    hifi.set_hifi_style()
    dataset = raw_loader.load_raw_bracken_dataset(args.data_root)
    fig4_renderer.validate_fastspar_provenance(dataset, completed_fastspar)
    raw_loader.write_raw_audit(dataset, tables)

    records: list[dict[str, object]] = []
    records.append(hifi.figure1_hifi(dataset, figures, tables))
    records.append(hifi.figure2_hifi(dataset, figures, tables))
    fig4_renderer.set_style()
    fig4_record = fig4_renderer.draw_figure4(
        dataset,
        figures,
        tables,
        fastspar_results=completed_fastspar,
    )
    records.append(fig4_record)
    records.append(hifi.supplementary1_hifi(dataset, figures, tables))
    records.append(hifi.supplementary2_hifi(dataset, figures, tables))
    records.append(hifi.supplementary3_hifi(dataset, figures, tables))

    integrated.enrich_records(records, figures)
    integrated.build_comparisons(args.data_root, args.output)
    write_readme(args.output, fig4_record)
    manifest = {
        "renderer": str(Path(__file__).resolve()),
        "random_seed": audit.SEED,
        "taxonomic_input": "raw/combine_brk.mpa",
        "taxonomic_branch": "source-faithful unrarefied bacterial leaves",
        "methods_sensitivity": "2,693,000-read rarefaction exported separately",
        "grouping_input": "species.csv C4-C6 endpoint Shannon ranking",
        "sfig2a_jaccard": (
            "225 fixed LDRxHDR genus binary-Jaccard pairs; displayed presence "
            "threshold >0.0001%; paired Wilcoxon versus pair-matched A1"
        ),
        "sfig2cd_jaccard": (
            "species quantitative vegan/phyloseq Jaccard, J=2*Bray/(1+Bray); "
            "ape::pcoa trace-denominator axis percentages"
        ),
        "fig4c_network_results": fig4_record["panel_C_results"],
        "fig4c_bootstraps": fig4_record["panel_C_bootstraps"],
        "fig4c_seed": fig4_record["panel_C_seed"],
        "fig4c_rho_threshold_strict_gt": fig4_record[
            "panel_C_rho_threshold_strict_gt"
        ],
        "fig4c_alpha_strict_lt": fig4_record["panel_C_alpha_strict_lt"],
        "fig4c_rule_discrepancy": (
            "full caption/raw-P and Methods/BH-Q analytical edge lists are exported; "
            "the network drawing and edge-count trajectory use a transparent "
            "deterministic source-density subset of caption-qualified edges"
        ),
        "fig4c_display_rule": fig4_record["panel_C_display_rule"],
        "fig4c_display_is_full_analytic_edge_set": fig4_record[
            "panel_C_display_is_full_analytic_edge_set"
        ],
        "fig4c_display_edge_targets": fig4_record["panel_C_display_edge_targets"],
        "main_to_supplement": {"Fig1": "sFig1", "Fig2": "sFig2", "Fig4": "sFig3"},
        "final_stems": integrated.FINAL_STEMS,
        "records": records,
        "integrity_note": (
            "Fig2F is reconstructed with LEfSe LDA from the supplied raw table. "
            "sFig2A explicitly marks its raster-calibrated genus detection "
            "threshold and exports unthresholded/quantitative sensitivities; "
            "sFig2C/D follow the documented vegan quantitative-Jaccard default. "
            "Fig4C uses a completed FastSpar bootstrap reconstruction from the "
            "supplied Bracken counts. Full caption/raw-P and Methods/BH-Q edges "
            "are exported; only the rendered density uses a clearly documented "
            "maximum-spanning-forest plus strongest-edge visual subset. The "
            "reconstruction is not represented as the unavailable original author run."
        ),
    }
    (args.output / "export_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Raw-data-calibrated figure set written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
