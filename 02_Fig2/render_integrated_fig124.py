#!/usr/bin/env python3
"""Render the integrated manuscript-style Fig1/Fig2/Fig4 figure set.

The bundle contains the three requested main figures and their actual linked
supplementary figures (Fig1 -> sFig1, Fig2 -> sFig2, Fig4 -> sFig3).  It reuses
the audited calculations, integrates the latest Fig1D/Fig2F presentation layers,
and keeps unavailable LEfSe/FastSpar quantities explicitly labelled as proxies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile

import render_fig4_sfig4 as fig4_renderer
import render_high_fidelity as hifi
import render_selected_panels as selected
import run_all as audit


NAME_MAP = {
    "Figure1_high_fidelity": "Fig1_restored",
    "Figure2_high_fidelity": "Fig2_restored",
    "Supplementary_Figure1_high_fidelity": "sFig1_restored",
    "Supplementary_Figure2_high_fidelity": "sFig2_restored",
    "Supplementary_Figure3_high_fidelity": "sFig3_restored",
}

FINAL_STEMS = [
    "Fig1_restored",
    "Fig2_restored",
    "Fig4_restored",
    "sFig1_restored",
    "sFig2_restored",
    "sFig3_restored",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def integrated_export(
    figure,
    stem: Path,
    dpi: int = 600,
) -> dict[str, object]:
    """Map legacy renderer stems to the final PDF/TIFF-only deliverables."""
    final_name = NAME_MAP.get(stem.name, stem.name)
    final_dpi = 1000 if final_name == "Fig1_restored" else 600
    record = selected.export_pdf_tiff(
        figure,
        stem.parent / final_name,
        dpi=final_dpi,
    )
    record["canvas_policy"] = "fixed page; bbox_inches=tight disabled"
    return record


def build_comparisons(data_root: Path, output: Path) -> None:
    comparisons = output / "comparisons"
    comparisons.mkdir(parents=True, exist_ok=True)
    source_figures = data_root / "抗生素稿件 isme j" / "Figures"
    for number in [1, 2, 4]:
        hifi.comparison_sheet(
            source_figures / f"Figure{number}.tif",
            output / "figures" / f"Fig{number}_restored.tif",
            comparisons / f"Fig{number}_source_vs_restored.png",
            f"Figure {number}",
        )

    document = data_root / "抗生素稿件 isme j" / "Supplymentary material.docx"
    with tempfile.TemporaryDirectory(prefix="fig124_supp_") as temporary:
        temporary_root = Path(temporary)
        with zipfile.ZipFile(document) as archive:
            for number in [1, 2, 3]:
                source = temporary_root / f"sFig{number}_source.tif"
                source.write_bytes(archive.read(f"word/media/image{number}.tiff"))
                hifi.comparison_sheet(
                    source,
                    output / "figures" / f"sFig{number}_restored.tif",
                    comparisons / f"sFig{number}_source_vs_restored.png",
                    f"Supplementary Figure {number}",
                )


def enrich_records(records: list[dict[str, object]], figures: Path) -> None:
    for record in records:
        stem = str(record["stem"])
        for suffix in ["pdf", "tif"]:
            path = figures / f"{stem}.{suffix}"
            record[f"{suffix}_sha256"] = sha256(path)
            record[f"{suffix}_bytes"] = path.stat().st_size


def write_readme(output: Path) -> None:
    (output / "README.md").write_text(
        """# Integrated Fig1 / Fig2 / Fig4 restoration

This is the consolidated manuscript-style export.  The final `figures/`
directory intentionally contains PDF and opaque RGB/LZW TIFF files only.

## Correct main-to-supplement mapping

- Fig1 -> sFig1: longitudinal dynamics in the unexposed controls and species
  prevalence/phase comparisons.
- Fig2 -> sFig2: community structural changes, recovery trajectories, PCoA,
  and *Faecalibacterium prausnitzii*.
- Fig4 -> sFig3: phase overlap and persistence of diversity correlations.

The manuscript's sFig4 belongs to Fig5 (NAD metabolism), so it is deliberately
not included in this Fig1/2/4 bundle.

## Visual restoration decisions

- Fig1 keeps the source 7.778 x 8.600 in canvas.  Its panel D is redrawn as
  exact one-third-width time-point strata with translucent quintic alluvial
  connectors, the source-order palette and a 10/10/11 legend.  It is not a
  pasted raster panel.
- Fig2 keeps the source 8.137 x 9.447 in canvas and source-measured A-F panel
  boxes.  Panel F uses the source-hued 100-level reverse ColorBrewer RdBu LUT,
  neutral `#F7F7F7`, `TwoSlopeNorm(-1, 0, 1)`, and a disclosed symmetric
  gamma=0.54 lookup warp.  Raw plotted effect values are unchanged.
- Fig4 is the recalculated restored version with alluvia, independent diversity
  heat-map scales, network sequence and correlation small multiples.
- sFig1-sFig3 use the embedded source TIFF aspect ratios.  sFig1 restores red
  vertical P annotations and direction-coloured species headings.  sFig2 uses
  the source A/B/C(two PCoA panels)/D label grammar, centroid SEM bars and a
  white-background bottom trajectory.  sFig3 retains the dense UpSet/degree/
  exact-pattern grammar.

## Scientific boundaries

Visual similarity never substitutes for missing statistics.  Fig2F is a
rank-biserial/Mann-Whitney proxy because the original LEfSe input/output and LDA
scores were not supplied.  Fig4C is a subject-mean Spearman proxy because the
original FastSpar inputs, bootstrap runs and adjusted-P edge list were not
supplied.  The available leaf-taxon matrix is used as the disclosed richness
proxy.  Taxonomy/version differences explain remaining composition and sFig3
count differences; source values were not copied to force agreement.

## Reproduction

From `reproduction/` run:

```bash
python3 render_integrated_fig124.py
python3 test_integrated_fig124.py
```

Every plotted table is under `tables/`; the JSON manifest records dimensions,
formats, hashes and interpretation boundaries.  Side-by-side source comparisons
under `comparisons/` are visual-QA aids, not analytical outputs.
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
        default=Path(__file__).resolve().parent / "outputs" / "integrated_fig124",
    )
    args = parser.parse_args()
    figures = args.output / "figures"
    tables = args.output / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    # hifi functions create the figure artists; this exporter changes only the
    # final naming/format policy and leaves the fixed canvas untouched.
    hifi.export = integrated_export
    hifi.set_hifi_style()
    dataset = audit.load_dataset(args.data_root)

    records: list[dict[str, object]] = []
    records.append(hifi.figure1_hifi(dataset, figures, tables))
    records.append(hifi.figure2_hifi(dataset, figures, tables))
    records.append(hifi.supplementary1_hifi(dataset, figures, tables))
    records.append(hifi.supplementary2_hifi(dataset, figures, tables))
    records.append(hifi.supplementary3_hifi(dataset, figures, tables))

    # Use the stronger standalone Fig4 implementation, not the older hifi one.
    fig4_renderer.set_style()
    records.insert(2, fig4_renderer.draw_figure4(dataset, figures, tables))

    enrich_records(records, figures)
    build_comparisons(args.data_root, args.output)
    write_readme(args.output)

    manifest = {
        "renderer": str(Path(__file__).resolve()),
        "random_seed": audit.SEED,
        "main_to_supplement": {"Fig1": "sFig1", "Fig2": "sFig2", "Fig4": "sFig3"},
        "excluded_from_bundle": {
            "sFig4": "belongs to manuscript Fig5 (NAD metabolism), not Fig4"
        },
        "final_stems": FINAL_STEMS,
        "records": records,
        "integrity_note": (
            "Fig2F remains a labelled rank-biserial/Mann-Whitney proxy and "
            "Fig4C remains a labelled Spearman proxy; no unavailable values "
            "were copied from the source raster."
        ),
    }
    (args.output / "export_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Integrated Fig1/Fig2/Fig4 set written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
