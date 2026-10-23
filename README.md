# Manuscript and reviewer-response analysis code

This directory is a curated, figure-oriented copy of the scripts used for the original manuscript reconstruction and subsequent reviewer-requested analyses. Original source files remain in `reproduction/` and `reviewer_response/`; this package does not move or delete them.

## Directory structure

- `00_shared/`: common raw-data loading, shared rendering, QA, and dependency files.
- `01_Fig1/`–`04_Fig4/`: original main-figure reconstruction scripts.
- `05_sFig1/`–`08_sFig4/`: Supplementary Fig. 1–4 scripts; `08_sFig4/reference_NAD/`
  contains the corrected NAD-pathway reference mapping.
- `09_Fig5/reference/`: supplied human MetaCyc/KO, correlation, enzyme and
  Sankey reference scripts for Fig. 5.
- `10_Fig6/`: preferred animal-level mouse 16S reanalysis plus legacy QIIME 2
  and plotting reference scripts for Fig. 6.
- `11_sFig5/`: corresponding mouse 16S reference scripts for Supplementary
  Fig. 5.
- `revisions/`: reviewer-requested analyses, grouped by comment and affected panel.
- `manuscript_build/`: utilities for assembling the revised manuscript and line locations.
- `SCRIPT_INDEX.tsv`: generated inventory with language, file size, and SHA256 for every curated script.
- `SOURCE_MAP.tsv`: source-to-curated-path provenance map.

## Important conventions

1. Raw data are not duplicated here. Scripts continue to resolve inputs relative to `/home/qqt/20061009response` unless otherwise stated.
2. Test files, caches, simulated-training outputs, and superseded exploratory scripts are excluded.
3. Several original figures share `load_raw_bracken.py` and `render_high_fidelity.py`; copies are placed in each relevant figure folder for clarity, while canonical shared copies remain in `00_shared/`.
4. Fig. 3 uses the no-qPCR MDSINE2 workflow. `03_Fig3/` contains input preparation, inference, posterior summarization, longitudinal statistics, rendering, finalization, and validation scripts.
5. The revised Fig. 3 renderer containing the complete Fig. 3E/3F radial legends is also copied to `revisions/Minor4_Fig3EF/`.
6. Files below a `reference/` or `legacy_reference/` directory are supplied as
   methodological references rather than guaranteed one-command workflows.
   Several preserve historical absolute paths, environment names, and plotting
   assumptions and therefore require path review before execution.
7. The mouse experiment contains 16S rRNA amplicon data only. Accordingly,
   Fig. 6 and Supplementary Fig. 5 contain taxonomic, alpha-diversity and
   beta-diversity code, but no mouse metagenomic KO/pathway code.
8. `PANEL_SCRIPT_MAP.tsv` maps every main and supplementary panel to its best
   available script and labels it as direct, shared, reference, or unavailable.

## Rebuild the package

```bash
bash /home/qqt/20061009response/code/assemble_code_package.sh
```

The command refreshes curated copies from their canonical source locations and regenerates `SCRIPT_INDEX.tsv`.

