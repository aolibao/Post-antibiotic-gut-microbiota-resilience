#!/usr/bin/env python3
"""Generate small, local-evidence tables for unresolved reviewer comments.

This is an audit helper, not a replacement analysis.  It records the
multiplicity families used by the reconstructed Fig. 2F, verifies the three
taxonomic-name examples against the supplied Bracken MPA paths, and lists the
biological unit/n that can (or cannot) be established for selected boxplots.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "reviewer_response" / "03_tables"
FIG2_TABLES = (
    ROOT / "reproduction" / "outputs" / "integrated_fig124_raw" / "tables"
)


def bh_adjust(values: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg q values with finite-value preservation."""
    p = values.to_numpy(dtype=float)
    result = np.full(p.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(p))
    if finite.size == 0:
        return result
    order = finite[np.argsort(p[finite], kind="mergesort")]
    ranked = p[order] * finite.size / np.arange(1, finite.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result[order] = np.minimum(ranked, 1.0)
    return result


def write_lefse_multiplicity_audit() -> None:
    eligible = pd.read_csv(FIG2_TABLES / "Fig2F_LEfSe_all_eligible_species.csv")
    displayed = pd.read_csv(FIG2_TABLES / "Fig2F_LEfSe_annotations.csv")
    eligible_global_q = bh_adjust(eligible["p_exact_mannwhitney"])
    displayed_global_q = bh_adjust(displayed["p_exact_mannwhitney"])

    rows = [
        {
            "family": "eligible_268_species_within_each_time",
            "n_tests": len(eligible),
            "correction_scope": "16 separate BH families of 268 species",
            "n_q_lt_0_05": int(eligible["q_BH_candidate_universe"].lt(0.05).sum()),
            "n_q_lt_0_05_and_abs_LDA_ge_2": int(
                (
                    eligible["q_BH_candidate_universe"].lt(0.05)
                    & eligible["passes_caption_lda"]
                ).sum()
            ),
            "used_for_current_fig2f_stars": False,
            "audit_note": "Controls species within a time point, not the 16-time family.",
        },
        {
            "family": "eligible_268_species_global_all_times",
            "n_tests": len(eligible),
            "correction_scope": "one BH family over 268 species x 16 times",
            "n_q_lt_0_05": int(np.count_nonzero(eligible_global_q < 0.05)),
            "n_q_lt_0_05_and_abs_LDA_ge_2": int(
                np.count_nonzero(
                    (eligible_global_q < 0.05)
                    & eligible["passes_caption_lda"].to_numpy(dtype=bool)
                )
            ),
            "used_for_current_fig2f_stars": False,
            "audit_note": (
                "A valid across-time family for this exploratory grid, but it does not "
                "remove outcome-group selection or repeated-time interpretive limits."
            ),
        },
        {
            "family": "displayed_45_species_within_each_time",
            "n_tests": len(displayed),
            "correction_scope": "16 separate BH families of 45 displayed species",
            "n_q_lt_0_05": int(displayed["q_BH_displayed_45"].lt(0.05).sum()),
            "n_q_lt_0_05_and_abs_LDA_ge_2": int(
                (
                    displayed["q_BH_displayed_45"].lt(0.05)
                    & displayed["passes_caption_lda"]
                ).sum()
            ),
            "used_for_current_fig2f_stars": True,
            "audit_note": (
                "This is the reconstructed visual star rule; the 45 displayed taxa were "
                "selected from the effect display, so it is not a confirmatory family."
            ),
        },
        {
            "family": "displayed_45_species_global_all_times",
            "n_tests": len(displayed),
            "correction_scope": "one BH family over 45 displayed species x 16 times",
            "n_q_lt_0_05": int(np.count_nonzero(displayed_global_q < 0.05)),
            "n_q_lt_0_05_and_abs_LDA_ge_2": int(
                np.count_nonzero(
                    (displayed_global_q < 0.05)
                    & displayed["passes_caption_lda"].to_numpy(dtype=bool)
                )
            ),
            "used_for_current_fig2f_stars": False,
            "audit_note": (
                "Across-time correction only within a post-selected display set; still "
                "unsuitable as primary inference."
            ),
        },
    ]
    pd.DataFrame(rows).to_csv(TABLES / "lefse_multiplicity_audit.csv", index=False)


def write_taxonomy_harmonization_audit() -> None:
    paths = pd.read_csv(
        ROOT / "raw" / "combine_brk.mpa", sep="\t", usecols=[0]
    ).iloc[:, 0].astype(str)
    rows = [
        {
            "recommended_display_name": "Phocaeicola vulgatus",
            "legacy_name_to_harmonize": "Bacteroides vulgatus",
            "expected_raw_terminal": "s__Phocaeicola_vulgatus",
        },
        {
            "recommended_display_name": "Segatella copri",
            "legacy_name_to_harmonize": "Prevotella copri",
            "expected_raw_terminal": "s__Segatella_copri",
        },
        {
            "recommended_display_name": "Mediterraneibacter gnavus",
            "legacy_name_to_harmonize": "[Ruminococcus] gnavus / Ruminococcus gnavus",
            "expected_raw_terminal": "s__Mediterraneibacter_gnavus",
        },
    ]
    for row in rows:
        matches = paths[paths.str.endswith("|" + row["expected_raw_terminal"])]
        row["raw_bracken_path_count"] = int(len(matches))
        row["raw_bracken_path"] = matches.iloc[0] if len(matches) == 1 else ""
        row["audit_boundary"] = (
            "Name is verified against the supplied MPA path; taxonomy snapshot/build "
            "date and accession list remain author input."
        )
    pd.DataFrame(rows).to_csv(TABLES / "taxonomy_harmonization_audit.csv", index=False)


def write_boxplot_reporting_audit() -> None:
    meta = pd.read_csv(ROOT / "raw" / "meta.csv")
    exposed_n = int(meta.loc[meta["abx"].eq("Y"), "subject"].nunique())
    control_n = int(meta.loc[meta["abx"].eq("N"), "subject"].nunique())
    rows = [
        {
            "panel": "Figure 1B",
            "locally_auditable_unit": "participant at each canonical stage",
            "locally_auditable_n": exposed_n,
            "status": "auditable",
            "required_legend_action": (
                "State n=60 participants at each stage and define median, IQR, whiskers "
                "and any displayed points."
            ),
        },
        {
            "panel": "Supplementary Figure 1A/B",
            "locally_auditable_unit": "unexposed participant at each canonical stage",
            "locally_auditable_n": control_n,
            "status": "auditable_descriptive_only",
            "required_legend_action": (
                "State n=5 and descriptive-reference status; A1 self-comparisons should "
                "be labelled reference or left unannotated, not NA/NK."
            ),
        },
        {
            "panel": "Supplementary Figure 1D (actual plotted taxon facets)",
            "locally_auditable_unit": "participant-level phase mean, paired A versus C",
            "locally_auditable_n": exposed_n,
            "status": "legend_mismatch",
            "required_legend_action": (
                "Replace the supplied duplicate Jaccard legend with the actual paired "
                "taxon-abundance comparison and define box/point elements."
            ),
        },
        {
            "panel": "Supplementary Figure 4B/C",
            "locally_auditable_unit": "should be one baseline aggregate per participant",
            "locally_auditable_n": "15/group only if original late-diversity strata are retained",
            "status": "original_plotting_unit_not_established",
            "required_legend_action": (
                "Do not count A1-A4 samples as independent replicates; report participant "
                "aggregation and n explicitly."
            ),
        },
        {
            "panel": "Supplementary Figure 5A/C/D",
            "locally_auditable_unit": "not established; animal/cage tables absent",
            "locally_auditable_n": "not auditable from local source files",
            "status": "author_input_required",
            "required_legend_action": (
                "Supply animal IDs, group allocation, cage IDs, exclusions and raw values; "
                "then state animal n and whether points are animals or technical values."
            ),
        },
    ]
    pd.DataFrame(rows).to_csv(TABLES / "boxplot_reporting_audit.csv", index=False)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    write_lefse_multiplicity_audit()
    write_taxonomy_harmonization_audit()
    write_boxplot_reporting_audit()
    print("Remaining-gap audit tables written to", TABLES)


if __name__ == "__main__":
    main()
