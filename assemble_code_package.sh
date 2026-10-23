#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/qqt/20061009response"
DEST="$ROOT/code"

copy_script() {
    local category="$1"
    local source_rel="$2"
    local target_name="${3:-$(basename "$source_rel")}" 
    local source="$ROOT/$source_rel"
    local target="$DEST/$category/$target_name"
    if [[ ! -f "$source" ]]; then
        echo "Missing required source: $source_rel" >&2
        exit 1
    fi
    install -D -m 0644 "$source" "$target"
}

# Shared original-manuscript data loading, rendering and validation.
copy_script "00_shared" "reproduction/run_all.py"
copy_script "00_shared" "reproduction/load_raw_bracken.py"
copy_script "00_shared" "reproduction/render_high_fidelity.py"
copy_script "00_shared" "reproduction/visual_qa.py"
copy_script "00_shared" "reproduction/requirements.txt"

# Original manuscript figures. Shared renderers are intentionally copied into
# each figure directory so that each directory remains self-describing.
for fig in 01_Fig1 02_Fig2 05_sFig1 06_sFig2; do
    copy_script "$fig" "reproduction/load_raw_bracken.py"
    copy_script "$fig" "reproduction/render_high_fidelity.py"
done
copy_script "01_Fig1" "reproduction/render_integrated_fig124_raw.py"
copy_script "01_Fig1" "reproduction/render_integrated_fig124.py"

copy_script "02_Fig2" "reproduction/lefse_reconstruction.py"
copy_script "02_Fig2" "reproduction/render_selected_panels.py"
copy_script "02_Fig2" "reproduction/render_integrated_fig124_raw.py"
copy_script "02_Fig2" "reproduction/render_integrated_fig124.py"

for script in prepare_inputs.py mdsine2_noqpcr_loader.py mdsine2_noqpcr_wrapper.py \
              run_inference.py summarize_posteriors.py longitudinal_stats.py \
              render_fig3.py validate_fig3.py finalize_when_ready.py; do
    copy_script "03_Fig3" "reproduction/fig3_mdsine2/$script"
done

copy_script "04_Fig4" "reproduction/fastspar_fig4.py"
copy_script "04_Fig4" "reproduction/render_fig4_sfig4.py"
copy_script "04_Fig4" "reproduction/render_integrated_fig124_raw.py"

copy_script "07_sFig3" "reproduction/audit_sfig3_variants.py"
copy_script "07_sFig3" "reproduction/render_high_fidelity.py"
# Remove two previously misclassified generated copies. Their canonical files
# remain under reproduction/ and their correct curated copies remain in Fig. 4.
rm -f "$DEST/08_sFig4/fastspar_fig4.py" "$DEST/08_sFig4/render_fig4_sfig4.py"

# The remaining manuscript figures were produced before the high-fidelity
# reconstruction pipeline.  Keep their available analysis programs as
# explicitly labelled reference implementations.  They retain the paths and
# assumptions of the supplied project and should be reviewed before rerunning.
for script in analyze_pathway_lefse_overall.R analyze_pathway_lefse.R \
              analyze_pathway_limma.R analyze_pathway_wilcox_log2fc.R \
              compare_pathway_groupings.R compare_pathway_lefse.R \
              correlate_shannon_all_pathways.R correlate_shannon_pathway.R \
              correlate_shannon_pathways_by_timepoint.R list_all_pathways.R \
              plot_combined_boxplots.R plot_pathway_volcano.R \
              plot_scatter_highlight.R plot_shannon_selected_pathways_corr.R \
              plot_shannon_selected_pathways_corr_square.R \
              plot_specific_pathway_boxplot.R; do
    copy_script "09_Fig5/reference" "ABX课题/ABX课题/Script/Fig5/$script"
done
for script in analyze_pwy7761.py analyze_pyridnucsyn.py \
              check_pyridnucsyn_details.py plot_pwy7761_subplots.py \
              plot_sankey_nad.py plot_sankey_p1_p4.py; do
    copy_script "09_Fig5/reference" "ABX课题/ABX课题/$script"
done

# Supplementary Fig. 4 uses the same baseline pathway/KO inputs as Fig. 5.
for script in correlate_shannon_pathway.R plot_specific_pathway_boxplot.R \
              plot_combined_boxplots.R; do
    copy_script "08_sFig4/reference_NAD" "ABX课题/ABX课题/Script/Fig5/$script"
done
for script in analyze_pyridnucsyn.py check_pyridnucsyn_details.py \
              plot_pwy7761_subplots.py; do
    copy_script "08_sFig4/reference_NAD" "ABX课题/ABX课题/$script"
done

# Fig. 6 and Supplementary Fig. 5 are supported by the supplied mouse 16S
# workflow.  No mouse shotgun-metagenomic pathway/KO script is included because
# that assay was not performed.  The reviewer reanalysis is the preferred
# animal-level statistical implementation; older QIIME/plotting programs are
# retained under legacy_reference for provenance only.
for category in 10_Fig6 11_sFig5; do
    copy_script "$category/preferred" "reviewer_response/R1.9_mouse_reanalysis/scripts/run_mouse_alpha_reanalysis.R"
    for script in setup_qiime2_data.py run_analysis.sh \
                  export_diversity_table.py export_diversity_table_new_meta.py \
                  export_evenness_data.py analyze_diversity.py \
                  analyze_diversity_new_meta.py plot_diversity_annotated.py \
                  plot_evenness_annotated.py; do
        copy_script "$category/legacy_reference" "ABX课题/Scripts/$script"
    done
done

# Reviewer-requested reanalyses.
for script in run_R1_1_reanalysis.R summarize_fig24_d0_vs_alternatives.R \
              summarize_fig3_across_definitions.R render_full_original_layout.py \
              run_fast_mdsine2.py validate_complete_package.py validate_R1_1.R; do
    copy_script "revisions/R1.1" "reviewer_response/R1.1/$script"
done
copy_script "revisions/R1.1" "reviewer_response/R1.1_focused_4panel/scripts/build_R1_1_focused_response.R"

for script in run_R1_2_analysis.R run_repeated_rarefaction.py validate_R1_2_outputs.R; do
    copy_script "revisions/R1.2" "reviewer_response/R1.2/$script"
done

copy_script "revisions/R1.3_DMM" "reviewer_response/R1.3_DMM/scripts/run_R1_3_DMM_analysis.R"
copy_script "revisions/R1.3_DMM" "reviewer_response/R1.3_DMM/scripts/render_R1_3_DMM_two_panel.R"
copy_script "revisions/R1.3" "reviewer_response/R1.3/scripts/run_R1_3_analysis.R"
copy_script "revisions/R1.4" "reviewer_response/R1.4/scripts/run_R1_4_flow.R"
copy_script "revisions/R1.5" "reviewer_response/R1.5/analyze_r1_5.py"

copy_script "revisions/R1.6/Fig1C" "reviewer_response/R1.6_Fig1C_permutation/scripts/reanalyse_Fig1C.R"
copy_script "revisions/R1.6/sFig2A" "reviewer_response/R1.6_sFig2A_subject_permutation/scripts/reanalyse_sFig2A.R"
copy_script "revisions/R1.6/sFig2C" "reviewer_response/R1.6_sFig2C_PERMDISP/scripts/reanalyse_sFig2C.R"
copy_script "revisions/R1.6/non_independence" "reviewer_response/R1.6_non_independence/scripts/render_R1_6_non_independence.R"
copy_script "revisions/R1.6/distance_based" "reviewer_response/R1.6_distance_based/scripts/make_R1_6_distance_based_figure.R"
copy_script "revisions/R1.6/dispersion_location" "reviewer_response/R1.6_dispersion_location/scripts/make_R1_6_dispersion_location_figure.R"
copy_script "revisions/R1.6/multiplicity" "reviewer_response/R1.6_multiplicity/scripts/make_R1_6_multiplicity_figure.R"
copy_script "revisions/R1.6/core" "reviewer_response/R1.6/scripts/run_R1_6_analysis.R"
copy_script "revisions/R1.6/core" "reviewer_response/R1.6/scripts/run_phase_clr_analysis.py"

copy_script "revisions/R1.7" "reviewer_response/R1.7/scripts/run_R1_7_analysis.R"
copy_script "revisions/R1.6_R1.7_joint" "reviewer_response/run_R1_6_R1_7_revised_analysis.R"
copy_script "revisions/R1.8" "reviewer_response/R1.8/scripts/run_R1_8_analysis.sh"
copy_script "revisions/R1.9_mouse" "reviewer_response/R1.9/scripts/audit_mouse_data_availability.R"
copy_script "revisions/R1.9_mouse" "reviewer_response/R1.9_mouse_reanalysis/scripts/run_mouse_alpha_reanalysis.R"

# Minor Comment 4 uses the revised Figure 3 renderer with complete radial keys.
copy_script "revisions/Minor4_Fig3EF" "reproduction/fig3_mdsine2/render_fig3.py"

# General reviewer-analysis pipeline and manuscript assembly utilities.
for script in core_reanalysis.py cohort_flow.py remaining_gap_audit.py; do
    copy_script "revisions/00_reviewer_core" "reviewer_response/02_analysis/$script"
done
copy_script "revisions/pathway_reanalysis" "reviewer_response/02_analysis/pathway_reanalysis/run_pathway_reanalysis.py"
copy_script "manuscript_build" "reviewer_response/manuscript_revision/scripts/build_revised_manuscript.py"
copy_script "manuscript_build" "reviewer_response/manuscript_revision/scripts/add_manuscript_line_locations.py"

# Rebuild a provenance index from the actual copied files.
python3 "$DEST/build_script_index.py"

echo "Code package assembled at $DEST"

