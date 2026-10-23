#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  root <- normalizePath(getwd(), mustWork = TRUE)
}
out <- file.path(root, "reviewer_response", "R1.1")
checks <- list()

add_check <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name,
    status = ifelse(isTRUE(pass), "PASS", "FAIL"),
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

membership <- read_csv(file.path(out, "R1_1_group_manifest_D0_D6.csv"), show_col_types = FALSE)
phase <- read_csv(file.path(out, "R1_1_phase_effects.csv"), show_col_types = FALSE)
stage <- read_csv(file.path(out, "R1_1_stagewise_effects.csv"), show_col_types = FALSE)
trajectory <- read_csv(file.path(out, "R1_1_trajectory_LRT.csv"), show_col_types = FALSE)
correlation <- read_csv(file.path(out, "R1_1_Fig2E_correlations.csv"), show_col_types = FALSE)
tax_phase <- read_csv(file.path(out, "R1_1_taxonomic_phase_effects.csv"), show_col_types = FALSE)
tax_stage <- read_csv(file.path(out, "R1_1_taxonomic_stagewise_effects.csv"), show_col_types = FALSE)
visual_qa <- read_csv(file.path(out, "R1_1_visual_QA.csv"), show_col_types = FALSE)
source_subject <- read_csv(
  file.path(root, "reviewer_response", "03_tables", "core_subject_recovery_metrics.csv"),
  show_col_types = FALSE
) |>
  filter(dataset_branch == "rarefied_primary")

add_check("definitions", n_distinct(membership$definition_id) == 7L,
          n_distinct(membership$definition_id), 7)
add_check("subjects per definition",
          all((membership |> count(definition_id))$n == 60L),
          paste(sort(unique((membership |> count(definition_id))$n)), collapse = ","), 60)
group_counts <- membership |>
  filter(recovery_group %in% c("HDR", "LDR")) |>
  count(definition_id, recovery_group)
add_check("15 HDR per definition", all(group_counts$n[group_counts$recovery_group == "HDR"] == 15L),
          paste(group_counts$n[group_counts$recovery_group == "HDR"], collapse = ","), 15)
add_check("15 LDR per definition", all(group_counts$n[group_counts$recovery_group == "LDR"] == 15L),
          paste(group_counts$n[group_counts$recovery_group == "LDR"], collapse = ","), 15)

d0_manifest <- membership |>
  filter(definition_id == "D0") |>
  select(subject, observed = recovery_group)
d0_source <- source_subject |>
  transmute(subject, expected = original_endpoint_group)
d0_compare <- inner_join(d0_manifest, d0_source, by = "subject")
add_check("D0 manuscript membership preserved",
          all(d0_compare$observed == d0_compare$expected),
          sum(d0_compare$observed == d0_compare$expected), 60)

add_check("phase result rows", nrow(phase) == 105L, nrow(phase), 105)
add_check("stage result rows", nrow(stage) == 560L, nrow(stage), 560)
add_check("trajectory model rows", nrow(trajectory) == 35L, nrow(trajectory), 35)
add_check("Fig2E correlation rows", nrow(correlation) == 14L, nrow(correlation), 14)
add_check("taxonomic phase rows", nrow(tax_phase) == 63L, nrow(tax_phase), 63)
add_check("taxonomic stage rows", nrow(tax_stage) == 336L, nrow(tax_stage), 336)

numeric_complete <- function(x, cols) all(vapply(x[cols], function(z) all(is.finite(z)), logical(1)))
add_check("phase statistics finite",
          numeric_complete(phase, c("standardized_effect", "ci_low", "ci_high", "p_value", "q_within_definition")),
          "checked", "all finite")
add_check("stage statistics finite",
          numeric_complete(stage, c("standardized_effect", "ci_low", "ci_high", "p_value", "q_within_definition_metric")),
          "checked", "all finite")
add_check("trajectory models complete", all(trajectory$status == "complete"),
          paste(table(trajectory$status), collapse = ","), "35 complete")
add_check("visual preview QA", all(visual_qa$preview_status == "PASS"),
          paste(visual_qa$preview_status, collapse = ","), "all PASS")

expected_files <- c(
  "R1_1_Fig2_4_cross_definition.pdf", "R1_1_Fig2_4_cross_definition.tif",
  "R1_1_Fig2_4_stagewise.pdf", "R1_1_Fig2_4_stagewise.tif",
  "R1_1_taxonomic_sensitivity.pdf", "R1_1_taxonomic_sensitivity.tif",
  "R1.1_response_bilingual.md", "R1.1_response_bilingual.docx",
  "R1_1_analysis_manifest.json"
)
expected_paths <- file.path(out, expected_files)
add_check("deliverables exist", all(file.exists(expected_paths)),
          sum(file.exists(expected_paths)), length(expected_paths))
add_check("deliverables non-empty", all(file.info(expected_paths)$size > 0),
          min(file.info(expected_paths)$size), ">0 bytes")

result <- bind_rows(checks)
write_tsv(result, file.path(out, "R1_1_validation.tsv"))

failed <- result |> filter(status != "PASS")
summary_text <- c(
  sprintf("R1.1 validation: %d PASS, %d FAIL", sum(result$status == "PASS"), nrow(failed)),
  if (nrow(failed) == 0) "All automated checks passed." else paste("Failed:", paste(failed$check, collapse = "; "))
)
writeLines(summary_text, file.path(out, "R1_1_validation_summary.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")

if (nrow(failed) > 0) quit(status = 1)
