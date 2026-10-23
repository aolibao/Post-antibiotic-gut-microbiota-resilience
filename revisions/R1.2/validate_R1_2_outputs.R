#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(digest)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve this script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
analysis_dir <- dirname(script_path)
root_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()
add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name,
    pass = isTRUE(pass),
    detail = as.character(detail)
  )
}

required <- c(
  file.path(analysis_dir, "R1_2_analysis_outline.md"),
  file.path(analysis_dir, "README.md"),
  file.path(analysis_dir, "run_R1_2_analysis.R"),
  file.path(analysis_dir, "run_repeated_rarefaction.py"),
  file.path(analysis_dir, "run_R1_2_pipeline.sh"),
  file.path(analysis_dir, "R1_2_results_summary_bilingual.md"),
  file.path(analysis_dir, "R1_2_figure_legend.md"),
  file.path(analysis_dir, "cache", "R1_2_repeated_rarefaction_subject_metrics.csv"),
  file.path(table_dir, "R1_2_primary_full_cohort_model.csv"),
  file.path(table_dir, "R1_2_D1_D6_continuous_models.csv"),
  file.path(table_dir, "R1_2_cutoff_sensitivity.csv"),
  file.path(table_dir, "R1_2_cutoff_group_assignments.csv"),
  file.path(table_dir, "R1_2_depth_sensitivity_models.csv"),
  file.path(table_dir, "R1_2_repeated_rarefaction_models.csv"),
  file.path(table_dir, "R1_2_repeated_rarefaction_summary.csv"),
  file.path(table_dir, "R1_2_repeated_rarefaction_validation.csv"),
  file.path(table_dir, "R1_2_baseline_Shannon_reliability.csv"),
  file.path(table_dir, "R1_2_primary_model_diagnostics.csv"),
  file.path(table_dir, "R1_2_subject_analysis_dataset.csv"),
  file.path(figure_dir, "R1_2_response_figure.pdf"),
  file.path(figure_dir, "R1_2_response_figure.tif"),
  file.path(figure_dir, "R1_2_response_figure_preview.png"),
  file.path(figure_dir, "R1_2_primary_model_diagnostics.pdf"),
  file.path(figure_dir, "R1_2_primary_model_diagnostics.tif"),
  file.path(log_dir, "R1_2_analysis_manifest.json"),
  file.path(log_dir, "R1_2_repeated_rarefaction_manifest.json")
)
for (path in required) {
  add_check(
    paste0("required_file:", basename(path)),
    file.exists(path) && file.info(path)$size > 0,
    if (file.exists(path)) paste0(file.info(path)$size, " bytes") else "missing"
  )
}

primary <- read_csv(file.path(table_dir, "R1_2_primary_full_cohort_model.csv"), show_col_types = FALSE)
add_check("primary_unique_D5", nrow(primary) == 1L && primary$metric_id == "D5" && isTRUE(primary$primary), "one prespecified D5 row")
add_check("primary_n60", primary$n_subjects == 60L, paste0("n=", primary$n_subjects))
add_check("primary_finite_effect", all(is.finite(unlist(primary[, c(
  "standardized_beta", "beta_HC3_SE", "beta_CI95_low", "beta_CI95_high",
  "beta_p_value", "R2", "participant_bootstrap_CI95_low",
  "participant_bootstrap_CI95_high", "permutation_p_value"
)]))), "all primary effect fields finite")
add_check("primary_interval_order", primary$beta_CI95_low < primary$standardized_beta &&
            primary$standardized_beta < primary$beta_CI95_high, "HC3 estimate lies inside CI")

models <- read_csv(file.path(table_dir, "R1_2_D1_D6_continuous_models.csv"), show_col_types = FALSE)
add_check("D1_D6_complete", identical(sort(models$metric_id), paste0("D", 1:6)), paste(sort(models$metric_id), collapse = ","))
add_check("only_D5_primary", sum(models$primary) == 1L && models$metric_id[models$primary] == "D5", "D5 only")
add_check("D1_D6_all_n60", all(models$n_subjects == 60L), paste(unique(models$n_subjects), collapse = ","))
add_check("secondary_Holm_only", all(is.na(models$p_Holm_secondary[models$primary])) &&
            all(is.finite(models$p_Holm_secondary[!models$primary])), "five adjusted secondary P values")

subject <- read_csv(file.path(table_dir, "R1_2_subject_analysis_dataset.csv"), show_col_types = FALSE)
add_check("subject_grid_n60", nrow(subject) == 60L && n_distinct(subject$subject) == 60L, "60 unique participants")
add_check("subject_primary_complete", all(subject$analysis_inclusion == "Included in full-cohort primary model"), "all included")
add_check("subject_no_missing_primary", all(is.finite(subject$baseline_shannon)) &&
            all(is.finite(subject$recovery_auc_bray_actual_day)), "complete predictor and D5")

cutoff <- read_csv(file.path(table_dir, "R1_2_cutoff_sensitivity.csv"), show_col_types = FALSE)
d5_cutoff <- cutoff |> filter(selection_definition == "D5_recovery")
add_check("cutoff_unique_12_15_20", identical(sort(d5_cutoff$n_per_tail), c(12, 15, 20)), paste(sort(d5_cutoff$n_per_tail), collapse = ","))
add_check("cutoff_no_boundary_ties", !any(d5_cutoff$high_boundary_tie | d5_cutoff$low_boundary_tie), "no ties at any D5 boundary")
add_check(
  "cutoff_full_cohort_reranking_bootstrap",
  all(grepl("Full-cohort pairs bootstrap.*membership recomputed", d5_cutoff$bootstrap_design)),
  paste(unique(d5_cutoff$bootstrap_design), collapse = "; ")
)
assignment <- read_csv(file.path(table_dir, "R1_2_cutoff_group_assignments.csv"), show_col_types = FALSE)
assignment_counts <- assignment |>
  filter(selection_definition == "D5_recovery") |>
  count(scheme, group) |>
  pivot_wider(names_from = group, values_from = n)
add_check("cutoff_assignment_counts",
          all(assignment_counts$High == c(15, 12, 20)[match(assignment_counts$scheme, c("quartile", "quintile", "tertile"))]) &&
            all(assignment_counts$Low == assignment_counts$High),
          paste(apply(assignment_counts, 1, paste, collapse = ":"), collapse = "; "))

rarefaction <- read_csv(file.path(table_dir, "R1_2_repeated_rarefaction_models.csv"), show_col_types = FALSE)
add_check("rarefaction_200_complete", nrow(rarefaction) == 200L &&
            isTRUE(all.equal(as.numeric(sort(rarefaction$replicate)), as.numeric(1:200))), "200 model rows")
add_check("rarefaction_all_n60", all(rarefaction$n_subjects == 60L), "all replicates n=60")
add_check(
  "rarefaction_seed_sequence",
  identical(as.integer(rarefaction$seed), as.integer(20261009L + 0:199)),
  sprintf("%d through %d", min(rarefaction$seed), max(rarefaction$seed))
)
rarefaction_cache <- read_csv(
  file.path(analysis_dir, "cache", "R1_2_repeated_rarefaction_subject_metrics.csv"),
  show_col_types = FALSE
)
rarefaction_cache_grid <- rarefaction_cache |>
  group_by(replicate) |>
  summarise(n_rows = n(), n_subjects = n_distinct(subject), .groups = "drop")
add_check(
  "rarefaction_cache_200x60_unique",
  nrow(rarefaction_cache_grid) == 200L &&
    all(rarefaction_cache_grid$n_rows == 60L) &&
    all(rarefaction_cache_grid$n_subjects == 60L),
  sprintf("%d replicates; row range %d-%d; subject range %d-%d",
          nrow(rarefaction_cache_grid), min(rarefaction_cache_grid$n_rows),
          max(rarefaction_cache_grid$n_rows), min(rarefaction_cache_grid$n_subjects),
          max(rarefaction_cache_grid$n_subjects))
)
rarefaction_validation <- read_csv(file.path(table_dir, "R1_2_repeated_rarefaction_validation.csv"), show_col_types = FALSE)
add_check("rarefaction_seed1_reproduces_primary",
          max(rarefaction_validation$maximum_absolute_difference) < 1e-10,
          sprintf("max delta %.3g", max(rarefaction_validation$maximum_absolute_difference)))

depth <- read_csv(file.path(table_dir, "R1_2_depth_sensitivity_models.csv"), show_col_types = FALSE)
depth_models <- depth |> filter(term == "z_baseline_shannon") |> pull(model)
add_check("depth_models_complete",
          identical(sort(depth_models), sort(c("rarefied_unadjusted", "unrarefied_unadjusted", "unrarefied_depth_adjusted"))),
          paste(depth_models, collapse = ","))
add_check("depth_models_all_n60", all(depth$n_subjects == 60L), paste(unique(depth$n_subjects), collapse = ","))

icc <- read_csv(file.path(table_dir, "R1_2_baseline_Shannon_reliability.csv"), show_col_types = FALSE)
add_check("ICC_complete", identical(sort(icc$reliability), c("ICC_A_1", "ICC_A_k")) &&
            all(icc$CI95_low <= icc$estimate & icc$estimate <= icc$CI95_high),
          "ICC(A,1) and ICC(A,4), ordered intervals")

manifest_path <- file.path(log_dir, "R1_2_analysis_manifest.json")
manifest <- fromJSON(manifest_path, simplifyVector = FALSE)
manifest_output_status <- vapply(manifest$outputs, function(entry) {
  output_path <- file.path(root_dir, entry$path)
  file.exists(output_path) &&
    identical(unname(digest(output_path, algo = "sha256", file = TRUE)), entry$sha256)
}, logical(1))
add_check(
  "manifest_output_hashes_current",
  length(manifest_output_status) > 0L && all(manifest_output_status),
  sprintf("%d/%d output hashes match", sum(manifest_output_status), length(manifest_output_status))
)
manifest_output_paths <- vapply(
  manifest$outputs, function(entry) file.path(root_dir, entry$path), character(1)
)
add_check(
  "manifest_newer_than_recorded_outputs",
  file.info(manifest_path)$mtime >= max(file.info(manifest_output_paths)$mtime),
  paste0("manifest mtime ", file.info(manifest_path)$mtime)
)

identify_bin <- Sys.which("identify")
if (nzchar(identify_bin)) {
  for (stem in c("R1_2_response_figure", "R1_2_primary_model_diagnostics")) {
    path <- file.path(figure_dir, paste0(stem, ".tif"))
    metadata <- system2(identify_bin,
                        c("-format", shQuote("%m|%w|%h|%x|%y|%U|%[compression]|%[colorspace]|%[channels]"),
                          shQuote(path)),
                        stdout = TRUE, stderr = TRUE)
    fields <- strsplit(metadata[1], "\\|")[[1]]
    add_check(paste0(stem, "_tiff_format"), length(fields) >= 9L && fields[1] == "TIFF", metadata[1])
    add_check(paste0(stem, "_tiff_600dpi"), length(fields) >= 6L &&
                abs(as.numeric(fields[4]) - 600) < 1 && abs(as.numeric(fields[5]) - 600) < 1,
              metadata[1])
    add_check(paste0(stem, "_tiff_LZW"), length(fields) >= 7L && fields[7] == "LZW", metadata[1])
    add_check(paste0(stem, "_tiff_sRGB"), length(fields) >= 9L &&
                fields[8] == "sRGB" && grepl("srgb", fields[9], ignore.case = TRUE), metadata[1])
  }
} else {
  add_check("identify_available", FALSE, "ImageMagick identify not found")
}

pdfinfo_bin <- Sys.which("pdfinfo")
if (nzchar(pdfinfo_bin)) {
  for (stem in c("R1_2_response_figure", "R1_2_primary_model_diagnostics")) {
    path <- file.path(figure_dir, paste0(stem, ".pdf"))
    info <- system2(pdfinfo_bin, path, stdout = TRUE, stderr = TRUE)
    page_line <- grep("^Pages:", info, value = TRUE)
    add_check(paste0(stem, "_pdf_one_page"), length(page_line) == 1L &&
                as.integer(sub("^Pages:\\s+", "", page_line)) == 1L,
              paste(info[grepl("^(Pages|Page size|File size):", info)], collapse = "; "))
  }
} else {
  add_check("pdfinfo_available", FALSE, "pdfinfo not found")
}

pdffonts_bin <- Sys.which("pdffonts")
if (nzchar(pdffonts_bin)) {
  for (stem in c("R1_2_response_figure", "R1_2_primary_model_diagnostics")) {
    path <- file.path(figure_dir, paste0(stem, ".pdf"))
    font_info <- system2(pdffonts_bin, path, stdout = TRUE, stderr = TRUE)
    font_rows <- font_info[-seq_len(min(2L, length(font_info)))]
    embedded <- length(font_rows) > 0L && all(grepl(
      "\\s+yes\\s+(yes|no)\\s+(yes|no)\\s+[0-9]+\\s+[0-9]+\\s*$",
      font_rows
    ))
    add_check(paste0(stem, "_pdf_fonts_embedded"), embedded,
              paste(font_rows, collapse = "; "))
  }
} else {
  add_check("pdffonts_available", FALSE, "pdffonts not found")
}

validation <- bind_rows(checks)
write_tsv(validation, file.path(log_dir, "R1_2_validation.tsv"))
summary_lines <- c(
  "# R1.2 validation summary",
  "",
  sprintf("- Checks passed: %d/%d", sum(validation$pass), nrow(validation)),
  sprintf("- Validation time: %s", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  "",
  "The gate verifies file presence, inferential hierarchy, cohort/cutoff grids,",
  "rarefaction completeness and seed-1 reproduction, model families, and",
  "machine-readable PDF/TIFF properties. It does not certify scientific or",
  "journal acceptance."
)
writeLines(summary_lines, file.path(log_dir, "R1_2_validation_summary.md"))

if (any(!validation$pass)) {
  print(validation |> filter(!pass), n = Inf)
  stop(sprintf("R1.2 validation failed: %d/%d checks failed", sum(!validation$pass), nrow(validation)))
}
cat(sprintf("R1.2 validation passed: %d/%d checks\n", sum(validation$pass), nrow(validation)))
