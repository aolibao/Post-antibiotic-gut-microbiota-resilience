#!/usr/bin/env Rscript

# Reviewer 1.1: cross-definition HDR/LDR sensitivity analysis
#
# All statistics and every visual artefact produced by this script use R.
# The primary branch is the methods-faithful rarefied dataset already audited
# in reviewer_response/03_tables.  D0 preserves the manuscript membership;
# D1-D6 independently rank the 60 exposed participants and compare the lowest
# versus highest 15 recovery-metric values.  Membership overlap is deliberately
# not calculated because it is not the inferential target of R1.1.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(scales)
  library(ragg)
  library(nlme)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20261012)

cli_args <- commandArgs(trailingOnly = TRUE)
unknown_args <- setdiff(cli_args, "--groups-only")
if (length(unknown_args) > 0L) {
  stop("Unknown command-line argument(s): ", paste(unknown_args, collapse = ", "))
}
groups_only <- "--groups-only" %in% cli_args

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
}
out_dir <- file.path(project_root, "reviewer_response", "R1.1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subject_path <- file.path(
  project_root, "reviewer_response", "03_tables",
  "core_subject_recovery_metrics.csv"
)
sample_path <- file.path(
  project_root, "reviewer_response", "03_tables",
  "core_sample_longitudinal_metrics.csv"
)

stopifnot(file.exists(subject_path), file.exists(sample_path))

subject_all <- readr::read_csv(subject_path, show_col_types = FALSE)
sample_all <- readr::read_csv(sample_path, show_col_types = FALSE)

subject <- subject_all |>
  filter(dataset_branch == "rarefied_primary") |>
  arrange(subject)
sample <- sample_all |>
  filter(dataset_branch == "rarefied_primary") |>
  arrange(subject, canons)

stopifnot(nrow(subject) == 60L)
stopifnot(n_distinct(subject$subject) == 60L)
stopifnot(nrow(sample) == 960L)
stopifnot(all(table(sample$subject) == 16L))
stopifnot(all(table(sample$canon) == 60L))
stopifnot(sum(subject$original_endpoint_group == "HDR") == 15L)
stopifnot(sum(subject$original_endpoint_group == "LDR") == 15L)

stage_order <- c(
  paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6)
)
sample$canon <- factor(sample$canon, levels = stage_order, ordered = TRUE)

definition_spec <- tibble::tribble(
  ~definition_id, ~metric_column, ~definition_short, ~definition_full,
  "D0", "grouping_endpoint_shannon", "Endpoint Shannon (reference)",
  "Original mean C4-C6 Shannon diversity; manuscript membership preserved",
  "D1", "endpoint_abs_shannon_deviation", "Endpoint abs. Shannon deviation",
  "Absolute C4-C6 Shannon deviation from the personal A1-A4 mean",
  "D2", "endpoint_mean_bray_to_baseline_centroid", "Endpoint Bray-Curtis",
  "Mean C4-C6 Bray-Curtis distance from the personal A1-A4 centroid",
  "D3", "endpoint_mean_jaccard_to_baseline_consensus", "Endpoint Jaccard",
  "Mean C4-C6 Jaccard distance from the personal A1-A4 consensus",
  "D4", "recovery_auc_abs_shannon_actual_day", "AUC abs. Shannon displacement",
  "Actual-day normalized AUC of absolute Shannon displacement",
  "D5", "recovery_auc_bray_actual_day", "AUC Bray-Curtis",
  "Actual-day normalized AUC of Bray-Curtis displacement",
  "D6", "recovery_auc_jaccard_actual_day", "AUC Jaccard",
  "Actual-day normalized AUC of Jaccard displacement"
) |>
  mutate(
    direction = if_else(definition_id == "D0", "larger", "smaller"),
    definition_display = paste0(definition_id, "  ", definition_short),
    definition_order = row_number()
  )

make_definition_manifest <- function(spec_row) {
  id <- spec_row$definition_id[[1]]
  metric_column <- spec_row$metric_column[[1]]

  if (id == "D0") {
    metric_value <- subject$grouping_endpoint_shannon
    metric_rank <- subject$original_endpoint_rank
    recovery_group <- subject$original_endpoint_group
  } else {
    metric_value <- subject[[metric_column]]
    stopifnot(!anyNA(metric_value))
    ordering <- order(metric_value, subject$subject, na.last = NA)
    recovery_group <- rep("Middle", nrow(subject))
    recovery_group[ordering[seq_len(15)]] <- "HDR"
    recovery_group[ordering[(nrow(subject) - 14):nrow(subject)]] <- "LDR"
    metric_rank <- rank(metric_value, ties.method = "first")
  }

  tibble(
    subject = subject$subject,
    definition_id = id,
    definition_short = spec_row$definition_short[[1]],
    definition_full = spec_row$definition_full[[1]],
    metric_column = metric_column,
    high_recovery_direction = spec_row$direction[[1]],
    metric_value = metric_value,
    metric_rank_ascending = metric_rank,
    recovery_group = recovery_group,
    n_high = sum(recovery_group == "HDR"),
    n_low = sum(recovery_group == "LDR"),
    tie_at_boundary = FALSE,
    dataset_branch = "rarefied_primary"
  )
}

membership <- bind_rows(lapply(
  seq_len(nrow(definition_spec)),
  function(i) make_definition_manifest(definition_spec[i, ])
)) |>
  mutate(
    recovery_group = factor(recovery_group, levels = c("LDR", "Middle", "HDR"))
  )

stopifnot(all(
  membership |>
    group_by(definition_id) |>
    summarise(
      n_hdr = sum(recovery_group == "HDR"),
      n_ldr = sum(recovery_group == "LDR"),
      .groups = "drop"
    ) |>
    transmute(ok = n_hdr == 15 & n_ldr == 15) |>
    pull(ok)
))

readr::write_csv(
  membership |>
    arrange(definition_id, metric_rank_ascending, subject),
  file.path(out_dir, "R1_1_group_manifest_D0_D6.csv"),
  na = ""
)
readr::write_csv(
  definition_spec |>
    select(
      definition_id, definition_short, definition_full, metric_column,
      high_recovery_direction = direction
    ) |>
    mutate(
      group_rule = if_else(
        definition_id == "D0",
        "Preserve manuscript HDR/LDR membership (15/15)",
        "Rank all 60; lowest 15 = HDR and highest 15 = LDR"
      ),
      downstream_target = "Within-definition HDR versus LDR phenotype contrasts"
    ),
  file.path(out_dir, "R1_1_definition_table.csv")
)

if (groups_only) {
  message(
    "PASS: wrote D0-D6 group manifest and definition table only; ",
    "no cross-definition summary figure was generated."
  )
  quit(save = "no", status = 0L)
}

metric_spec <- tibble::tribble(
  ~metric, ~metric_display, ~orientation, ~figure_alignment,
  "shannon", "Shannon", 1, "Fig. 2D / Fig. 4B",
  "richness", "Richness", 1, "Fig. 2C / Fig. 4B",
  "simpson", "Simpson", 1, "Fig. 4B",
  "bray_to_personal_baseline_centroid", "Bray similarity (-distance)", -1,
  "Observed perturbation/recovery phenotype; not MDSINE2",
  "jaccard_to_personal_baseline_consensus", "Jaccard similarity (-distance)", -1,
  "Observed perturbation/recovery phenotype; not MDSINE2"
) |>
  mutate(metric_order = row_number())

sample_analysis <- sample |>
  select(
    ids, subject, canon, day, shannon, richness, simpson,
    bray_to_personal_baseline_centroid,
    jaccard_to_personal_baseline_consensus
  ) |>
  pivot_longer(
    cols = all_of(metric_spec$metric),
    names_to = "metric", values_to = "value"
  ) |>
  left_join(metric_spec, by = "metric") |>
  inner_join(
    membership |>
      filter(recovery_group %in% c("HDR", "LDR")) |>
      select(subject, definition_id, recovery_group),
    by = "subject",
    relationship = "many-to-many"
  ) |>
  left_join(
    definition_spec |>
      select(definition_id, definition_short, definition_display, definition_order),
    by = "definition_id"
  ) |>
  mutate(
    recovery_group = factor(as.character(recovery_group), levels = c("LDR", "HDR")),
    canon = factor(as.character(canon), levels = stage_order, ordered = TRUE)
  )

bootstrap_standardized_difference <- function(value, group, orientation, B = 2000L) {
  keep <- is.finite(value) & group %in% c("HDR", "LDR")
  value <- value[keep]
  group <- as.character(group[keep])
  hdr <- value[group == "HDR"]
  ldr <- value[group == "LDR"]
  global_sd <- stats::sd(value)
  if (!is.finite(global_sd) || global_sd == 0 || length(hdr) < 2 || length(ldr) < 2) {
    return(c(estimate_raw = NA, standardized_effect = NA, ci_low = NA, ci_high = NA,
             p_value = NA, n_high = length(hdr), n_low = length(ldr)))
  }
  estimate_raw <- mean(hdr) - mean(ldr)
  standardized_effect <- orientation * estimate_raw / global_sd
  boot_hdr <- rowMeans(matrix(
    sample(hdr, B * length(hdr), replace = TRUE), nrow = B, ncol = length(hdr)
  ))
  boot_ldr <- rowMeans(matrix(
    sample(ldr, B * length(ldr), replace = TRUE), nrow = B, ncol = length(ldr)
  ))
  boot_effect <- orientation * (boot_hdr - boot_ldr) / global_sd
  ci <- unname(stats::quantile(boot_effect, c(0.025, 0.975), na.rm = TRUE))
  p <- suppressWarnings(stats::wilcox.test(hdr, ldr, exact = FALSE)$p.value)
  c(
    estimate_raw = estimate_raw,
    standardized_effect = standardized_effect,
    ci_low = ci[[1]], ci_high = ci[[2]], p_value = p,
    n_high = length(hdr), n_low = length(ldr)
  )
}

stage_keys <- sample_analysis |>
  distinct(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment, canon
  ) |>
  arrange(definition_order, metric_order, canon)

stage_results <- vector("list", nrow(stage_keys))
for (i in seq_len(nrow(stage_keys))) {
  key <- stage_keys[i, ]
  dat <- sample_analysis |>
    filter(
      definition_id == key$definition_id,
      metric == key$metric,
      canon == key$canon
    )
  set.seed(20261012 + i)
  stat <- bootstrap_standardized_difference(
    dat$value, dat$recovery_group, key$orientation, B = 2000L
  )
  stage_results[[i]] <- bind_cols(key, tibble::as_tibble_row(stat))
}
stage_results <- bind_rows(stage_results) |>
  group_by(definition_id, metric) |>
  mutate(q_within_definition_metric = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(q_across_stage_grid = p.adjust(p_value, method = "BH")) |>
  arrange(definition_order, metric_order, canon)

phase_lookup <- tibble(
  canon = factor(stage_order, levels = stage_order, ordered = TRUE),
  phase = factor(
    c(rep("Baseline", 4), rep("Exposure", 5), "Early recovery", rep("Late recovery", 6)),
    levels = c("Baseline", "Exposure", "Early recovery", "Late recovery")
  )
)

# The three reviewer-facing summaries use A1-A4, B1-B5 and C4-C6.
phase_sample <- sample_analysis |>
  mutate(
    phase = case_when(
      canon %in% paste0("A", 1:4) ~ "Baseline",
      canon %in% paste0("B", 1:5) ~ "Exposure",
      canon %in% paste0("C", 4:6) ~ "Late recovery",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(phase)) |>
  group_by(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment,
    recovery_group, subject, phase
  ) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    phase = factor(phase, levels = c("Baseline", "Exposure", "Late recovery"))
  )

phase_keys <- phase_sample |>
  distinct(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment, phase
  ) |>
  arrange(definition_order, phase, metric_order)

phase_results <- vector("list", nrow(phase_keys))
for (i in seq_len(nrow(phase_keys))) {
  key <- phase_keys[i, ]
  dat <- phase_sample |>
    filter(
      definition_id == key$definition_id,
      metric == key$metric,
      phase == key$phase
    )
  set.seed(20271012 + i)
  stat <- bootstrap_standardized_difference(
    dat$value, dat$recovery_group, key$orientation, B = 5000L
  )
  phase_results[[i]] <- bind_cols(key, tibble::as_tibble_row(stat))
}
phase_results <- bind_rows(phase_results) |>
  group_by(definition_id) |>
  mutate(q_within_definition = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(q_across_phase_grid = p.adjust(p_value, method = "BH")) |>
  arrange(definition_order, phase, metric_order)

# Random-intercept trajectory tests.  These test the joint group-by-time term;
# the stage-specific bootstrap estimates above remain the plotted estimands.
trajectory_tests <- list()
test_index <- 0L
for (id in definition_spec$definition_id) {
  for (metric_name in metric_spec$metric) {
    dat <- sample_analysis |>
      filter(definition_id == id, metric == metric_name) |>
      mutate(
        subject = factor(subject),
        recovery_group = relevel(factor(recovery_group), ref = "LDR"),
        canon = factor(canon, levels = stage_order)
      )
    test_index <- test_index + 1L
    fit_full <- tryCatch(
      nlme::lme(
        value ~ recovery_group * canon,
        random = ~1 | subject,
        data = dat, method = "ML", na.action = na.omit,
        control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 200)
      ),
      error = function(e) NULL
    )
    fit_add <- tryCatch(
      nlme::lme(
        value ~ recovery_group + canon,
        random = ~1 | subject,
        data = dat, method = "ML", na.action = na.omit,
        control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 200)
      ),
      error = function(e) NULL
    )
    fit_time <- tryCatch(
      nlme::lme(
        value ~ canon,
        random = ~1 | subject,
        data = dat, method = "ML", na.action = na.omit,
        control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 200)
      ),
      error = function(e) NULL
    )
    if (is.null(fit_full) || is.null(fit_add) || is.null(fit_time)) {
      p_interaction <- NA_real_
      lrt_interaction <- NA_real_
      df_interaction <- NA_real_
      p_any_group <- NA_real_
      lrt_any_group <- NA_real_
      df_any_group <- NA_real_
      status <- "model_failed"
    } else {
      comparison_interaction <- anova(fit_add, fit_full)
      p_interaction <- comparison_interaction$`p-value`[[2]]
      lrt_interaction <- comparison_interaction$L.Ratio[[2]]
      df_interaction <- comparison_interaction$df[[2]] - comparison_interaction$df[[1]]
      comparison_any_group <- anova(fit_time, fit_full)
      p_any_group <- comparison_any_group$`p-value`[[2]]
      lrt_any_group <- comparison_any_group$L.Ratio[[2]]
      df_any_group <- comparison_any_group$df[[2]] - comparison_any_group$df[[1]]
      status <- "complete"
    }
    trajectory_tests[[test_index]] <- tibble(
      definition_id = id,
      metric = metric_name,
      n_subjects = n_distinct(dat$subject),
      n_observations = nrow(dat),
      lrt_any_group_related = lrt_any_group,
      df_any_group_related = df_any_group,
      p_any_group_related = p_any_group,
      lrt_group_by_time = lrt_interaction,
      df_group_by_time = df_interaction,
      p_group_by_time = p_interaction,
      status = status
    )
  }
}
trajectory_tests <- bind_rows(trajectory_tests) |>
  group_by(definition_id) |>
  mutate(
    q_any_group_within_definition = p.adjust(p_any_group_related, method = "BH"),
    q_interaction_within_definition = p.adjust(p_group_by_time, method = "BH")
  ) |>
  ungroup() |>
  mutate(
    q_any_group_across_models = p.adjust(p_any_group_related, method = "BH"),
    q_interaction_across_models = p.adjust(p_group_by_time, method = "BH")
  ) |>
  left_join(metric_spec, by = "metric") |>
  left_join(
    definition_spec |>
      select(definition_id, definition_short, definition_order),
    by = "definition_id"
  ) |>
  arrange(definition_order, metric_order)

# Fig. 2E-aligned within-group baseline-to-endpoint correlations.  These are
# descriptive phenotype characterisations, not the cross-definition endpoint.
correlation_subject <- sample |>
  mutate(
    period = case_when(
      canon %in% paste0("A", 1:4) ~ "baseline",
      canon %in% paste0("C", 4:6) ~ "endpoint",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period)) |>
  group_by(subject, period) |>
  summarise(shannon = mean(shannon), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = shannon) |>
  inner_join(
    membership |>
      filter(recovery_group %in% c("HDR", "LDR")) |>
      select(subject, definition_id, recovery_group),
    by = "subject",
    relationship = "many-to-many"
  )

bootstrap_spearman <- function(x, y, B = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  estimate <- suppressWarnings(cor(x, y, method = "spearman"))
  p_value <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  idx <- matrix(sample.int(n, B * n, replace = TRUE), nrow = B, ncol = n)
  boot <- apply(idx, 1, function(ii) suppressWarnings(cor(x[ii], y[ii], method = "spearman")))
  ci <- unname(quantile(boot[is.finite(boot)], c(0.025, 0.975), na.rm = TRUE))
  c(rho = estimate, ci_low = ci[[1]], ci_high = ci[[2]], p_value = p_value, n = n)
}

correlation_results <- list()
idx <- 0L
for (id in definition_spec$definition_id) {
  for (grp in c("HDR", "LDR")) {
    idx <- idx + 1L
    dat <- correlation_subject |>
      filter(definition_id == id, recovery_group == grp)
    set.seed(20281012 + idx)
    stat <- bootstrap_spearman(dat$baseline, dat$endpoint)
    correlation_results[[idx]] <- tibble(
      definition_id = id,
      recovery_group = grp,
      rho = stat[["rho"]],
      ci_low = stat[["ci_low"]],
      ci_high = stat[["ci_high"]],
      p_value = stat[["p_value"]],
      n = as.integer(stat[["n"]])
    )
  }
}
correlation_results <- bind_rows(correlation_results) |>
  group_by(definition_id) |>
  mutate(q_within_definition = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(q_across_grid = p.adjust(p_value, method = "BH")) |>
  left_join(
    definition_spec |>
      select(definition_id, definition_short, definition_order),
    by = "definition_id"
  ) |>
  arrange(definition_order, recovery_group)

# Classify alternative-definition effects against D0 by estimates and CIs.
d0_reference <- phase_results |>
  filter(definition_id == "D0") |>
  select(
    phase, metric, d0_effect = standardized_effect,
    d0_ci_low = ci_low, d0_ci_high = ci_high
  )

consistency <- phase_results |>
  filter(definition_id != "D0") |>
  left_join(d0_reference, by = c("phase", "metric")) |>
  mutate(
    same_direction = sign(standardized_effect) == sign(d0_effect),
    d0_interval_excludes_zero = d0_ci_low > 0 | d0_ci_high < 0,
    interval_wholly_opposite = case_when(
      d0_effect > 0 ~ d0_interval_excludes_zero & ci_high < 0,
      d0_effect < 0 ~ d0_interval_excludes_zero & ci_low > 0,
      TRUE ~ FALSE
    ),
    consistency_class = case_when(
      interval_wholly_opposite ~ "Definition-sensitive",
      same_direction & q_within_definition < 0.05 ~ "Statistically corroborated",
      same_direction ~ "Directionally consistent",
      TRUE ~ "Inconclusive"
    )
  ) |>
  arrange(phase, metric_order, definition_order)

consistency_summary <- consistency |>
  count(phase, metric, metric_display, metric_order, consistency_class, name = "n_definitions") |>
  complete(
    nesting(phase, metric, metric_display, metric_order),
    consistency_class = c(
      "Statistically corroborated", "Directionally consistent",
      "Inconclusive", "Definition-sensitive"
    ),
    fill = list(n_definitions = 0L)
  ) |>
  mutate(
    proportion = n_definitions / 6,
    consistency_class = factor(
      consistency_class,
      levels = c(
        "Statistically corroborated", "Directionally consistent",
        "Inconclusive", "Definition-sensitive"
      )
    )
  ) |>
  arrange(phase, metric_order, consistency_class)

readr::write_csv(stage_results, file.path(out_dir, "R1_1_stagewise_effects.csv"))
readr::write_csv(phase_results, file.path(out_dir, "R1_1_phase_effects.csv"))
readr::write_csv(trajectory_tests, file.path(out_dir, "R1_1_trajectory_LRT.csv"))
readr::write_csv(correlation_results, file.path(out_dir, "R1_1_Fig2E_correlations.csv"))
readr::write_csv(consistency, file.path(out_dir, "R1_1_consistency_classification.csv"))
readr::write_csv(consistency_summary, file.path(out_dir, "R1_1_consistency_summary.csv"))

# Focused dominant-genus sensitivity analysis aligned to Fig. 2F and Fig. 4A/D.
# The supplied taxonomic relative-abundance matrix is a separate source-data
# branch, so these results are explicitly secondary and are not called LEfSe,
# FastSpar, or a rarefied count analysis.
taxa_path <- file.path(project_root, "taxa.csv")
abundance_path <- file.path(project_root, "species.csv")
stopifnot(file.exists(taxa_path), file.exists(abundance_path))

taxa_annotation <- read.csv(taxa_path, check.names = FALSE, row.names = 1)
taxon_abundance <- read.csv(abundance_path, check.names = FALSE, row.names = 1)
stopifnot(identical(rownames(taxa_annotation), rownames(taxon_abundance)))

genus_annotation <- sub("^g__", "", taxa_annotation$g)
sample_ids_taxon <- colnames(taxon_abundance)
taxon_matrix <- as.matrix(taxon_abundance)
storage.mode(taxon_matrix) <- "double"

sum_genus <- function(genus_name) {
  selected <- genus_annotation == genus_name
  stopifnot(any(selected))
  colSums(taxon_matrix[selected, , drop = FALSE], na.rm = TRUE) / 100
}

taxon_sample <- tibble(
  ids = sample_ids_taxon,
  bacteroides_prop = sum_genus("Bacteroides"),
  faecalibacterium_prop = sum_genus("Faecalibacterium")
) |>
  mutate(
    bacteroides_asin = asin(sqrt(pmin(pmax(bacteroides_prop, 0), 1))),
    faecalibacterium_asin = asin(sqrt(pmin(pmax(faecalibacterium_prop, 0), 1))),
    faecalibacterium_bacteroides_logratio = log10(
      (faecalibacterium_prop + 1e-4) / (bacteroides_prop + 1e-4)
    )
  ) |>
  inner_join(sample |> select(ids, subject, canon, day), by = "ids")

stopifnot(nrow(taxon_sample) == 960L)

taxon_metric_spec <- tibble::tribble(
  ~metric, ~metric_display, ~orientation, ~figure_alignment,
  "bacteroides_asin", "Lower Bacteroides", -1, "Fig. 2F / Fig. 4A/D",
  "faecalibacterium_asin", "Higher Faecalibacterium", 1, "Fig. 2F / Fig. 4A/D",
  "faecalibacterium_bacteroides_logratio", "Higher F:B balance", 1, "Fig. 4A/D"
) |>
  mutate(metric_order = row_number())

taxon_analysis <- taxon_sample |>
  pivot_longer(
    cols = all_of(taxon_metric_spec$metric),
    names_to = "metric", values_to = "value"
  ) |>
  left_join(taxon_metric_spec, by = "metric") |>
  inner_join(
    membership |>
      filter(recovery_group %in% c("HDR", "LDR")) |>
      select(subject, definition_id, recovery_group),
    by = "subject",
    relationship = "many-to-many"
  ) |>
  left_join(
    definition_spec |>
      select(definition_id, definition_short, definition_display, definition_order),
    by = "definition_id"
  ) |>
  mutate(
    recovery_group = factor(as.character(recovery_group), levels = c("LDR", "HDR")),
    canon = factor(as.character(canon), levels = stage_order, ordered = TRUE)
  )

taxon_stage_keys <- taxon_analysis |>
  distinct(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment, canon
  ) |>
  arrange(definition_order, metric_order, canon)

taxon_stage_results <- vector("list", nrow(taxon_stage_keys))
for (i in seq_len(nrow(taxon_stage_keys))) {
  key <- taxon_stage_keys[i, ]
  dat <- taxon_analysis |>
    filter(
      definition_id == key$definition_id,
      metric == key$metric,
      canon == key$canon
    )
  set.seed(20291012 + i)
  stat <- bootstrap_standardized_difference(
    dat$value, dat$recovery_group, key$orientation, B = 2000L
  )
  taxon_stage_results[[i]] <- bind_cols(key, tibble::as_tibble_row(stat))
}
taxon_stage_results <- bind_rows(taxon_stage_results) |>
  group_by(definition_id, metric) |>
  mutate(q_within_definition_metric = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(q_across_stage_grid = p.adjust(p_value, method = "BH")) |>
  arrange(definition_order, metric_order, canon)

taxon_phase_sample <- taxon_analysis |>
  mutate(
    phase = case_when(
      canon %in% paste0("A", 1:4) ~ "Baseline",
      canon %in% paste0("B", 1:5) ~ "Exposure",
      canon %in% paste0("C", 4:6) ~ "Late recovery",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(phase)) |>
  group_by(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment,
    recovery_group, subject, phase
  ) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(phase = factor(phase, levels = c("Baseline", "Exposure", "Late recovery")))

taxon_phase_keys <- taxon_phase_sample |>
  distinct(
    definition_id, definition_short, definition_display, definition_order,
    metric, metric_display, metric_order, orientation, figure_alignment, phase
  ) |>
  arrange(definition_order, phase, metric_order)

taxon_phase_results <- vector("list", nrow(taxon_phase_keys))
for (i in seq_len(nrow(taxon_phase_keys))) {
  key <- taxon_phase_keys[i, ]
  dat <- taxon_phase_sample |>
    filter(
      definition_id == key$definition_id,
      metric == key$metric,
      phase == key$phase
    )
  set.seed(20301012 + i)
  stat <- bootstrap_standardized_difference(
    dat$value, dat$recovery_group, key$orientation, B = 5000L
  )
  taxon_phase_results[[i]] <- bind_cols(key, tibble::as_tibble_row(stat))
}
taxon_phase_results <- bind_rows(taxon_phase_results) |>
  group_by(definition_id) |>
  mutate(q_within_definition = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(q_across_phase_grid = p.adjust(p_value, method = "BH")) |>
  arrange(definition_order, phase, metric_order)

taxon_d0_reference <- taxon_phase_results |>
  filter(definition_id == "D0") |>
  select(phase, metric, d0_effect = standardized_effect)

taxon_key_results <- taxon_phase_results |>
  filter(definition_id != "D0") |>
  left_join(taxon_d0_reference, by = c("phase", "metric")) |>
  mutate(same_direction = sign(standardized_effect) == sign(d0_effect)) |>
  group_by(phase, metric, metric_display) |>
  summarise(
    n_same_direction = sum(same_direction),
    n_statistically_supported = sum(same_direction & q_within_definition < 0.05),
    .groups = "drop"
  ) |>
  arrange(
    factor(phase, levels = c("Baseline", "Exposure", "Late recovery")),
    match(metric, taxon_metric_spec$metric)
  )

readr::write_csv(
  taxon_stage_results,
  file.path(out_dir, "R1_1_taxonomic_stagewise_effects.csv")
)
readr::write_csv(
  taxon_phase_results,
  file.path(out_dir, "R1_1_taxonomic_phase_effects.csv")
)
readr::write_csv(
  taxon_key_results,
  file.path(out_dir, "R1_1_taxonomic_key_results.csv")
)

theme_response <- function(base_size = 8.2) {
  theme_minimal(base_family = "sans", base_size = base_size) +
    theme(
      text = element_text(color = "#202124"),
      plot.title = element_text(face = "bold", size = rel(1.18), margin = margin(b = 5)),
      plot.subtitle = element_text(color = "#4B5158", size = rel(0.92), margin = margin(b = 6)),
      plot.caption = element_text(color = "#5F6368", size = rel(0.72), hjust = 0),
      strip.text = element_text(face = "bold", color = "#202124", size = rel(0.94)),
      strip.background = element_rect(fill = "#ECE9E2", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "#E5E2DC"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(3.2, "mm"),
      plot.margin = margin(7, 8, 7, 8)
    )
}

warm <- "#B44A3A"
warm_light <- "#E2A18D"
cool <- "#335C7D"
neutral <- "#F5F1E8"
ink <- "#202124"
muted <- "#697178"

# Panel A: the inferential workflow.  No overlap analysis appears here.
flow <- tibble(
  x = 1:4,
  label = c(
    "Six baseline-referenced\nrecovery metrics",
    "Rank all 60 exposed\nparticipants",
    "New HDR (n = 15)\nvs new LDR (n = 15)",
    "Re-test common\nlongitudinal phenotypes"
  ),
  fill = c("#E7EEF3", "#EEEAE3", "#F3DFD8", "#E2EFE8")
)
p_flow <- ggplot(flow, aes(x = x, y = 0)) +
  geom_segment(
    data = tibble(x = c(1.38, 2.38, 3.38), xend = c(1.62, 2.62, 3.62)),
    aes(x = x, xend = xend, y = 0, yend = 0),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(2.2, "mm"), type = "closed"),
    linewidth = 0.55, color = "#6B7075"
  ) +
  geom_label(
    aes(label = label, fill = fill),
    size = 3.0, linewidth = 0.25, label.r = grid::unit(2.1, "mm"),
    color = ink, lineheight = 0.95, family = "sans"
  ) +
  scale_fill_identity() +
  annotate(
    "text", x = 2.5, y = -0.43,
    label = "Target: HDR-LDR phenotype contrasts within every new definition",
    size = 2.5, color = muted, family = "sans"
  ) +
  coord_cartesian(xlim = c(0.55, 4.45), ylim = c(-0.56, 0.36), clip = "off") +
  labs(title = "Analysis logic") +
  theme_void(base_family = "sans", base_size = 8.2) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, color = ink, margin = margin(b = 2)),
    plot.margin = margin(5, 8, 0, 8)
  )

definition_levels <- rev(definition_spec$definition_display)
metric_levels <- metric_spec$metric_display
phase_levels <- c("Baseline", "Exposure", "Late recovery")

phase_plot_data <- phase_results |>
  mutate(
    definition_display = factor(definition_display, levels = definition_levels),
    metric_display = factor(metric_display, levels = metric_levels),
    phase = factor(as.character(phase), levels = phase_levels),
    effect_label = if_else(
      q_within_definition < 0.05,
      sprintf("%+.2f*", standardized_effect),
      sprintf("%+.2f", standardized_effect)
    ),
    label_color = if_else(abs(standardized_effect) >= 1.15, "white", ink)
  )

p_phase <- ggplot(
  phase_plot_data,
  aes(x = metric_display, y = definition_display, fill = standardized_effect)
) +
  geom_tile(color = "white", linewidth = 0.65, width = 0.98, height = 0.96) +
  geom_tile(
    data = phase_plot_data |> filter(definition_id == "D0"),
    fill = NA, color = "#151515", linewidth = 0.75, width = 0.98, height = 0.96
  ) +
  geom_text(
    aes(label = effect_label, color = label_color,
        fontface = if_else(q_within_definition < 0.05, "bold", "plain")),
    size = 2.25, family = "sans", show.legend = FALSE
  ) +
  facet_grid(. ~ phase, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = cool, mid = neutral, high = warm, midpoint = 0,
    limits = c(-2.6, 2.6), oob = scales::squish,
    name = "Oriented HDR-LDR\neffect (SD units)"
  ) +
  scale_color_identity() +
  scale_x_discrete(position = "top", labels = function(x) gsub("- ", "-", x, fixed = TRUE)) +
  labs(
    title = "HDR-LDR contrasts under each recovery definition",
    subtitle = paste0(
      "Positive values reproduce the original HDR direction: higher alpha diversity or smaller distance to baseline. ",
      "* BH q < 0.05 within each definition across 15 phase-by-metric tests."
    ),
    x = NULL, y = NULL
  ) +
  theme_response(7.9) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 0, vjust = 0, size = 6.5),
    axis.text.y = element_text(size = 6.7, color = ink),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(2.4, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(18, "mm"),
    plot.margin = margin(4, 8, 4, 8)
  )

late_forest_data <- phase_results |>
  filter(phase == "Late recovery") |>
  mutate(
    definition_display = factor(definition_display, levels = definition_levels),
    metric_display = factor(metric_display, levels = metric_levels),
    adjusted_support = q_within_definition < 0.05,
    point_fill = if_else(adjusted_support, warm, "#EEEAE3"),
    point_color = if_else(adjusted_support, "#7F2E25", "#4D5358")
  )

p_late_forest <- ggplot(late_forest_data, aes(y = definition_display)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#5F6368", linetype = 2) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = definition_display),
    linewidth = 0.62, color = "#4F555A", lineend = "round"
  ) +
  geom_point(
    aes(x = standardized_effect, fill = point_fill, color = point_color),
    shape = 21, size = 2.15, stroke = 0.48
  ) +
  facet_grid(. ~ metric_display) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    limits = c(-1.0, 2.65), breaks = c(-1, 0, 1, 2),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    title = "Late-recovery effects and 95% bootstrap intervals",
    subtitle = "Filled points indicate BH q < 0.05 within the corresponding definition.",
    x = "Oriented HDR-LDR effect (SD units)", y = NULL
  ) +
  theme_response(7.7) +
  theme(
    axis.text.y = element_text(size = 6.4),
    axis.text.x = element_text(size = 6.2),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = grid::unit(2.3, "mm"),
    strip.text = element_text(size = 6.5),
    plot.margin = margin(3, 8, 3, 8)
  )

consistency_colors <- c(
  "Statistically corroborated" = warm,
  "Directionally consistent" = warm_light,
  "Inconclusive" = "#C8C5BF",
  "Definition-sensitive" = cool
)

p_consistency <- consistency_summary |>
  mutate(
    phase = factor(as.character(phase), levels = phase_levels),
    metric_display = factor(metric_display, levels = rev(metric_levels))
  ) |>
  ggplot(aes(x = n_definitions, y = metric_display, fill = consistency_class)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 6, linewidth = 0.3, color = "#404040") +
  facet_grid(. ~ phase) +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6.08), expand = c(0, 0)) +
  scale_fill_manual(values = consistency_colors, drop = FALSE, name = NULL) +
  labs(
    title = "How many alternative definitions reproduced each D0 effect?",
    subtitle = "Counts are based on D1-D6 estimates, intervals and multiplicity-adjusted support.",
    x = "Number of alternative definitions (of 6)", y = NULL
  ) +
  theme_response(7.8) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.y = element_text(size = 6.8),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = grid::unit(4, "mm"),
    plot.margin = margin(4, 8, 5, 8)
  )

main_figure <- p_flow / p_phase / p_late_forest / p_consistency +
  plot_layout(heights = c(1.0, 3.25, 2.25, 2.05)) +
  plot_annotation(
    title = "Alternative recovery definitions: cross-figure phenotype sensitivity",
    subtitle = paste0(
      "Each definition independently creates new high- and low-recovery groups. ",
      "The comparison is the phenotype difference within that definition."
    ),
    caption = paste0(
      "Primary rarefied branch; 60 exposed participants; 15 HDR and 15 LDR per definition. ",
      "Effects are mean differences standardized by the 30 included participants; 95% percentile bootstrap CIs use 5,000 resamples. ",
      "The distance panels are observed longitudinal phenotypes, not MDSINE2 estimates."
    ),
    tag_levels = "A",
    theme = theme(
      text = element_text(family = "sans", color = ink),
      plot.title = element_text(face = "bold", size = 12.2, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 8.7, color = "#4B5158", margin = margin(b = 5)),
      plot.caption = element_text(size = 6.3, color = muted, hjust = 0, lineheight = 0.95),
      plot.tag = element_text(face = "bold", size = 11.5)
    )
  )

stage_plot_data <- stage_results |>
  mutate(
    definition_display = factor(definition_display, levels = definition_levels),
    canon = factor(as.character(canon), levels = stage_order),
    sig = q_within_definition_metric < 0.05
  )

stage_metric_plot <- function(metric_name, title_text, show_x = TRUE) {
  dat <- stage_plot_data |> filter(metric == metric_name)
  ggplot(dat, aes(x = canon, y = definition_display, fill = standardized_effect)) +
    geom_vline(xintercept = c(4.5, 10.5), linewidth = 0.45, color = "#777777") +
    geom_tile(color = "white", linewidth = 0.43, width = 0.98, height = 0.92) +
    geom_tile(
      data = dat |> filter(definition_id == "D0"),
      fill = NA, color = "#151515", linewidth = 0.6, width = 0.98, height = 0.92
    ) +
    geom_point(
      data = dat |> filter(sig),
      shape = 21, size = 1.25, stroke = 0.25,
      fill = "white", color = "#242424"
    ) +
    scale_fill_gradient2(
      low = cool, mid = neutral, high = warm, midpoint = 0,
      limits = c(-2.6, 2.6), oob = scales::squish,
      name = "Oriented HDR-LDR\neffect (SD units)"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.25, 0.25))) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_response(7.7) +
    theme(
      plot.title = element_text(size = 8.9, margin = margin(b = 4)),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6.3),
      axis.text.x = if (show_x) element_text(size = 6.4, color = ink) else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      legend.position = "none",
      plot.margin = margin(7, 7, 1, 7)
    )
}

p_s1 <- stage_metric_plot("shannon", "Fig. 2D / Fig. 4B | Shannon", FALSE)
p_s2 <- stage_metric_plot("richness", "Fig. 2C / Fig. 4B | Richness", FALSE)
p_s3 <- stage_metric_plot("simpson", "Fig. 4B | Simpson", FALSE)
p_s4 <- stage_metric_plot(
  "bray_to_personal_baseline_centroid",
  "Observed perturbation/recovery phenotype | Bray similarity (-distance; not MDSINE2)",
  FALSE
)
p_s5 <- stage_metric_plot(
  "jaccard_to_personal_baseline_consensus",
  "Observed perturbation/recovery phenotype | Jaccard similarity (-distance; not MDSINE2)",
  TRUE
)

legend_plot <- ggplot(
  tibble(x = seq(-2.6, 2.6, length.out = 100), y = 1),
  aes(x = x, y = y, fill = x)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = cool, mid = neutral, high = warm, midpoint = 0,
    limits = c(-2.6, 2.6), name = "Oriented HDR-LDR effect (SD units)"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", barwidth = grid::unit(55, "mm"),
    barheight = grid::unit(3.3, "mm")
  )) +
  theme_void(base_family = "sans") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold", size = 7.2),
        legend.text = element_text(size = 6.5))

stage_figure <- p_s1 / p_s2 / p_s3 / p_s4 / p_s5 / legend_plot +
  plot_layout(heights = c(1, 1, 1, 1, 1.08, 0.34)) +
  plot_annotation(
    title = "Sixteen-stage HDR-LDR effects under alternative recovery definitions",
    subtitle = paste0(
      "Positive values indicate the original HDR-like direction. Open circles mark BH q < 0.05 ",
      "within a definition-metric family across 16 stages."
    ),
    caption = paste0(
      "D0 is outlined and preserves manuscript membership. D1-D6 are newly ranked groups. ",
      "Observed distance trajectories provide an empirical perturbation phenotype only; they do not replace formal MDSINE2 reruns."
    ),
    theme = theme(
      text = element_text(family = "sans", color = ink),
      plot.title = element_text(face = "bold", size = 12, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 8.5, color = "#4B5158", margin = margin(b = 3)),
      plot.caption = element_text(size = 6.3, color = muted, hjust = 0)
    )
  )

taxon_phase_plot_data <- taxon_phase_results |>
  mutate(
    definition_display = factor(definition_display, levels = definition_levels),
    metric_display = factor(metric_display, levels = taxon_metric_spec$metric_display),
    phase = factor(as.character(phase), levels = phase_levels),
    effect_label = if_else(
      q_within_definition < 0.05,
      sprintf("%+.2f*", standardized_effect),
      sprintf("%+.2f", standardized_effect)
    ),
    label_color = if_else(abs(standardized_effect) >= 1.15, "white", ink)
  )

p_taxon_phase <- ggplot(
  taxon_phase_plot_data,
  aes(x = metric_display, y = definition_display, fill = standardized_effect)
) +
  geom_tile(color = "white", linewidth = 0.65, width = 0.98, height = 0.96) +
  geom_tile(
    data = taxon_phase_plot_data |> filter(definition_id == "D0"),
    fill = NA, color = "#151515", linewidth = 0.75, width = 0.98, height = 0.96
  ) +
  geom_text(
    aes(label = effect_label, color = label_color,
        fontface = if_else(q_within_definition < 0.05, "bold", "plain")),
    size = 2.35, family = "sans", show.legend = FALSE
  ) +
  facet_grid(. ~ phase, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = cool, mid = neutral, high = warm, midpoint = 0,
    limits = c(-2.6, 2.6), oob = scales::squish,
    name = "Oriented HDR-LDR\neffect (SD units)"
  ) +
  scale_color_identity() +
  scale_x_discrete(position = "top") +
  labs(
    title = "Dominant-genus effects across recovery definitions",
    subtitle = paste0(
      "Positive values indicate the original HDR-like direction. ",
      "* BH q < 0.05 within each definition across nine phase-by-metric tests."
    ),
    x = NULL, y = NULL
  ) +
  theme_response(7.9) +
  theme(
    axis.text.x = element_text(angle = 34, hjust = 0, vjust = 0, size = 6.6),
    axis.text.y = element_text(size = 6.7),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(2.4, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(20, "mm"),
    plot.margin = margin(4, 8, 3, 8)
  )

taxon_stage_plot_data <- taxon_stage_results |>
  mutate(
    definition_display = factor(definition_display, levels = definition_levels),
    canon = factor(as.character(canon), levels = stage_order),
    sig = q_within_definition_metric < 0.05
  )

taxon_stage_metric_plot <- function(metric_name, title_text, show_x = TRUE) {
  dat <- taxon_stage_plot_data |> filter(metric == metric_name)
  ggplot(dat, aes(x = canon, y = definition_display, fill = standardized_effect)) +
    geom_vline(xintercept = c(4.5, 10.5), linewidth = 0.45, color = "#777777") +
    geom_tile(color = "white", linewidth = 0.43, width = 0.98, height = 0.92) +
    geom_tile(
      data = dat |> filter(definition_id == "D0"),
      fill = NA, color = "#151515", linewidth = 0.6, width = 0.98, height = 0.92
    ) +
    geom_point(
      data = dat |> filter(sig),
      shape = 21, size = 1.25, stroke = 0.25,
      fill = "white", color = "#242424"
    ) +
    scale_fill_gradient2(
      low = cool, mid = neutral, high = warm, midpoint = 0,
      limits = c(-2.6, 2.6), oob = scales::squish,
      name = "Oriented HDR-LDR effect"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.25, 0.25))) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_response(7.7) +
    theme(
      plot.title = element_text(size = 8.8, margin = margin(b = 3)),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6.3),
      axis.text.x = if (show_x) element_text(size = 6.4) else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
      legend.position = "none",
      plot.margin = margin(3, 7, 1, 7)
    )
}

p_tax_s1 <- taxon_stage_metric_plot(
  "bacteroides_asin", "Fig. 2F / Fig. 4A-D | Lower Bacteroides", FALSE
)
p_tax_s2 <- taxon_stage_metric_plot(
  "faecalibacterium_asin", "Fig. 2F / Fig. 4A-D | Higher Faecalibacterium", FALSE
)
p_tax_s3 <- taxon_stage_metric_plot(
  "faecalibacterium_bacteroides_logratio",
  "Fig. 4A-D | Higher Faecalibacterium:Bacteroides balance", TRUE
)

taxon_figure <- p_taxon_phase / p_tax_s1 / p_tax_s2 / p_tax_s3 / legend_plot +
  plot_layout(heights = c(2.65, 1, 1, 1.08, 0.34)) +
  plot_annotation(
    title = "Dominant-genus sensitivity under alternative recovery definitions",
    subtitle = paste0(
      "Focused taxonomic readouts test whether the Bacteroides/Faecalibacterium pattern ",
      "recurs within each newly defined HDR-LDR contrast."
    ),
    caption = paste0(
      "Secondary source-data branch from the supplied species relative-abundance matrix. ",
      "Genus proportions use arcsine-square-root transformation; F:B uses log10((F+10^-4)/(B+10^-4)). ",
      "This focused sensitivity analysis is not a LEfSe or FastSpar rerun."
    ),
    theme = theme(
      text = element_text(family = "sans", color = ink),
      plot.title = element_text(face = "bold", size = 12, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 8.5, color = "#4B5158", margin = margin(b = 3)),
      plot.caption = element_text(size = 6.3, color = muted, hjust = 0)
    )
  )

render_plot <- function(plot_object, stem, width_in, height_in) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  tif_path <- file.path(out_dir, paste0(stem, ".tif"))
  png_path <- file.path(out_dir, paste0(stem, "_preview.png"))

  grDevices::pdf(
    file = pdf_path, width = width_in, height = height_in,
    family = "Helvetica", onefile = TRUE, bg = "white", useDingbats = FALSE
  )
  print(plot_object)
  grDevices::dev.off()

  ragg::agg_tiff(
    filename = tif_path, width = width_in, height = height_in, units = "in",
    res = 600, compression = "lzw", background = "white",
    scaling = 1, bitsize = 8
  )
  print(plot_object)
  grDevices::dev.off()

  ragg::agg_png(
    filename = png_path, width = width_in, height = height_in, units = "in",
    res = 300, background = "white", scaling = 1
  )
  print(plot_object)
  grDevices::dev.off()

  tibble(
    stem = stem,
    format = c("PDF", "TIFF", "PNG preview"),
    path = c(pdf_path, tif_path, png_path),
    width_in = width_in,
    height_in = height_in,
    resolution_dpi = c(NA_integer_, 600L, 300L),
    size_bytes = file.info(c(pdf_path, tif_path, png_path))$size
  )
}

render_manifest <- bind_rows(
  render_plot(main_figure, "R1_1_Fig2_4_cross_definition", 7.205, 10.65),
  render_plot(stage_figure, "R1_1_Fig2_4_stagewise", 7.205, 10.4),
  render_plot(taxon_figure, "R1_1_taxonomic_sensitivity", 7.205, 9.2)
)
readr::write_csv(render_manifest, file.path(out_dir, "R1_1_render_manifest.csv"))

# R-only numerical visual QA on the actual PNG previews.
qa_rows <- list()
if (requireNamespace("png", quietly = TRUE)) {
  for (i in seq_len(nrow(render_manifest |> filter(format == "PNG preview")))) {
    row <- (render_manifest |> filter(format == "PNG preview"))[i, ]
    img <- png::readPNG(row$path)
    rgb <- img[, , seq_len(min(3, dim(img)[3])), drop = FALSE]
    ink_mask <- apply(rgb, c(1, 2), function(z) any(z < 0.97))
    edge_fraction <- c(
      top = mean(ink_mask[1:3, ]),
      bottom = mean(ink_mask[(nrow(ink_mask) - 2):nrow(ink_mask), ]),
      left = mean(ink_mask[, 1:3]),
      right = mean(ink_mask[, (ncol(ink_mask) - 2):ncol(ink_mask)])
    )
    qa_rows[[length(qa_rows) + 1L]] <- tibble(
      figure = row$stem,
      width_px = dim(img)[2],
      height_px = dim(img)[1],
      top_edge_ink_fraction = edge_fraction[["top"]],
      bottom_edge_ink_fraction = edge_fraction[["bottom"]],
      left_edge_ink_fraction = edge_fraction[["left"]],
      right_edge_ink_fraction = edge_fraction[["right"]],
      edge_clipping_flag = any(edge_fraction > 0.10),
      preview_status = if_else(any(edge_fraction > 0.10), "REVIEW", "PASS")
    )
  }
}
visual_qa <- bind_rows(qa_rows)
readr::write_csv(visual_qa, file.path(out_dir, "R1_1_visual_QA.csv"))

input_hashes <- tibble(
  file = c(subject_path, sample_path, taxa_path, abundance_path),
  sha256 = c(
    digest::digest(file = subject_path, algo = "sha256"),
    digest::digest(file = sample_path, algo = "sha256"),
    digest::digest(file = taxa_path, algo = "sha256"),
    digest::digest(file = abundance_path, algo = "sha256")
  )
)

analysis_manifest <- list(
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = "R",
  r_version = R.version.string,
  seed = 20261012,
  dataset_branch = "rarefied_primary",
  biological_unit = "participant",
  n_exposed = 60,
  n_high_per_definition = 15,
  n_low_per_definition = 15,
  definitions = definition_spec$definition_id,
  excluded_from_matched_grid = list(
    D7 = "Exploratory composite has a boundary tie; omitted from the confirmatory D0-D6 grid.",
    TTR = "Sustained-return time is right-censored and cannot support a matched 15/15 split without a prespecified censoring rule."
  ),
  estimand = paste0(
    "HDR minus LDR mean difference, standardized by the SD among the 30 compared participants; ",
    "distance outcomes multiplied by -1 so positive values consistently indicate the original HDR-like direction."
  ),
  confidence_intervals = "Participant-level percentile bootstrap; 2,000 stagewise and 5,000 phase-level resamples.",
  tests = list(
    cell_level = "Two-sided Wilcoxon rank-sum test.",
    longitudinal = paste0(
      "Random-intercept nlme model; likelihood-ratio tests of any group-related term ",
      "and of the group-by-categorical-time interaction."
    ),
    multiplicity = "BH within the prespecified panel family and BH across the complete sensitivity grid."
  ),
  evidence_boundary = list(
    complete = c(
      "Fig. 2/4 alpha-diversity trajectories",
      "Bray-Curtis and Jaccard observed distance trajectories",
      "Focused Bacteroides/Faecalibacterium sensitivity analysis"
    ),
    descriptive = "Fig. 2E baseline-endpoint within-group Spearman correlations.",
    not_claimed = c(
      "Alternative-definition MDSINE2 posterior growth, perturbation and interaction results",
      "Full alternative-definition LEfSe results",
      "Alternative-definition FastSpar networks"
    )
  ),
  input_hashes = split(input_hashes$sha256, input_hashes$file),
  outputs = render_manifest$path
)
jsonlite::write_json(
  analysis_manifest,
  file.path(out_dir, "R1_1_analysis_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)

capture.output(sessionInfo(), file = file.path(out_dir, "R1_1_R_sessionInfo.txt"))

# Machine-readable key results used to draft the reviewer response.
key_results <- consistency |>
  group_by(phase, metric, metric_display) |>
  summarise(
    n_same_direction = sum(same_direction),
    n_statistically_corroborated = sum(consistency_class == "Statistically corroborated"),
    n_definition_sensitive = sum(consistency_class == "Definition-sensitive"),
    .groups = "drop"
  ) |>
  arrange(factor(phase, levels = phase_levels), match(metric, metric_spec$metric))
readr::write_csv(key_results, file.path(out_dir, "R1_1_key_results.csv"))

cat("R1.1 R reanalysis complete.\n")
cat("Output directory:", out_dir, "\n")
cat("Main figure:", file.path(out_dir, "R1_1_Fig2_4_cross_definition.pdf"), "\n")
cat("Stagewise figure:", file.path(out_dir, "R1_1_Fig2_4_stagewise.pdf"), "\n")
cat("Taxonomic sensitivity:", file.path(out_dir, "R1_1_taxonomic_sensitivity.pdf"), "\n")
