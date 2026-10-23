#!/usr/bin/env Rscript

# Reviewer 1.2: full-cohort continuous analysis, cutoff sensitivity, and
# sequencing-depth/rarefaction robustness.  All inferential cross-sectional
# models use the participant as the sampling unit.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(digest)
})

SEED <- 20261012L
BOOTSTRAP_R <- 5000L
ICC_BOOTSTRAP_R <- 2000L
PRIMARY_BOOTSTRAP_R <- 10000L
PRIMARY_PERMUTATIONS <- 9999L
RAREFACTION_DEPTH <- 2693000L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve this script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
analysis_dir <- dirname(script_path)
root_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")
cache_dir <- file.path(analysis_dir, "cache")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

save_embedded_pdf <- function(filename, plot, width_mm, height_mm) {
  raw_pdf <- tempfile(pattern = "R1_2_raw_", fileext = ".pdf")
  embedded_pdf <- tempfile(pattern = "R1_2_embedded_", fileext = ".pdf")
  on.exit(unlink(c(raw_pdf, embedded_pdf)), add = TRUE)
  ggsave(
    raw_pdf, plot, device = "pdf", width = width_mm, height = height_mm,
    units = "mm", bg = "white", family = "sans", useDingbats = FALSE
  )
  grDevices::embedFonts(raw_pdf, outfile = embedded_pdf)
  if (!file.copy(embedded_pdf, filename, overwrite = TRUE)) {
    stop("Unable to install embedded-font PDF: ", filename)
  }
}

subject_path <- file.path(root_dir, "reviewer_response", "03_tables", "core_subject_recovery_metrics.csv")
sample_path <- file.path(root_dir, "reviewer_response", "03_tables", "core_sample_longitudinal_metrics.csv")
rarefaction_path <- file.path(cache_dir, "R1_2_repeated_rarefaction_subject_metrics.csv")
required_inputs <- c(subject_path, sample_path, rarefaction_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "),
       ". Run run_repeated_rarefaction.py before this script.")
}

z_score <- function(x) {
  scale_value <- sd(x, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value == 0) stop("Cannot standardize a constant/non-finite vector")
  as.numeric((x - mean(x, na.rm = TRUE)) / scale_value)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))
}

hc3_vcov <- function(fit) {
  x <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  if (any(leverage >= 1 - 1e-12)) stop("HC3 undefined because leverage approaches one")
  bread <- solve(crossprod(x))
  weights <- residual^2 / (1 - leverage)^2
  meat <- crossprod(x, x * weights)
  unname(bread %*% meat %*% bread)
}

hc3_tidy <- function(fit, conf_level = 0.95) {
  covariance <- hc3_vcov(fit)
  estimate <- coef(fit)
  standard_error <- sqrt(diag(covariance))
  statistic <- estimate / standard_error
  degrees_freedom <- df.residual(fit)
  alpha <- 1 - conf_level
  critical <- qt(1 - alpha / 2, df = degrees_freedom)
  tibble(
    term = names(estimate),
    estimate = as.numeric(estimate),
    HC3_SE = as.numeric(standard_error),
    statistic = as.numeric(statistic),
    df_residual = degrees_freedom,
    p_value = 2 * pt(abs(statistic), df = degrees_freedom, lower.tail = FALSE),
    CI95_low = as.numeric(estimate - critical * standard_error),
    CI95_high = as.numeric(estimate + critical * standard_error)
  )
}

fit_recovery_model <- function(data, outcome_column, metric_id, metric_label, primary = FALSE) {
  analysis <- data |>
    transmute(
      subject,
      baseline_shannon,
      displacement = .data[[outcome_column]],
      recovery_score = -.data[[outcome_column]]
    ) |>
    filter(if_all(everything(), ~ !is.na(.x))) |>
    mutate(
      z_baseline_shannon = z_score(baseline_shannon),
      z_recovery_score = z_score(recovery_score)
    )
  if (nrow(analysis) != n_distinct(analysis$subject)) stop("Participant duplication in model data")
  raw_fit <- lm(recovery_score ~ z_baseline_shannon, data = analysis)
  standardized_fit <- lm(z_recovery_score ~ z_baseline_shannon, data = analysis)
  raw_term <- hc3_tidy(raw_fit) |> filter(term == "z_baseline_shannon")
  standardized_term <- hc3_tidy(standardized_fit) |> filter(term == "z_baseline_shannon")
  tibble(
    metric_id = metric_id,
    metric_label = metric_label,
    outcome_column = outcome_column,
    primary = primary,
    score_definition = paste0("-", outcome_column, " (higher = better recovery)"),
    n_subjects = nrow(analysis),
    raw_score_slope_per_1SD_baseline = raw_term$estimate,
    raw_score_HC3_SE = raw_term$HC3_SE,
    raw_score_CI95_low = raw_term$CI95_low,
    raw_score_CI95_high = raw_term$CI95_high,
    standardized_beta = standardized_term$estimate,
    beta_HC3_SE = standardized_term$HC3_SE,
    beta_statistic = standardized_term$statistic,
    beta_df = standardized_term$df_residual,
    beta_p_value = standardized_term$p_value,
    beta_CI95_low = standardized_term$CI95_low,
    beta_CI95_high = standardized_term$CI95_high,
    R2 = summary(standardized_fit)$r.squared,
    adjusted_R2 = summary(standardized_fit)$adj.r.squared,
    predictor_mean = mean(analysis$baseline_shannon),
    predictor_SD = sd(analysis$baseline_shannon),
    displacement_mean = mean(analysis$displacement),
    displacement_SD = sd(analysis$displacement)
  )
}

hedges_g <- function(high, low) {
  n_high <- length(high)
  n_low <- length(low)
  pooled_sd <- sqrt(((n_high - 1) * var(high) + (n_low - 1) * var(low)) /
                      (n_high + n_low - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  cohen_d <- (mean(high) - mean(low)) / pooled_sd
  correction <- 1 - 3 / (4 * (n_high + n_low) - 9)
  correction * cohen_d
}

bootstrap_ranked_tail_contrast <- function(data, selection_column,
                                           baseline_column, n_tail,
                                           high_recovery_at_low_value,
                                           repetitions, seed) {
  set.seed(seed)
  n_subjects <- nrow(data)
  draws <- replicate(repetitions, {
    index <- sample.int(n_subjects, n_subjects, replace = TRUE)
    bootstrap_data <- data[index, , drop = FALSE]
    ordering <- order(bootstrap_data[[selection_column]], seq_len(n_subjects), na.last = NA)
    if (!high_recovery_at_low_value) ordering <- rev(ordering)
    high_star <- bootstrap_data[[baseline_column]][ordering[seq_len(n_tail)]]
    low_star <- bootstrap_data[[baseline_column]][ordering[(n_subjects - n_tail + 1):n_subjects]]
    c(mean_difference = mean(high_star) - mean(low_star), hedges_g = hedges_g(high_star, low_star))
  })
  tibble(
    statistic = rownames(draws),
    CI95_low = apply(draws, 1, quantile, probs = 0.025, na.rm = TRUE, type = 6),
    CI95_high = apply(draws, 1, quantile, probs = 0.975, na.rm = TRUE, type = 6)
  )
}

icc_absolute_agreement <- function(matrix_values) {
  x <- as.matrix(matrix_values)
  if (any(!is.finite(x))) stop("ICC matrix contains non-finite values")
  n <- nrow(x)
  k <- ncol(x)
  grand <- mean(x)
  row_mean <- rowMeans(x)
  column_mean <- colMeans(x)
  ss_subject <- k * sum((row_mean - grand)^2)
  ss_time <- n * sum((column_mean - grand)^2)
  residual <- sweep(sweep(x, 1, row_mean, "-"), 2, column_mean, "-") + grand
  ss_error <- sum(residual^2)
  ms_subject <- ss_subject / (n - 1)
  ms_time <- ss_time / (k - 1)
  ms_error <- ss_error / ((n - 1) * (k - 1))
  icc_a1 <- (ms_subject - ms_error) /
    (ms_subject + (k - 1) * ms_error + k * (ms_time - ms_error) / n)
  icc_ak <- (ms_subject - ms_error) /
    (ms_subject + (ms_time - ms_error) / n)
  c(ICC_A_1 = icc_a1, ICC_A_k = icc_ak)
}

bootstrap_icc <- function(matrix_values, repetitions, seed) {
  set.seed(seed)
  n <- nrow(matrix_values)
  draws <- replicate(repetitions, {
    index <- sample.int(n, n, replace = TRUE)
    icc_absolute_agreement(matrix_values[index, , drop = FALSE])
  })
  tibble(
    reliability = rownames(draws),
    estimate = icc_absolute_agreement(matrix_values),
    CI95_low = apply(draws, 1, quantile, probs = 0.025, type = 6, na.rm = TRUE),
    CI95_high = apply(draws, 1, quantile, probs = 0.975, type = 6, na.rm = TRUE),
    bootstrap_repetitions = repetitions
  )
}

subject_all <- read_csv(subject_path, show_col_types = FALSE)
sample_all <- read_csv(sample_path, show_col_types = FALSE)
subject <- subject_all |> filter(dataset_branch == "rarefied_primary") |> arrange(subject)
subject_unrarefied <- subject_all |> filter(dataset_branch == "unrarefied_sensitivity") |> arrange(subject)
sample_primary <- sample_all |> filter(dataset_branch == "rarefied_primary") |> arrange(subject, canons)
sample_unrarefied <- sample_all |> filter(dataset_branch == "unrarefied_sensitivity") |> arrange(subject, canons)

if (nrow(subject) != 60L || n_distinct(subject$subject) != 60L) stop("Expected 60 unique exposed participants")
if (nrow(subject_unrarefied) != 60L || !identical(subject$subject, subject_unrarefied$subject)) {
  stop("Rarefied and unrarefied participant grids differ")
}
if (nrow(sample_primary) != 960L || any(table(sample_primary$subject) != 16L)) {
  stop("Expected a complete 60 x 16 rarefied sample grid")
}
if (any(sample_primary$read_depth != RAREFACTION_DEPTH)) stop("Primary sample table is not at fixed depth")

metric_definitions <- tribble(
  ~metric_id, ~metric_label, ~outcome_column, ~primary, ~metric_family,
  "D1", "Endpoint |Shannon deviation|", "endpoint_abs_shannon_deviation", FALSE, "Shannon-derived",
  "D2", "Endpoint Bray-Curtis distance", "endpoint_mean_bray_to_baseline_centroid", FALSE, "Composition-based",
  "D3", "Endpoint Jaccard distance", "endpoint_mean_jaccard_to_baseline_consensus", FALSE, "Composition-based",
  "D4", "Shannon displacement AUC", "recovery_auc_abs_shannon_actual_day", FALSE, "Shannon-derived",
  "D5", "Bray-Curtis displacement AUC", "recovery_auc_bray_actual_day", TRUE, "Composition-based",
  "D6", "Jaccard displacement AUC", "recovery_auc_jaccard_actual_day", FALSE, "Composition-based"
)

missing_metric_columns <- setdiff(metric_definitions$outcome_column, names(subject))
if (length(missing_metric_columns) > 0L) stop("Missing recovery metric columns: ", paste(missing_metric_columns, collapse = ", "))
if (any(!is.finite(as.matrix(subject[, c("baseline_shannon", metric_definitions$outcome_column)])))) {
  stop("Primary subject data contain non-finite predictor/outcome values")
}

model_rows <- lapply(seq_len(nrow(metric_definitions)), function(index) {
  definition <- metric_definitions[index, ]
  fit_recovery_model(
    subject,
    definition$outcome_column,
    definition$metric_id,
    definition$metric_label,
    definition$primary
  )
})
recovery_models <- bind_rows(model_rows) |>
  left_join(metric_definitions |> select(metric_id, metric_family), by = "metric_id") |>
  mutate(
    p_Holm_secondary = NA_real_,
    inference_role = if_else(primary, "Primary", "Secondary")
  )
secondary_index <- which(!recovery_models$primary)
recovery_models$p_Holm_secondary[secondary_index] <- p.adjust(
  recovery_models$beta_p_value[secondary_index], method = "holm"
)
primary_model <- recovery_models |> filter(primary)
if (nrow(primary_model) != 1L || primary_model$metric_id != "D5") stop("D5 must be the unique primary model")

# Reference analysis that directly replaces the original within-tail
# baseline--endpoint correlations.  It is reported as descriptive because
# absolute endpoint Shannon is not itself a recovery metric.
reference_data <- subject |>
  transmute(
    subject,
    baseline_shannon,
    endpoint_shannon,
    endpoint_signed_shannon_change,
    z_baseline_shannon = z_score(baseline_shannon),
    z_endpoint_shannon = z_score(endpoint_shannon),
    z_signed_change = z_score(endpoint_signed_shannon_change)
  )
reference_endpoint_fit <- lm(z_endpoint_shannon ~ z_baseline_shannon, data = reference_data)
reference_change_fit <- lm(z_signed_change ~ z_baseline_shannon, data = reference_data)
reference_endpoint_term <- hc3_tidy(reference_endpoint_fit) |> filter(term == "z_baseline_shannon")
reference_change_term <- hc3_tidy(reference_change_fit) |> filter(term == "z_baseline_shannon")
reference_endpoint_spearman <- suppressWarnings(cor.test(
  reference_data$baseline_shannon, reference_data$endpoint_shannon,
  method = "spearman", exact = FALSE
))
reference_analysis <- bind_rows(
  tibble(
    outcome = "Mean C4-C6 endpoint Shannon (absolute endpoint; not recovery)",
    standardized_beta = reference_endpoint_term$estimate,
    beta_HC3_SE = reference_endpoint_term$HC3_SE,
    CI95_low = reference_endpoint_term$CI95_low,
    CI95_high = reference_endpoint_term$CI95_high,
    p_value = reference_endpoint_term$p_value,
    R2 = summary(reference_endpoint_fit)$r.squared,
    spearman_rho = unname(reference_endpoint_spearman$estimate),
    spearman_p = reference_endpoint_spearman$p.value,
    interpretation_boundary = "Full-cohort replacement for selection-coupled within-tail correlations; endpoint level is not a recovery metric."
  ),
  tibble(
    outcome = "Endpoint signed Shannon change (change score; mathematical coupling noted)",
    standardized_beta = reference_change_term$estimate,
    beta_HC3_SE = reference_change_term$HC3_SE,
    CI95_low = reference_change_term$CI95_low,
    CI95_high = reference_change_term$CI95_high,
    p_value = reference_change_term$p_value,
    R2 = summary(reference_change_fit)$r.squared,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    interpretation_boundary = "Sensitivity only; change scores are mathematically coupled to baseline values."
  )
)

# Nonparametric sensitivity for all prespecified continuous recovery metrics.
spearman_sensitivity <- lapply(seq_len(nrow(metric_definitions)), function(index) {
  definition <- metric_definitions[index, ]
  recovery_score <- -subject[[definition$outcome_column]]
  test <- suppressWarnings(cor.test(subject$baseline_shannon, recovery_score, method = "spearman", exact = FALSE))
  tibble(
    metric_id = definition$metric_id,
    metric_label = definition$metric_label,
    primary = definition$primary,
    n_subjects = length(recovery_score),
    spearman_rho = unname(test$estimate),
    statistic_S = unname(test$statistic),
    p_value = test$p.value
  )
}) |> bind_rows() |>
  mutate(p_Holm_secondary = NA_real_)
spearman_secondary_index <- which(!spearman_sensitivity$primary)
spearman_sensitivity$p_Holm_secondary[spearman_secondary_index] <- p.adjust(
  spearman_sensitivity$p_value[spearman_secondary_index], method = "holm"
)

# Primary-model diagnostics, quadratic sensitivity, and leave-one-participant-out stability.
primary_data <- subject |>
  transmute(
    subject,
    baseline_shannon,
    d5_displacement_auc = recovery_auc_bray_actual_day,
    d5_recovery_similarity = 1 - recovery_auc_bray_actual_day,
    z_baseline_shannon = z_score(baseline_shannon),
    z_recovery_score = z_score(-recovery_auc_bray_actual_day)
  )
primary_fit <- lm(z_recovery_score ~ z_baseline_shannon, data = primary_data)
primary_robust <- hc3_tidy(primary_fit)

# Participant-resampling and permutation sensitivities are included because
# the residual screen is non-Gaussian. They supplement, rather than replace,
# the prespecified HC3 slope and confidence interval.
set.seed(SEED + 200L)
primary_bootstrap_beta <- replicate(PRIMARY_BOOTSTRAP_R, {
  index <- sample.int(nrow(primary_data), nrow(primary_data), replace = TRUE)
  boot_data <- primary_data[index, ] |>
    transmute(
      z_baseline_shannon = z_score(baseline_shannon),
      z_recovery_score = z_score(-d5_displacement_auc)
    )
  unname(coef(lm(z_recovery_score ~ z_baseline_shannon, data = boot_data))["z_baseline_shannon"])
})
set.seed(SEED + 201L)
observed_beta <- unname(coef(primary_fit)["z_baseline_shannon"])
permuted_beta <- replicate(PRIMARY_PERMUTATIONS, {
  permuted_outcome <- sample(primary_data$z_recovery_score, replace = FALSE)
  unname(coef(lm(permuted_outcome ~ primary_data$z_baseline_shannon))[2])
})
primary_resampling <- tibble(
  analysis = c("Participant bootstrap", "Outcome-label permutation"),
  repetitions = c(PRIMARY_BOOTSTRAP_R, PRIMARY_PERMUTATIONS),
  observed_standardized_beta = observed_beta,
  estimate_or_null_median = c(median(primary_bootstrap_beta), median(permuted_beta)),
  interval_2.5_percent = c(quantile(primary_bootstrap_beta, 0.025, type = 6),
                           quantile(permuted_beta, 0.025, type = 6)),
  interval_97.5_percent = c(quantile(primary_bootstrap_beta, 0.975, type = 6),
                            quantile(permuted_beta, 0.975, type = 6)),
  two_sided_p_value = c(NA_real_,
                        (1 + sum(abs(permuted_beta) >= abs(observed_beta))) /
                          (PRIMARY_PERMUTATIONS + 1)),
  interpretation = c(
    "Percentile CI from participant-level resampling; all transformations refit within each draw.",
    "Two-sided randomization sensitivity under exchangeability of participant outcomes."
  )
)
primary_model <- primary_model |>
  mutate(
    participant_bootstrap_repetitions = PRIMARY_BOOTSTRAP_R,
    participant_bootstrap_CI95_low = primary_resampling$interval_2.5_percent[1],
    participant_bootstrap_CI95_high = primary_resampling$interval_97.5_percent[1],
    permutation_repetitions = PRIMARY_PERMUTATIONS,
    permutation_p_value = primary_resampling$two_sided_p_value[2]
  )
primary_residual <- residuals(primary_fit)
primary_standardized_residual <- rstandard(primary_fit)
primary_cook <- cooks.distance(primary_fit)
primary_leverage <- hatvalues(primary_fit)
shapiro_result <- shapiro.test(primary_residual)
bp_auxiliary <- lm(I(primary_residual^2) ~ fitted(primary_fit))
bp_statistic <- nrow(primary_data) * summary(bp_auxiliary)$r.squared
bp_p <- pchisq(bp_statistic, df = 1, lower.tail = FALSE)
quadratic_fit <- lm(z_recovery_score ~ z_baseline_shannon + I(z_baseline_shannon^2), data = primary_data)
quadratic_term <- hc3_tidy(quadratic_fit) |> filter(term == "I(z_baseline_shannon^2)")

diagnostics <- tibble(
  diagnostic = c(
    "Residual Shapiro-Wilk W", "Residual Shapiro-Wilk P",
    "Koenker-style BP nR2", "Koenker-style BP P",
    "Maximum Cook's distance", "Cook threshold 4/n", "Subjects above Cook threshold",
    "Maximum leverage", "Leverage threshold 2p/n", "Subjects above leverage threshold",
    "Absolute standardized residual > 3", "Quadratic HC3 coefficient",
    "Quadratic HC3 95% CI low", "Quadratic HC3 95% CI high", "Quadratic HC3 P"
  ),
  value = c(
    unname(shapiro_result$statistic), shapiro_result$p.value,
    bp_statistic, bp_p,
    max(primary_cook), 4 / nrow(primary_data), sum(primary_cook > 4 / nrow(primary_data)),
    max(primary_leverage), 2 * length(coef(primary_fit)) / nrow(primary_data),
    sum(primary_leverage > 2 * length(coef(primary_fit)) / nrow(primary_data)),
    sum(abs(primary_standardized_residual) > 3), quadratic_term$estimate,
    quadratic_term$CI95_low, quadratic_term$CI95_high, quadratic_term$p_value
  ),
  note = c(
    "Distributional screen; inference uses HC3 covariance.",
    "A P value is not proof of normality.",
    "Auxiliary regression of squared residuals on fitted values.",
    "A P value is not proof of homoscedasticity; HC3 remains primary.",
    primary_data$subject[which.max(primary_cook)],
    "Screening threshold, not an exclusion rule.",
    paste(primary_data$subject[primary_cook > 4 / nrow(primary_data)], collapse = ";"),
    primary_data$subject[which.max(primary_leverage)],
    "Screening threshold, not an exclusion rule.",
    paste(primary_data$subject[primary_leverage > 2 * length(coef(primary_fit)) / nrow(primary_data)], collapse = ";"),
    "No participant is removed on this basis.",
    "Prespecified quadratic sensitivity for nonlinearity.",
    "Prespecified quadratic sensitivity for nonlinearity.",
    "Prespecified quadratic sensitivity for nonlinearity.",
    "P for quadratic term with HC3 covariance."
  )
)

leave_one_out <- lapply(seq_len(nrow(primary_data)), function(index) {
  analysis <- primary_data[-index, ] |>
    transmute(
      subject,
      z_baseline_shannon = z_score(baseline_shannon),
      z_recovery_score = z_score(-d5_displacement_auc)
    )
  fit <- lm(z_recovery_score ~ z_baseline_shannon, data = analysis)
  term <- hc3_tidy(fit) |> filter(term == "z_baseline_shannon")
  tibble(
    omitted_subject = primary_data$subject[index],
    n_subjects = nrow(analysis),
    standardized_beta = term$estimate,
    HC3_SE = term$HC3_SE,
    CI95_low = term$CI95_low,
    CI95_high = term$CI95_high,
    p_value = term$p_value
  )
}) |> bind_rows()
leave_one_out_summary <- leave_one_out |>
  summarise(
    repetitions = n(),
    beta_min = min(standardized_beta),
    beta_median = median(standardized_beta),
    beta_max = max(standardized_beta),
    proportion_positive = mean(standardized_beta > 0),
    minimum_p = min(p_value),
    maximum_p = max(p_value)
  )

# Baseline Shannon reliability across A1--A4.
baseline_matrix_frame <- sample_primary |>
  filter(canon %in% paste0("A", 1:4)) |>
  select(subject, canon, shannon) |>
  pivot_wider(names_from = canon, values_from = shannon) |>
  arrange(subject)
if (nrow(baseline_matrix_frame) != 60L || any(is.na(baseline_matrix_frame))) stop("Incomplete A1--A4 matrix")
baseline_matrix <- as.matrix(baseline_matrix_frame[, paste0("A", 1:4)])
baseline_reliability <- bootstrap_icc(baseline_matrix, ICC_BOOTSTRAP_R, SEED + 100L) |>
  mutate(
    label = recode(reliability,
                   ICC_A_1 = "Single baseline measurement: ICC(A,1)",
                   ICC_A_k = "Mean of A1-A4: ICC(A,4)"),
    n_subjects = nrow(baseline_matrix),
    n_baseline_measurements = ncol(baseline_matrix)
  )

# Cutoff sensitivity.  D5-defined groups are the reviewer-facing sensitivity;
# the historical endpoint-Shannon definition is retained only as a clearly
# selection-coupled reference audit.
cutoff_specs <- tribble(
  ~scheme, ~fraction_each_tail, ~n_per_tail, ~scheme_label,
  "quintile", 0.20, 12L, "Quintile tails (20%; n=12 each)",
  "quartile", 0.25, 15L, "Quartile tails (25%; n=15 each)",
  "tertile", 1 / 3, 20L, "Tertile tails (33.3%; n=20 each)"
)
selection_specs <- tribble(
  ~selection_definition, ~selection_column, ~high_recovery_at_low_value, ~selection_label,
  "D5_recovery", "recovery_auc_bray_actual_day", TRUE,
  "D5 Bray-Curtis recovery AUC (primary recovery definition)",
  "historical_endpoint_shannon", "grouping_endpoint_shannon", FALSE,
  "Historical absolute endpoint Shannon (selection-coupled reference only)"
)

cutoff_results_list <- list()
cutoff_assignments_list <- list()
counter <- 1L
for (selection_index in seq_len(nrow(selection_specs))) {
  selection <- selection_specs[selection_index, ]
  ordering <- order(subject[[selection$selection_column]], subject$subject, na.last = NA)
  if (!selection$high_recovery_at_low_value) ordering <- rev(ordering)
  sorted_values <- subject[[selection$selection_column]][ordering]
  sorted_subjects <- subject$subject[ordering]
  for (cutoff_index in seq_len(nrow(cutoff_specs))) {
    cutoff <- cutoff_specs[cutoff_index, ]
    n_tail <- cutoff$n_per_tail
    high_ids <- sorted_subjects[seq_len(n_tail)]
    low_ids <- sorted_subjects[(length(sorted_subjects) - n_tail + 1):length(sorted_subjects)]
    high_values <- subject$baseline_shannon[match(high_ids, subject$subject)]
    low_values <- subject$baseline_shannon[match(low_ids, subject$subject)]
    welch <- t.test(high_values, low_values, var.equal = FALSE)
    bootstrap <- bootstrap_ranked_tail_contrast(
      subject, selection$selection_column, "baseline_shannon", n_tail,
      selection$high_recovery_at_low_value, BOOTSTRAP_R,
      SEED + selection_index * 1000L + n_tail
    )
    mean_ci <- bootstrap |> filter(statistic == "mean_difference")
    g_ci <- bootstrap |> filter(statistic == "hedges_g")
    high_boundary_tie <- isTRUE(all.equal(sorted_values[n_tail], sorted_values[n_tail + 1], tolerance = 1e-12))
    low_boundary <- length(sorted_values) - n_tail
    low_boundary_tie <- isTRUE(all.equal(sorted_values[low_boundary], sorted_values[low_boundary + 1], tolerance = 1e-12))
    cutoff_results_list[[counter]] <- tibble(
      selection_definition = selection$selection_definition,
      selection_label = selection$selection_label,
      scheme = cutoff$scheme,
      scheme_label = cutoff$scheme_label,
      fraction_each_tail = cutoff$fraction_each_tail,
      n_per_tail = n_tail,
      n_middle = nrow(subject) - 2 * n_tail,
      high_group_label = if_else(selection$selection_definition == "D5_recovery", "Higher recovery", "Higher endpoint Shannon"),
      low_group_label = if_else(selection$selection_definition == "D5_recovery", "Lower recovery", "Lower endpoint Shannon"),
      baseline_shannon_high_mean = mean(high_values),
      baseline_shannon_high_SD = sd(high_values),
      baseline_shannon_low_mean = mean(low_values),
      baseline_shannon_low_SD = sd(low_values),
      mean_difference_high_minus_low = mean(high_values) - mean(low_values),
      mean_difference_boot_CI95_low = mean_ci$CI95_low,
      mean_difference_boot_CI95_high = mean_ci$CI95_high,
      hedges_g = hedges_g(high_values, low_values),
      hedges_g_boot_CI95_low = g_ci$CI95_low,
      hedges_g_boot_CI95_high = g_ci$CI95_high,
      welch_t = unname(welch$statistic),
      welch_df = unname(welch$parameter),
      welch_p_value = welch$p.value,
      high_boundary_tie = high_boundary_tie,
      low_boundary_tie = low_boundary_tie,
      bootstrap_repetitions = BOOTSTRAP_R,
      bootstrap_design = "Full-cohort pairs bootstrap; selection ranking and tail membership recomputed within every resample.",
      inference_role = "Secondary cutoff sensitivity; effect direction/magnitude and CI are primary, not threshold-specific significance."
    )
    all_subject_ids <- subject$subject
    all_selection_values <- subject[[selection$selection_column]]
    all_baseline_shannon <- subject$baseline_shannon
    all_d5_values <- subject$recovery_auc_bray_actual_day
    cutoff_assignments_list[[counter]] <- tibble(
      subject = all_subject_ids,
      selection_definition = selection$selection_definition,
      scheme = cutoff$scheme,
      n_per_tail = n_tail,
      group = case_when(
        all_subject_ids %in% high_ids ~ "High",
        all_subject_ids %in% low_ids ~ "Low",
        TRUE ~ "Middle"
      ),
      selection_value = all_selection_values,
      baseline_shannon = all_baseline_shannon,
      D5_recovery_auc_bray = all_d5_values
    )
    counter <- counter + 1L
  }
}
cutoff_results <- bind_rows(cutoff_results_list) |>
  group_by(selection_definition) |>
  mutate(welch_p_Holm_across_cutoffs = p.adjust(welch_p_value, method = "holm")) |>
  ungroup()
cutoff_assignments <- bind_rows(cutoff_assignments_list)
if (any(cutoff_results$high_boundary_tie | cutoff_results$low_boundary_tie)) {
  stop("At least one cutoff boundary contains a tie; a prespecified tie rule is required")
}
historical_n15 <- cutoff_assignments |>
  filter(selection_definition == "historical_endpoint_shannon", scheme == "quartile") |>
  left_join(subject |> select(subject, original_endpoint_group), by = "subject") |>
  mutate(expected = case_when(
    original_endpoint_group == "HDR" ~ "High",
    original_endpoint_group == "LDR" ~ "Low",
    TRUE ~ "Middle"
  ))
if (any(historical_n15$group != historical_n15$expected)) {
  stop("Historical n=15 reconstruction does not match original HDR/LDR membership")
}

# Sequencing-depth sensitivity: Methods-faithful fixed-depth branch, an
# unrarefied branch, and an unrarefied model adjusted for average baseline and
# recovery-period read depths.
unrarefied_depth <- sample_unrarefied |>
  group_by(subject) |>
  summarise(
    baseline_mean_depth_recomputed = mean(read_depth[canon %in% paste0("A", 1:4)]),
    recovery_mean_depth = mean(read_depth[canon %in% c("B6", paste0("C", 1:6))]),
    .groups = "drop"
  )
depth_data <- subject_unrarefied |>
  select(subject, baseline_shannon, recovery_auc_bray_actual_day, baseline_mean_read_depth) |>
  left_join(unrarefied_depth, by = "subject") |>
  mutate(
    recovery_score = -recovery_auc_bray_actual_day,
    z_recovery_score = z_score(recovery_score),
    z_baseline_shannon = z_score(baseline_shannon),
    z_log_baseline_depth = z_score(log10(baseline_mean_read_depth)),
    z_log_recovery_depth = z_score(log10(recovery_mean_depth))
  )
if (max(abs(depth_data$baseline_mean_read_depth - depth_data$baseline_mean_depth_recomputed)) > 1e-6) {
  stop("Subject/sample baseline depth summaries do not agree")
}
depth_models <- list(
  rarefied_unadjusted = primary_fit,
  unrarefied_unadjusted = lm(z_recovery_score ~ z_baseline_shannon, data = depth_data),
  unrarefied_depth_adjusted = lm(
    z_recovery_score ~ z_baseline_shannon + z_log_baseline_depth + z_log_recovery_depth,
    data = depth_data
  )
)
depth_sensitivity <- lapply(names(depth_models), function(model_name) {
  fit <- depth_models[[model_name]]
  hc3_tidy(fit) |>
    mutate(
      model = model_name,
      n_subjects = nobs(fit),
      R2 = summary(fit)$r.squared,
      adjusted_R2 = summary(fit)$adj.r.squared,
      dataset_branch = if_else(grepl("^rarefied", model_name), "rarefied_primary", "unrarefied_sensitivity")
    ) |>
    select(dataset_branch, model, n_subjects, term, everything())
}) |> bind_rows()
depth_predictor_association <- suppressWarnings(cor.test(
  depth_data$baseline_shannon, log10(depth_data$baseline_mean_read_depth),
  method = "spearman", exact = FALSE
))
depth_audit <- tibble(
  n_unrarefied_samples = nrow(sample_unrarefied),
  n_subjects = n_distinct(sample_unrarefied$subject),
  read_depth_minimum = min(sample_unrarefied$read_depth),
  read_depth_median = median(sample_unrarefied$read_depth),
  read_depth_maximum = max(sample_unrarefied$read_depth),
  read_depth_CV = sd(sample_unrarefied$read_depth) / mean(sample_unrarefied$read_depth),
  baseline_Shannon_vs_log10_baseline_depth_Spearman_rho = unname(depth_predictor_association$estimate),
  baseline_Shannon_vs_log10_baseline_depth_p = depth_predictor_association$p.value,
  primary_fixed_depth = RAREFACTION_DEPTH
)

# Repeated exact rarefaction. Replicate 1 must reproduce the primary branch.
rarefaction_subject <- read_csv(rarefaction_path, show_col_types = FALSE)
expected_replicates <- sort(unique(rarefaction_subject$replicate))
if (!isTRUE(all.equal(as.numeric(expected_replicates), as.numeric(seq_len(max(expected_replicates)))))) {
  stop("Rarefaction replicate indices are incomplete")
}
if (any(table(rarefaction_subject$replicate) != 60L)) stop("Each rarefaction replicate must contain 60 subjects")
replicate_one <- rarefaction_subject |> filter(replicate == 1L) |> arrange(subject)
primary_sorted <- subject |> arrange(subject)
validation_delta <- tibble(
  variable = c("baseline_shannon", "recovery_auc_bray_actual_day"),
  maximum_absolute_difference = c(
    max(abs(replicate_one$baseline_shannon - primary_sorted$baseline_shannon)),
    max(abs(replicate_one$recovery_auc_bray_actual_day - primary_sorted$recovery_auc_bray_actual_day))
  )
)
if (any(validation_delta$maximum_absolute_difference > 1e-10)) {
  stop("Repeated rarefaction replicate 1 does not reproduce the primary branch")
}

rarefaction_model_results <- rarefaction_subject |>
  group_split(replicate, .keep = TRUE) |>
  lapply(function(data) {
    analysis <- data |>
      mutate(
        z_baseline_shannon = z_score(baseline_shannon),
        z_recovery_score = z_score(-recovery_auc_bray_actual_day)
      )
    fit <- lm(z_recovery_score ~ z_baseline_shannon, data = analysis)
    term <- hc3_tidy(fit) |> filter(term == "z_baseline_shannon")
    tibble(
      replicate = unique(data$replicate),
      seed = unique(data$seed),
      n_subjects = nrow(data),
      standardized_beta = term$estimate,
      HC3_SE = term$HC3_SE,
      statistic = term$statistic,
      df_residual = term$df_residual,
      p_value = term$p_value,
      CI95_low = term$CI95_low,
      CI95_high = term$CI95_high,
      R2 = summary(fit)$r.squared,
      adjusted_R2 = summary(fit)$adj.r.squared
    )
  }) |>
  bind_rows() |>
  arrange(replicate)
rarefaction_summary <- rarefaction_model_results |>
  summarise(
    rarefaction_replicates = n(),
    seed_first = min(seed),
    seed_last = max(seed),
    beta_median = median(standardized_beta),
    beta_empirical_CI95_low = quantile(standardized_beta, 0.025, type = 6),
    beta_empirical_CI95_high = quantile(standardized_beta, 0.975, type = 6),
    beta_min = min(standardized_beta),
    beta_max = max(standardized_beta),
    proportion_positive = mean(standardized_beta > 0),
    proportion_nominal_p_below_0_05 = mean(p_value < 0.05),
    replicate_1_beta = standardized_beta[replicate == 1L],
    replicate_1_p = p_value[replicate == 1L]
  )

# Participant-level analysis table with all transformations and the D5 quartile
# display group.  The middle 30 remain in the full-cohort model.
d5_quartile_assignment <- cutoff_assignments |>
  filter(selection_definition == "D5_recovery", scheme == "quartile") |>
  select(subject, D5_quartile_display_group = group)
subject_analysis_table <- subject |>
  select(
    subject, baseline_shannon, baseline_shannon_sd, baseline_mean_read_depth,
    endpoint_shannon, grouping_endpoint_shannon, original_endpoint_group,
    all_of(metric_definitions$outcome_column)
  ) |>
  left_join(d5_quartile_assignment, by = "subject") |>
  mutate(
    D5_recovery_similarity = 1 - recovery_auc_bray_actual_day,
    z_baseline_shannon = z_score(baseline_shannon),
    z_D5_recovery_score = z_score(-recovery_auc_bray_actual_day),
    analysis_inclusion = "Included in full-cohort primary model"
  )

descriptive_statistics <- bind_rows(
  tibble(variable = "baseline_shannon", values = list(subject$baseline_shannon)),
  lapply(seq_len(nrow(metric_definitions)), function(index) {
    tibble(variable = metric_definitions$outcome_column[index], values = list(subject[[metric_definitions$outcome_column[index]]]))
  }) |> bind_rows()
) |>
  rowwise() |>
  transmute(
    variable,
    n = sum(is.finite(values)),
    mean = mean(values, na.rm = TRUE),
    SD = sd(values, na.rm = TRUE),
    median = median(values, na.rm = TRUE),
    Q1 = quantile(values, 0.25, na.rm = TRUE),
    Q3 = quantile(values, 0.75, na.rm = TRUE),
    minimum = min(values, na.rm = TRUE),
    maximum = max(values, na.rm = TRUE),
    missing = sum(!is.finite(values))
  ) |>
  ungroup()

# Write all numerical outputs before plotting.
write_csv(subject_analysis_table, file.path(table_dir, "R1_2_subject_analysis_dataset.csv"))
write_csv(descriptive_statistics, file.path(table_dir, "R1_2_descriptive_statistics.csv"))
write_csv(primary_model, file.path(table_dir, "R1_2_primary_full_cohort_model.csv"))
write_csv(primary_resampling, file.path(table_dir, "R1_2_primary_resampling_sensitivity.csv"))
write_csv(recovery_models, file.path(table_dir, "R1_2_D1_D6_continuous_models.csv"))
write_csv(reference_analysis, file.path(table_dir, "R1_2_reference_endpoint_full_cohort.csv"))
write_csv(spearman_sensitivity, file.path(table_dir, "R1_2_spearman_sensitivity.csv"))
write_csv(diagnostics, file.path(table_dir, "R1_2_primary_model_diagnostics.csv"))
write_csv(leave_one_out, file.path(table_dir, "R1_2_leave_one_out.csv"))
write_csv(leave_one_out_summary, file.path(table_dir, "R1_2_leave_one_out_summary.csv"))
write_csv(baseline_reliability, file.path(table_dir, "R1_2_baseline_Shannon_reliability.csv"))
write_csv(cutoff_results, file.path(table_dir, "R1_2_cutoff_sensitivity.csv"))
write_csv(cutoff_assignments, file.path(table_dir, "R1_2_cutoff_group_assignments.csv"))
write_csv(depth_sensitivity, file.path(table_dir, "R1_2_depth_sensitivity_models.csv"))
write_csv(depth_audit, file.path(table_dir, "R1_2_sequencing_depth_audit.csv"))
write_csv(rarefaction_model_results, file.path(table_dir, "R1_2_repeated_rarefaction_models.csv"))
write_csv(rarefaction_summary, file.path(table_dir, "R1_2_repeated_rarefaction_summary.csv"))
write_csv(validation_delta, file.path(table_dir, "R1_2_repeated_rarefaction_validation.csv"))

# Publication-style general response figure. Publisher-specific settings are
# intentionally not claimed because no final submission specification was
# supplied; PDF remains vector and TIFF is 600 dpi with an opaque white ground.
palette <- c(
  blue = "#0072B2", vermillion = "#D55E00", green = "#009E73",
  purple = "#7A5195", grey = "#6B7280", light_blue = "#D9EDF7",
  light_orange = "#F6D7C3", black = "#1F2937"
)
theme_r12 <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.title = element_text(size = 8, colour = palette["black"]),
    axis.text = element_text(size = 7, colour = palette["black"]),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0,
                              margin = margin(b = 2, l = 11)),
    plot.subtitle = element_text(size = 6.7, colour = palette["grey"], hjust = 0,
                                 margin = margin(b = 4, l = 11)),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.8),
    strip.text = element_text(size = 7.5, face = "bold"),
    plot.margin = margin(5, 6, 5, 5),
    panel.grid.major.x = element_line(colour = "#ECECEC", linewidth = 0.25),
    panel.grid.minor = element_blank()
  )

primary_subtitle <- sprintf(
  paste0(
    "All n = 60 retained; OLS fit with pointwise HC3 95%% CI\n",
    "Standardized beta = %.2f [%.2f, %.2f]; HC3 %s; R2 = %.02f"
  ),
  primary_model$standardized_beta,
  primary_model$beta_CI95_low,
  primary_model$beta_CI95_high,
  fmt_p(primary_model$beta_p_value),
  primary_model$R2
)
primary_display_fit <- lm(d5_recovery_similarity ~ baseline_shannon, data = primary_data)
primary_display_vcov <- hc3_vcov(primary_display_fit)
primary_band <- tibble(
  baseline_shannon = seq(
    min(primary_data$baseline_shannon),
    max(primary_data$baseline_shannon),
    length.out = 200L
  )
)
primary_band_matrix <- model.matrix(~ baseline_shannon, data = primary_band)
primary_band <- primary_band |>
  mutate(
    fitted = as.numeric(primary_band_matrix %*% coef(primary_display_fit)),
    robust_SE = sqrt(rowSums((primary_band_matrix %*% primary_display_vcov) * primary_band_matrix)),
    CI95_low = fitted - qt(0.975, df.residual(primary_display_fit)) * robust_SE,
    CI95_high = fitted + qt(0.975, df.residual(primary_display_fit)) * robust_SE
  )
panel_a <- ggplot(primary_data, aes(x = baseline_shannon, y = d5_recovery_similarity)) +
  geom_ribbon(
    data = primary_band,
    aes(x = baseline_shannon, ymin = CI95_low, ymax = CI95_high),
    inherit.aes = FALSE, fill = palette["light_orange"], alpha = 0.55
  ) +
  geom_line(
    data = primary_band, aes(x = baseline_shannon, y = fitted),
    inherit.aes = FALSE, colour = palette["vermillion"], linewidth = 0.8
  ) +
  geom_point(shape = 21, size = 2.1, stroke = 0.45,
             fill = palette["blue"], colour = "white", alpha = 0.88) +
  labs(
    title = "Requested full-cohort D5 sensitivity (n = 60)",
    subtitle = primary_subtitle,
    x = "Baseline Shannon diversity (mean A1-A4)",
    y = "D5 recovery similarity (1 - Bray-Curtis AUC)"
  ) +
  theme_r12

forest_data <- recovery_models |>
  mutate(
    plot_label = recode(metric_id,
      D1 = "D1  Endpoint |Shannon deviation|",
      D2 = "D2  Endpoint Bray-Curtis",
      D3 = "D3  Endpoint Jaccard",
      D4 = "D4  Shannon AUC",
      D5 = "D5  Bray-Curtis AUC",
      D6 = "D6  Jaccard AUC"
    ),
    plot_label = factor(plot_label, levels = rev(plot_label)),
    role = if_else(primary, "D5 continuous sensitivity", "Other recovery metrics")
  )
panel_b <- ggplot(forest_data, aes(x = standardized_beta, y = plot_label)) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.45, linetype = "dashed") +
  geom_errorbar(aes(xmin = beta_CI95_low, xmax = beta_CI95_high, colour = role),
                orientation = "y", width = 0.16, linewidth = 0.65) +
  geom_point(aes(fill = role, shape = role), colour = "white", stroke = 0.45, size = 2.6) +
  scale_colour_manual(values = c("D5 continuous sensitivity" = unname(palette["vermillion"]),
                                 "Other recovery metrics" = unname(palette["blue"]))) +
  scale_fill_manual(values = c("D5 continuous sensitivity" = unname(palette["vermillion"]),
                               "Other recovery metrics" = unname(palette["blue"]))) +
  scale_shape_manual(values = c("D5 continuous sensitivity" = 23, "Other recovery metrics" = 21)) +
  guides(colour = "none", fill = "none", shape = "none") +
  labs(
    title = "Recovery-metric sensitivity",
    subtitle = "D5 highlighted in orange; other definitions in blue\nHC3 95% CI; multiplicity-adjusted P values in table",
    x = "Standardized beta\n(positive = better recovery)", y = NULL
  ) +
  theme_r12

cutoff_plot_data <- cutoff_results |>
  filter(selection_definition == "D5_recovery") |>
  mutate(
    short_label = case_when(
      n_per_tail == 12L ~ "Quintile (n = 12/tail)",
      n_per_tail == 15L ~ "Quartile (n = 15/tail)",
      n_per_tail == 20L ~ "Tertile (n = 20/tail)"
    ),
    scheme_label = factor(short_label,
                          levels = rev(c("Quintile (n = 12/tail)",
                                         "Quartile (n = 15/tail)",
                                         "Tertile (n = 20/tail)"))),
    reference_size = if_else(n_per_tail == 15L, "Same n=15 size", "Alternative cutoff")
  )
panel_c <- ggplot(cutoff_plot_data, aes(x = hedges_g, y = scheme_label)) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.45, linetype = "dashed") +
  geom_errorbar(aes(xmin = hedges_g_boot_CI95_low, xmax = hedges_g_boot_CI95_high,
                    colour = reference_size), orientation = "y", width = 0.16, linewidth = 0.7) +
  geom_point(aes(fill = reference_size, shape = reference_size),
             colour = "white", stroke = 0.45, size = 2.8) +
  scale_colour_manual(values = c("Same n=15 size" = unname(palette["vermillion"]),
                                 "Alternative cutoff" = unname(palette["blue"]))) +
  scale_fill_manual(values = c("Same n=15 size" = unname(palette["vermillion"]),
                               "Alternative cutoff" = unname(palette["blue"]))) +
  scale_shape_manual(values = c("Same n=15 size" = 23, "Alternative cutoff" = 21)) +
  guides(colour = "none", fill = "none", shape = "none") +
  labs(
    title = "Extreme-cutoff sensitivity",
    subtitle = paste0(
      "D5 tails re-ranked in 5,000 full-cohort bootstraps\n",
      "n=15 highlighted (same size, not historical members)"
    ),
    x = "Hedges' g for baseline Shannon\n(higher- minus lower-recovery tail)", y = NULL
  ) +
  theme_r12

depth_baseline_terms <- depth_sensitivity |>
  filter(term == "z_baseline_shannon")
depth_plot_data <- bind_rows(
  tibble(
    display = "Repeated rarefaction\n(seed variability only)",
    estimate = rarefaction_summary$beta_median,
    low = rarefaction_summary$beta_empirical_CI95_low,
    high = rarefaction_summary$beta_empirical_CI95_high,
    interval_type = "Middle 95% of seed estimates",
    role = "Repeated rarefaction"
  ),
  depth_baseline_terms |>
    transmute(
      display = recode(model,
        rarefied_unadjusted = "Fixed-depth primary model",
        unrarefied_unadjusted = "Unrarefied model",
        unrarefied_depth_adjusted = "Unrarefied + depth covariates"
      ),
      estimate,
      low = CI95_low,
      high = CI95_high,
      interval_type = "HC3 95% confidence interval",
      role = recode(model,
        rarefied_unadjusted = "Primary seed",
        unrarefied_unadjusted = "Unrarefied",
        unrarefied_depth_adjusted = "Depth-adjusted"
      )
    )
) |>
  mutate(display = factor(display, levels = rev(c(
    "Repeated rarefaction\n(seed variability only)", "Fixed-depth primary model",
    "Unrarefied model", "Unrarefied + depth covariates"
  ))))
panel_d <- ggplot(depth_plot_data, aes(x = estimate, y = display)) +
  geom_vline(xintercept = 0, colour = "#7A7A7A", linewidth = 0.45, linetype = "dashed") +
  geom_errorbar(data = depth_plot_data |> filter(role != "Repeated rarefaction"),
                aes(xmin = low, xmax = high, colour = role),
                orientation = "y", width = 0.16, linewidth = 0.7) +
  geom_segment(data = depth_plot_data |> filter(role == "Repeated rarefaction"),
               aes(x = low, xend = high, y = display, yend = display),
               inherit.aes = FALSE, colour = palette["purple"],
               linewidth = 2.2, alpha = 0.45, lineend = "round") +
  geom_point(aes(fill = role, shape = role), colour = "white", stroke = 0.45, size = 2.8) +
  scale_colour_manual(values = c(
    "Repeated rarefaction" = unname(palette["purple"]),
    "Primary seed" = unname(palette["vermillion"]),
    "Unrarefied" = unname(palette["blue"]),
    "Depth-adjusted" = unname(palette["green"])
  )) +
  scale_fill_manual(values = c(
    "Repeated rarefaction" = unname(palette["purple"]),
    "Primary seed" = unname(palette["vermillion"]),
    "Unrarefied" = unname(palette["blue"]),
    "Depth-adjusted" = unname(palette["green"])
  )) +
  scale_shape_manual(values = c(
    "Repeated rarefaction" = 22, "Primary seed" = 23,
    "Unrarefied" = 21, "Depth-adjusted" = 24
  )) +
  guides(colour = "none", fill = "none", shape = "none") +
  labs(
    title = "Depth and rarefaction sensitivity",
    subtitle = "Model rows: HC3 95% CI\nPurple row: 95% seed range (Monte Carlo, not CI)",
    x = "Standardized beta\n(positive = better recovery)", y = NULL
  ) +
  theme_r12

response_figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(widths = c(1, 1.08), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10), plot.tag.position = c(0, 1))

pdf_path <- file.path(figure_dir, "R1_2_response_figure.pdf")
tif_path <- file.path(figure_dir, "R1_2_response_figure.tif")
png_path <- file.path(figure_dir, "R1_2_response_figure_preview.png")
save_embedded_pdf(pdf_path, response_figure, width_mm = 210, height_mm = 180)
ggsave(tif_path, response_figure, device = "tiff", width = 210, height = 180,
       units = "mm", dpi = 600, compression = "lzw", bg = "white")
ggsave(png_path, response_figure, device = "png", width = 210, height = 180,
       units = "mm", dpi = 300, bg = "white")

# Diagnostic figure for model auditing, not intended as the main response figure.
diagnostic_frame <- primary_data |>
  mutate(
    fitted = fitted(primary_fit),
    residual = primary_residual,
    standardized_residual = primary_standardized_residual,
    sqrt_abs_standardized_residual = sqrt(abs(primary_standardized_residual)),
    cook_distance = primary_cook,
    leverage = primary_leverage,
    observation = row_number(),
    theoretical_quantile = qnorm((rank(primary_standardized_residual, ties.method = "first") - 0.5) / n())
  )
top_cook_subjects <- diagnostic_frame |>
  slice_max(cook_distance, n = 3, with_ties = FALSE)
diag_a <- ggplot(diagnostic_frame, aes(fitted, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888", linewidth = 0.4) +
  geom_point(shape = 21, fill = palette["blue"], colour = "white", size = 1.8, stroke = 0.35) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = palette["vermillion"], linewidth = 0.65) +
  labs(title = "Residuals vs fitted", x = "Fitted standardized recovery", y = "Residual") + theme_r12
diag_b <- ggplot(diagnostic_frame, aes(sample = standardized_residual)) +
  stat_qq(shape = 21, fill = palette["blue"], colour = "white", size = 1.8, stroke = 0.35) +
  stat_qq_line(linetype = "dashed", colour = "#888888", linewidth = 0.4) +
  labs(title = "Normal Q-Q screen", x = "Theoretical quantile", y = "Ordered standardized residual") + theme_r12
diag_c <- ggplot(diagnostic_frame, aes(fitted, sqrt_abs_standardized_residual)) +
  geom_point(shape = 21, fill = palette["blue"], colour = "white", size = 1.8, stroke = 0.35) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = palette["vermillion"], linewidth = 0.65) +
  labs(title = "Scale-location", x = "Fitted standardized recovery", y = expression(sqrt("|standardized residual|"))) + theme_r12
diag_d <- ggplot(diagnostic_frame, aes(observation, cook_distance)) +
  geom_hline(yintercept = 4 / nrow(diagnostic_frame), linetype = "dashed",
             colour = palette["vermillion"], linewidth = 0.45) +
  annotate("text", x = Inf, y = 4 / nrow(diagnostic_frame),
           label = "4/n screening threshold", hjust = 1.03, vjust = -0.45,
           colour = palette["vermillion"], size = 2.05) +
  geom_segment(aes(xend = observation, y = 0, yend = cook_distance), colour = palette["blue"], linewidth = 0.45) +
  geom_point(colour = palette["blue"], size = 1.2) +
  geom_text(data = top_cook_subjects, aes(label = subject), size = 2.2,
            nudge_y = 0.004, check_overlap = TRUE) +
  labs(title = "Influence screen", x = "Participant index", y = "Cook's distance") + theme_r12
diagnostic_figure <- (diag_a | diag_b) / (diag_c | diag_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0, 1))
diagnostic_pdf <- file.path(figure_dir, "R1_2_primary_model_diagnostics.pdf")
diagnostic_tif <- file.path(figure_dir, "R1_2_primary_model_diagnostics.tif")
save_embedded_pdf(diagnostic_pdf, diagnostic_figure, width_mm = 180, height_mm = 145)
ggsave(diagnostic_tif, diagnostic_figure, device = "tiff", width = 180, height = 145,
       units = "mm", dpi = 600, compression = "lzw", bg = "white")

figure_source <- bind_rows(
  primary_data |> transmute(panel = "A_full_cohort", subject, x = baseline_shannon,
                            y = d5_recovery_similarity, label = NA_character_),
  primary_band |> transmute(
    panel = "A_HC3_fit_and_band", subject = NA_character_, x = baseline_shannon,
    y = fitted,
    label = sprintf("HC3_pointwise_CI95_low=%.8f;high=%.8f", CI95_low, CI95_high)
  ),
  forest_data |> transmute(panel = "B_metric_forest", subject = NA_character_,
                           x = standardized_beta, y = beta_CI95_low,
                           label = paste(metric_id, metric_label, "CI_high", beta_CI95_high)),
  cutoff_plot_data |> transmute(panel = "C_cutoff_forest", subject = NA_character_,
                                x = hedges_g, y = hedges_g_boot_CI95_low,
                                label = paste(scheme_label, "CI_high", hedges_g_boot_CI95_high)),
  depth_plot_data |> transmute(panel = "D_depth_rarefaction", subject = NA_character_,
                               x = estimate, y = low,
                               label = paste(as.character(display), interval_type, "high", high))
)
write_csv(figure_source, file.path(table_dir, "R1_2_response_figure_source_data.csv"))

# Dynamic result summary and figure legend.  This is analysis-ready language,
# not a claim that manuscript page/line changes have already been made.
primary_spearman <- spearman_sensitivity |> filter(primary)
d5_cutoffs <- cutoff_results |> filter(selection_definition == "D5_recovery")
icc_average <- baseline_reliability |> filter(reliability == "ICC_A_k")
d3_result <- recovery_models |> filter(metric_id == "D3")
d6_result <- recovery_models |> filter(metric_id == "D6")
endpoint_reference <- reference_analysis |> slice(1)
results_lines <- c(
  "# R1.2 analysis results / R1.2分析结果",
  "",
  sprintf("- Full cohort: all %d exposed participants were retained; no middle participants were discarded.", primary_model$n_subjects),
  sprintf("- Primary D5 model: standardized beta = %.3f (HC3 95%% CI %.3f to %.3f; %s; R2 = %.3f). Positive values mean that higher baseline Shannon predicts better Bray-Curtis recovery.",
          primary_model$standardized_beta, primary_model$beta_CI95_low,
          primary_model$beta_CI95_high, fmt_p(primary_model$beta_p_value), primary_model$R2),
  sprintf("- On the original D5 recovery-score scale, a 1-SD higher baseline Shannon corresponded to +%.4f recovery units (HC3 95%% CI %.5f to %.5f).",
          primary_model$raw_score_slope_per_1SD_baseline,
          primary_model$raw_score_CI95_low, primary_model$raw_score_CI95_high),
  "- The D5 confidence interval includes zero; the primary model therefore does not provide clear evidence that baseline Shannon predicts Bray-Curtis AUC recovery.",
  sprintf("- Spearman sensitivity for D5: rho = %.3f (%s).", primary_spearman$spearman_rho, fmt_p(primary_spearman$p_value)),
  sprintf("- Participant-resampling sensitivity: 10,000-bootstrap 95%% CI %.3f to %.3f; 9,999-permutation P = %.3f.",
          primary_model$participant_bootstrap_CI95_low,
          primary_model$participant_bootstrap_CI95_high,
          primary_model$permutation_p_value),
  sprintf("- Repeated rarefaction (%d seeds): median beta = %.4f, middle 95%% Monte Carlo seed range %.4f to %.4f; %.1f%% of estimates were positive. This seed range is not a participant-level confidence interval.",
          rarefaction_summary$rarefaction_replicates, rarefaction_summary$beta_median,
          rarefaction_summary$beta_empirical_CI95_low, rarefaction_summary$beta_empirical_CI95_high,
          100 * rarefaction_summary$proportion_positive),
  sprintf("- Mean A1-A4 Shannon reliability: ICC(A,4) = %.3f (participant-bootstrap 95%% CI %.3f to %.3f).",
          icc_average$estimate, icc_average$CI95_low, icc_average$CI95_high),
  sprintf("- Secondary recovery metrics were definition-dependent: D3 endpoint Jaccard beta = %.3f (Holm P %.3g) and D6 Jaccard AUC beta = %.3f (Holm P %.3g), whereas D1, D2, and D4 did not show adjusted support.",
          d3_result$standardized_beta, d3_result$p_Holm_secondary,
          d6_result$standardized_beta, d6_result$p_Holm_secondary),
  sprintf("- Absolute endpoint Shannon remained strongly associated with baseline Shannon in the full cohort (beta = %.3f, 95%% CI %.3f to %.3f), but this endpoint level is not interpreted as recovery.",
          endpoint_reference$standardized_beta, endpoint_reference$CI95_low,
          endpoint_reference$CI95_high),
  sprintf("- D5-defined cutoff effects (Hedges' g, higher-minus-lower recovery): n=12 %.3f [%.3f, %.3f]; n=15 %.3f [%.3f, %.3f]; n=20 %.3f [%.3f, %.3f].",
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 20],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 20],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 20]),
  "- Cutoff intervals use full-cohort pairs bootstraps with D5 ranks and tail membership recomputed in every draw; n=15 denotes the same sample size as the historical comparison, not the same participants.",
  "- Interpretation boundary: the continuous D5 result is primary; tail comparisons are secondary illustrations. Absolute endpoint-Shannon correlations are not described as recovery evidence.",
  "",
  "## 中文要点",
  "",
  sprintf("- 全队列主分析保留全部%d名抗生素暴露参与者，不排除中间个体。", primary_model$n_subjects),
  sprintf("- D5主模型：标准化beta = %.3f（HC3 95%% CI %.3f至%.3f；%s；R2 = %.3f）。正值表示较高基线Shannon预测较好的Bray-Curtis恢复。",
          primary_model$standardized_beta, primary_model$beta_CI95_low,
          primary_model$beta_CI95_high, fmt_p(primary_model$beta_p_value), primary_model$R2),
  sprintf("- 在D5恢复分数的原始尺度上，基线Shannon每升高1个标准差，恢复分数增加%.4f（HC3 95%% CI %.5f至%.5f）。",
          primary_model$raw_score_slope_per_1SD_baseline,
          primary_model$raw_score_CI95_low, primary_model$raw_score_CI95_high),
  "- D5置信区间跨越0，因此主要模型未提供明确证据表明基线Shannon能够预测Bray-Curtis AUC恢复。",
  sprintf("- D5的Spearman敏感性分析：rho = %.3f（%s）。", primary_spearman$spearman_rho, fmt_p(primary_spearman$p_value)),
  sprintf("- 参与者重抽样敏感性：10,000次bootstrap 95%% CI为%.3f至%.3f；9,999次置换检验P = %.3f。",
          primary_model$participant_bootstrap_CI95_low,
          primary_model$participant_bootstrap_CI95_high,
          primary_model$permutation_p_value),
  sprintf("- %d次精确重复稀释中，beta中位数为%.4f，种子估计值中95%%范围为%.4f至%.4f；该范围不是参与者层面置信区间。",
          rarefaction_summary$rarefaction_replicates, rarefaction_summary$beta_median,
          rarefaction_summary$beta_empirical_CI95_low, rarefaction_summary$beta_empirical_CI95_high),
  sprintf("- A1-A4平均Shannon的可靠性为ICC(A,4) = %.3f（参与者bootstrap 95%% CI %.3f至%.3f）。",
          icc_average$estimate, icc_average$CI95_low, icc_average$CI95_high),
  sprintf("- 次要恢复指标表现出定义依赖性：D3终点Jaccard beta = %.3f（Holm P %.3g），D6 Jaccard AUC beta = %.3f（Holm P %.3g）；D1、D2和D4经校正后无支持。",
          d3_result$standardized_beta, d3_result$p_Holm_secondary,
          d6_result$standardized_beta, d6_result$p_Holm_secondary),
  sprintf("- 全队列中绝对终点Shannon仍与基线Shannon高度相关（beta = %.3f，95%% CI %.3f至%.3f），但该终点水平不被解释为恢复。",
          endpoint_reference$standardized_beta, endpoint_reference$CI95_low,
          endpoint_reference$CI95_high),
  sprintf("- D5定义的截断效应（高恢复减低恢复的Hedges' g）：n=12为%.3f [%.3f, %.3f]；n=15为%.3f [%.3f, %.3f]；n=20为%.3f [%.3f, %.3f]。",
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 12],
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 15],
          d5_cutoffs$hedges_g[d5_cutoffs$n_per_tail == 20],
          d5_cutoffs$hedges_g_boot_CI95_low[d5_cutoffs$n_per_tail == 20],
          d5_cutoffs$hedges_g_boot_CI95_high[d5_cutoffs$n_per_tail == 20]),
  "- 截断区间采用全队列配对bootstrap，每次均重新计算D5排序和尾部成员；n=15仅表示样本量与历史比较相同，并非相同参与者。",
  "- 主要推断仅来自全队列D5连续模型；五分位、四分位和三分位两端分组仅用于次要图示。"
)
writeLines(results_lines, file.path(analysis_dir, "R1_2_results_summary_bilingual.md"), useBytes = TRUE)

legend_lines <- c(
  "# Response Figure R1.2 legend",
  "",
  "**Response Figure R1.2 | Supplementary analyses addressing regression to the mean and cutoff choice.**",
  "(A) All 60 exposed participants are retained in the reviewer-requested continuous D5 sensitivity model. (B) Standardized slopes for D1-D6, with D5 highlighted as the prespecified continuous sensitivity metric. (C) D5-defined tails with 12, 15 or 20 participants per tail. (D) Sequencing-depth and repeated-rarefaction sensitivity. The original endpoint-Shannon 15+15 comparison remains the main phenotypic analysis of Figs. 2-4.",
  "",
  "**Accessibility/alt text.** Four-panel sensitivity figure. Panel A shows a shallow full-cohort association. Panel B shows dependence on recovery definition. Panel C shows all three cutoff confidence intervals crossing zero. Panel D shows similar estimates across sequencing-depth treatments and rarefaction seeds.",
  "",
  "**Primary-model diagnostic figure.** Panels show residuals versus fitted values, the normal Q-Q screen, scale-location, and Cook's-distance influence screening for the full-cohort D5 model. The orange dashed line is the conventional 4/n screening threshold, not an exclusion rule. All 60 participants remain in the primary analysis.",
  "",
  "Note: the original manuscript raster and legend disagree about the panel letter for the baseline-endpoint correlation. The response should refer to the 'baseline-endpoint correlation panel' until the revised figure legend is corrected."
)
writeLines(legend_lines, file.path(analysis_dir, "R1_2_figure_legend.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(log_dir, "R1_2_R_sessionInfo.txt"))

output_files <- c(
  list.files(table_dir, full.names = TRUE),
  list.files(figure_dir, full.names = TRUE),
  file.path(analysis_dir, "R1_2_results_summary_bilingual.md"),
  file.path(analysis_dir, "R1_2_figure_legend.md"),
  file.path(analysis_dir, "R1_2_analysis_outline.md"),
  file.path(analysis_dir, "README.md"),
  script_path,
  file.path(analysis_dir, "run_repeated_rarefaction.py"),
  file.path(analysis_dir, "validate_R1_2_outputs.R"),
  file.path(analysis_dir, "run_R1_2_pipeline.sh"),
  file.path(log_dir, "R1_2_R_sessionInfo.txt"),
  file.path(log_dir, "R1_2_repeated_rarefaction_manifest.json")
)
manifest <- list(
  analysis = "Reviewer 1.2 full-cohort continuous analysis",
  created_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  prespecified_primary = list(
    predictor = "mean A1-A4 Shannon diversity",
    outcome = "D5 actual-day-normalized B6-C6 Bray-Curtis displacement AUC",
    presentation_transform = "1 - D5 (higher = better recovery)",
    sampling_unit = "participant",
    n = 60,
    model = "OLS slope with HC3 covariance and residual degrees-of-freedom t reference"
  ),
  sensitivity = list(
    secondary_metrics = c("D1", "D2", "D3", "D4", "D6"),
    secondary_multiplicity = "Holm across five secondary baseline-Shannon slopes",
    cutoffs = c(12, 15, 20),
    cutoff_bootstrap_repetitions = BOOTSTRAP_R,
    cutoff_bootstrap_design = "Full-cohort pairs bootstrap with selection ranking and tail membership recomputed within each draw",
    rarefaction_replicates = rarefaction_summary$rarefaction_replicates,
    rarefaction_depth = RAREFACTION_DEPTH,
    participant_bootstrap_repetitions = PRIMARY_BOOTSTRAP_R,
    outcome_label_permutations = PRIMARY_PERMUTATIONS,
    leave_one_out_repetitions = nrow(leave_one_out),
    random_seeds = list(
      global_R_seed = SEED,
      ICC_bootstrap = SEED + 100L,
      primary_participant_bootstrap = SEED + 200L,
      primary_outcome_permutation = SEED + 201L,
      D5_cutoff_bootstrap = setNames(
        as.list(SEED + 1000L + cutoff_specs$n_per_tail),
        paste0("n", cutoff_specs$n_per_tail)
      ),
      historical_cutoff_audit_bootstrap = setNames(
        as.list(SEED + 2000L + cutoff_specs$n_per_tail),
        paste0("n", cutoff_specs$n_per_tail)
      )
    )
  ),
  input_sha256 = setNames(lapply(required_inputs, function(path) digest(path, algo = "sha256", file = TRUE)),
                          vapply(required_inputs, function(path) sub(paste0(root_dir, "/"), "", path, fixed = TRUE), character(1))),
  outputs = lapply(output_files[file.exists(output_files)], function(path) {
    list(
      path = sub(paste0(root_dir, "/"), "", path, fixed = TRUE),
      bytes = unname(file.info(path)$size),
      sha256 = digest(path, algo = "sha256", file = TRUE)
    )
  }),
  interpretive_boundary = c(
    "The D5 full-cohort model is primary.",
    "Extreme tails are secondary illustrations, not independent primary evidence.",
    "Absolute endpoint Shannon is not treated as a recovery outcome.",
    "No outlier or influence screen is used to delete participants."
  )
)
write_json(manifest, file.path(log_dir, "R1_2_analysis_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat(sprintf(
  "R1.2 analysis complete: n=%d, D5 beta=%.4f, 95%% CI [%.4f, %.4f], p=%.6g, rarefaction replicates=%d\n",
  primary_model$n_subjects, primary_model$standardized_beta,
  primary_model$beta_CI95_low, primary_model$beta_CI95_high,
  primary_model$beta_p_value, rarefaction_summary$rarefaction_replicates
))
