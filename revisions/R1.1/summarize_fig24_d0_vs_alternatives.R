#!/usr/bin/env Rscript

# Audit D0 versus the six alternative R1.1 groupings using only finalized
# Figure 2/Figure 4 bottom tables.  This script does not refit a model or run a
# new statistical test.  In particular:
#   * Figure 2 D/E group means come from Fig4B_paired_tests_BH.csv and the
#     already-computed Welch P values come from
#     figure2_trajectory_tests_high_fidelity.csv;
#   * all paired P/q values, correlations and network counts are copied from
#     their corresponding finalized tables;
#   * Figure 2F agreement is a deterministic comparison of the signs and
#     displayed-significance flags already present in the LEfSe annotation
#     tables.

options(stringsAsFactors = FALSE, warn = 1)

time_order <- c(
  "A1", "A2", "A3", "A4",
  "B1", "B2", "B3", "B4", "B5", "B6",
  "C1", "C2", "C3", "C4", "C5", "C6"
)
metrics <- c("shannon", "richness", "simpson")
groups <- c("HDR", "LDR")
stages <- c("Baseline", "B6", "C1", "C2", "C3", "End")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Cannot resolve script path from Rscript --file argument", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
output_root <- file.path(
  project_root, "reviewer_response", "R1.1", "full_original_layout"
)

definition_slugs <- c(
  D0 = NA_character_,
  D1 = "D1_endpoint_abs_shannon",
  D2 = "D2_endpoint_bray",
  D3 = "D3_endpoint_jaccard",
  D4 = "D4_auc_abs_shannon",
  D5 = "D5_auc_bray",
  D6 = "D6_auc_jaccard"
)
definition_labels <- c(
  D0 = "Original endpoint-Shannon grouping",
  D1 = "Endpoint absolute Shannon displacement",
  D2 = "Endpoint Bray-Curtis displacement",
  D3 = "Endpoint Jaccard displacement",
  D4 = "AUC absolute Shannon displacement",
  D5 = "AUC Bray-Curtis displacement",
  D6 = "AUC Jaccard displacement"
)

fail <- function(...) stop(sprintf(...), call. = FALSE)

assert_true <- function(condition, ...) {
  if (length(condition) != 1L || is.na(condition) || !condition) fail(...)
}

assert_columns <- function(frame, required, label) {
  missing <- setdiff(required, names(frame))
  if (length(missing)) {
    fail("%s lacks required columns: %s", label, paste(missing, collapse = ", "))
  }
}

key_string <- function(frame, keys) {
  do.call(paste, c(frame[keys], sep = "\u001f"))
}

assert_unique_keys <- function(frame, keys, label) {
  values <- key_string(frame, keys)
  duplicated_values <- unique(values[duplicated(values)])
  if (length(duplicated_values)) {
    fail(
      "%s has duplicated key(s) for [%s]; first duplicate: %s",
      label, paste(keys, collapse = ", "), duplicated_values[[1L]]
    )
  }
}

assert_exact_grid <- function(frame, keys, expected, label) {
  assert_unique_keys(frame, keys, label)
  observed_keys <- sort(key_string(frame, keys))
  expected_keys <- sort(key_string(expected, keys))
  if (!identical(observed_keys, expected_keys)) {
    missing <- setdiff(expected_keys, observed_keys)
    extra <- setdiff(observed_keys, expected_keys)
    fail(
      paste0(
        "%s does not match its required key grid; observed=%d, expected=%d, ",
        "missing(first 5)=%s, extra(first 5)=%s"
      ),
      label, length(observed_keys), length(expected_keys),
      paste(utils::head(missing, 5L), collapse = " | "),
      paste(utils::head(extra, 5L), collapse = " | ")
    )
  }
}

assert_finite <- function(values, label, allow_na = FALSE) {
  if (!allow_na && anyNA(values)) fail("%s contains missing values", label)
  present <- values[!is.na(values)]
  if (any(!is.finite(present))) fail("%s contains non-finite values", label)
}

parse_logical_strict <- function(values, label) {
  normalized <- tolower(trimws(as.character(values)))
  out <- rep(NA, length(normalized))
  out[normalized %in% c("true", "t", "1")] <- TRUE
  out[normalized %in% c("false", "f", "0")] <- FALSE
  if (anyNA(out)) {
    fail(
      "%s contains non-Boolean values: %s",
      label, paste(unique(values[is.na(out)]), collapse = ", ")
    )
  }
  out
}

relative_path <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(project_root, .Platform$file.sep)
  if (!startsWith(normalized, prefix)) fail("Path is outside project root: %s", path)
  substring(normalized, nchar(prefix) + 1L)
}

read_bottom_table <- function(path, label) {
  if (!file.exists(path)) fail("Missing %s: %s", label, path)
  frame <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(frame)) fail("%s is empty: %s", label, path)
  attr(frame, "source_path") <- relative_path(path)
  frame
}

table_dir_for <- function(definition_id) {
  if (definition_id == "D0") {
    return(file.path(
      project_root, "reproduction", "outputs", "integrated_fig124_raw", "tables"
    ))
  }
  file.path(output_root, unname(definition_slugs[[definition_id]]), "tables")
}

load_bundle <- function(definition_id) {
  table_dir <- table_dir_for(definition_id)
  paths <- c(
    fig2_trajectory = file.path(table_dir, "figure2_trajectory_tests_high_fidelity.csv"),
    fig2_correlation = file.path(table_dir, "figure2_panelC_high_fidelity.csv"),
    fig4b = file.path(table_dir, "Fig4B_paired_tests_BH.csv"),
    fig4c = file.path(table_dir, "Fig4C_fastspar_summary.csv"),
    fig4d = file.path(table_dir, "Fig4D_dominant_genus_correlations.csv"),
    fig2f = file.path(table_dir, "Fig2F_LEfSe_annotations.csv")
  )
  tables <- lapply(names(paths), function(name) {
    read_bottom_table(paths[[name]], sprintf("%s %s", definition_id, name))
  })
  names(tables) <- names(paths)

  trajectory <- tables$fig2_trajectory
  assert_columns(
    trajectory, c("metric", "canon", "p_welch_raw", "symbol"),
    sprintf("%s Figure 2 trajectory", definition_id)
  )
  expected_trajectory <- expand.grid(
    metric = c("richness", "shannon"), canon = time_order,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  assert_exact_grid(
    trajectory, c("metric", "canon"), expected_trajectory,
    sprintf("%s Figure 2 trajectory", definition_id)
  )
  assert_finite(
    trajectory$p_welch_raw,
    sprintf("%s Figure 2 Welch P values", definition_id)
  )
  assert_true(
    all(trajectory$p_welch_raw >= 0 & trajectory$p_welch_raw <= 1),
    "%s Figure 2 Welch P values fall outside [0,1]", definition_id
  )

  corr <- tables$fig2_correlation
  assert_columns(
    corr, c("group", "rho", "p", "n"),
    sprintf("%s Figure 2C", definition_id)
  )
  assert_exact_grid(
    corr, "group", data.frame(group = groups),
    sprintf("%s Figure 2C", definition_id)
  )
  assert_true(all(corr$n == 15L), "%s Figure 2C does not have n=15 per group", definition_id)
  assert_finite(corr$rho, sprintf("%s Figure 2C rho", definition_id))
  assert_finite(corr$p, sprintf("%s Figure 2C P", definition_id))

  fig4b <- tables$fig4b
  assert_columns(
    fig4b,
    c(
      "metric", "group", "canon", "n_pairs", "baseline_A1_mean",
      "current_mean", "paired_mean_difference", "p_paired_vs_A1_raw",
      "q_BH_within_metric_group"
    ),
    sprintf("%s Figure 4B", definition_id)
  )
  expected_fig4b <- expand.grid(
    metric = metrics, group = groups, canon = time_order,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  assert_exact_grid(
    fig4b, c("metric", "group", "canon"), expected_fig4b,
    sprintf("%s Figure 4B", definition_id)
  )
  assert_true(all(fig4b$n_pairs == 15L), "%s Figure 4B does not have 15 pairs", definition_id)
  assert_finite(fig4b$baseline_A1_mean, sprintf("%s Figure 4B A1 means", definition_id))
  assert_finite(fig4b$current_mean, sprintf("%s Figure 4B current means", definition_id))
  assert_finite(fig4b$paired_mean_difference, sprintf("%s Figure 4B differences", definition_id))
  non_a1 <- fig4b$canon != "A1"
  assert_finite(
    fig4b$p_paired_vs_A1_raw[non_a1],
    sprintf("%s Figure 4B non-A1 paired P", definition_id)
  )
  assert_finite(
    fig4b$q_BH_within_metric_group[non_a1],
    sprintf("%s Figure 4B non-A1 paired q", definition_id)
  )
  assert_true(
    all(abs(
      fig4b$paired_mean_difference -
        (fig4b$current_mean - fig4b$baseline_A1_mean)
    ) < 1e-9),
    "%s Figure 4B difference is inconsistent with its tabulated means", definition_id
  )

  fig4c <- tables$fig4c
  assert_columns(
    fig4c,
    c(
      "group", "stage", "canons", "n_subjects", "n_samples",
      "raw_p_strong_edges", "raw_p_strong_nodes", "rho_threshold_strict_gt",
      "alpha_strict_lt", "bootstraps", "eligible_nodes_caption_rule",
      "eligible_edges_caption_rule"
    ),
    sprintf("%s Figure 4C", definition_id)
  )
  expected_fig4c <- expand.grid(
    group = groups, stage = stages,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  assert_exact_grid(
    fig4c, c("group", "stage"), expected_fig4c,
    sprintf("%s Figure 4C", definition_id)
  )
  assert_true(all(fig4c$n_subjects == 15L), "%s Figure 4C has non-15 subject count", definition_id)
  assert_true(all(fig4c$bootstraps == 500L), "%s Figure 4C is not the 500-bootstrap run", definition_id)
  assert_true(
    all(abs(fig4c$rho_threshold_strict_gt - 0.6) < 1e-12) &&
      all(abs(fig4c$alpha_strict_lt - 0.05) < 1e-12),
    "%s Figure 4C network qualification thresholds drifted", definition_id
  )
  assert_true(
    all(fig4c$eligible_nodes_caption_rule == fig4c$raw_p_strong_nodes) &&
      all(fig4c$eligible_edges_caption_rule == fig4c$raw_p_strong_edges),
    "%s Figure 4C eligible-count columns disagree", definition_id
  )
  if (all(c("trajectory_nodes", "trajectory_edges") %in% names(fig4c))) {
    assert_true(
      all(fig4c$trajectory_nodes == fig4c$eligible_nodes_caption_rule) &&
        all(fig4c$trajectory_edges == fig4c$eligible_edges_caption_rule),
      "%s Figure 4C trajectory counts are not full eligible-network counts",
      definition_id
    )
  }

  fig4d <- tables$fig4d
  assert_columns(
    fig4d,
    c(
      "group", "metric", "canon", "genus", "n", "rho", "p_raw",
      "q_BH_within_curve_sensitivity"
    ),
    sprintf("%s Figure 4D", definition_id)
  )
  expected_fig4d <- expand.grid(
    group = groups, metric = metrics, canon = time_order,
    genus = c("Bacteroides", "Faecalibacterium"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  assert_exact_grid(
    fig4d, c("group", "metric", "canon", "genus"), expected_fig4d,
    sprintf("%s Figure 4D", definition_id)
  )
  assert_true(all(fig4d$n == 15L), "%s Figure 4D does not have n=15", definition_id)
  assert_finite(fig4d$rho, sprintf("%s Figure 4D rho", definition_id))
  assert_finite(fig4d$p_raw, sprintf("%s Figure 4D raw P", definition_id))
  assert_finite(
    fig4d$q_BH_within_curve_sensitivity,
    sprintf("%s Figure 4D BH q", definition_id)
  )

  fig2f <- tables$fig2f
  assert_columns(
    fig2f,
    c(
      "canon", "species", "lda_score_LDR_over_HDR", "p_exact_mannwhitney",
      "q_BH_displayed_45", "displayed_significant"
    ),
    sprintf("%s Figure 2F", definition_id)
  )
  assert_true(nrow(fig2f) == 720L, "%s Figure 2F does not have 720 cells", definition_id)
  assert_true(
    length(unique(fig2f$species)) == 45L,
    "%s Figure 2F does not have 45 displayed species", definition_id
  )
  assert_exact_grid(
    fig2f, c("canon", "species"),
    expand.grid(
      canon = time_order, species = sort(unique(fig2f$species)),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ),
    sprintf("%s Figure 2F", definition_id)
  )
  assert_finite(
    fig2f$lda_score_LDR_over_HDR,
    sprintf("%s Figure 2F LDA scores", definition_id)
  )
  assert_true(
    all(fig2f$lda_score_LDR_over_HDR != 0),
    "%s Figure 2F contains zero-direction LDA cells", definition_id
  )
  assert_finite(fig2f$q_BH_displayed_45, sprintf("%s Figure 2F q", definition_id))
  fig2f$displayed_significant <- parse_logical_strict(
    fig2f$displayed_significant,
    sprintf("%s Figure 2F displayed_significant", definition_id)
  )
  assert_true(
    identical(fig2f$displayed_significant, fig2f$q_BH_displayed_45 < 0.05),
    "%s Figure 2F significance flags disagree with q<0.05", definition_id
  )
  tables$fig2f <- fig2f
  tables$paths <- vapply(paths, relative_path, character(1L))
  tables
}

bundles <- lapply(names(definition_slugs), load_bundle)
names(bundles) <- names(definition_slugs)

# The displayed Figure 2F universe must be exactly the same across all runs.
d0_fig2f_keys <- sort(key_string(bundles$D0$fig2f, c("canon", "species")))
for (definition_id in names(bundles)[-1L]) {
  candidate_keys <- sort(key_string(bundles[[definition_id]]$fig2f, c("canon", "species")))
  assert_true(
    identical(candidate_keys, d0_fig2f_keys),
    "%s Figure 2F displayed universe differs from D0", definition_id
  )
}

numeric_direction <- function(value, positive_label, negative_label, tolerance = 1e-12) {
  if (is.na(value)) return(NA_character_)
  if (value > tolerance) return(positive_label)
  if (value < -tolerance) return(negative_label)
  "equal_or_zero"
}

rows <- list()
append_row <- function(
  figure_panel, item_id, definition_id, comparison, metric = NA_character_,
  canon = NA_character_, group = NA_character_, value_1_label = NA_character_,
  value_1 = NA_real_, value_2_label = NA_character_, value_2 = NA_real_,
  effect_label = "value_1_minus_value_2", effect_value = NA_real_,
  statistic_label = NA_character_, statistic_value = NA_real_,
  boundary_rule = NA_character_, boundary_pass = NA,
  observed_direction = NA_character_, evidence_scope,
  source_file_1, source_file_2 = NA_character_, note
) {
  rows[[length(rows) + 1L]] <<- data.frame(
    figure_panel = figure_panel,
    item_id = item_id,
    definition_id = definition_id,
    definition_label = unname(definition_labels[[definition_id]]),
    comparison = comparison,
    metric = metric,
    canon = canon,
    group = group,
    value_1_label = value_1_label,
    value_1 = as.numeric(value_1),
    value_2_label = value_2_label,
    value_2 = as.numeric(value_2),
    effect_label = effect_label,
    effect_value = as.numeric(effect_value),
    statistic_label = statistic_label,
    statistic_value = as.numeric(statistic_value),
    boundary_rule = boundary_rule,
    boundary_pass = as.logical(boundary_pass),
    observed_direction = observed_direction,
    evidence_scope = evidence_scope,
    source_file_1 = source_file_1,
    source_file_2 = source_file_2,
    note = note,
    stringsAsFactors = FALSE
  )
}

for (definition_id in names(bundles)) {
  bundle <- bundles[[definition_id]]

  # Figure 2 D/E: only the three prespecified direct contrasts audited in the
  # response text.  Means and existing Welch P values come from separate,
  # finalized bottom tables and are joined by metric/canon.
  direct_specs <- data.frame(
    metric = c("richness", "shannon", "shannon"),
    canon = c("B3", "B3", "B4")
  )
  for (index in seq_len(nrow(direct_specs))) {
    metric <- direct_specs$metric[[index]]
    canon <- direct_specs$canon[[index]]
    means <- bundle$fig4b[
      bundle$fig4b$metric == metric & bundle$fig4b$canon == canon,
    ]
    hdr_mean <- means$current_mean[means$group == "HDR"]
    ldr_mean <- means$current_mean[means$group == "LDR"]
    p_value <- bundle$fig2_trajectory$p_welch_raw[
      bundle$fig2_trajectory$metric == metric &
        bundle$fig2_trajectory$canon == canon
    ]
    effect <- hdr_mean - ldr_mean
    append_row(
      figure_panel = "Fig2D/E",
      item_id = sprintf("Fig2_direct_HDR_vs_LDR__%s__%s", metric, canon),
      definition_id = definition_id,
      comparison = "HDR versus LDR",
      metric = metric,
      canon = canon,
      value_1_label = "HDR_mean",
      value_1 = hdr_mean,
      value_2_label = "LDR_mean",
      value_2 = ldr_mean,
      effect_label = "HDR_minus_LDR_mean",
      effect_value = effect,
      statistic_label = "Welch_P_unadjusted",
      statistic_value = p_value,
      boundary_rule = "Welch_P_unadjusted < 0.05",
      boundary_pass = p_value < 0.05,
      observed_direction = numeric_direction(effect, "HDR>LDR", "HDR<LDR"),
      evidence_scope = "inferential_existing_bottom_table",
      source_file_1 = bundle$paths[["fig4b"]],
      source_file_2 = bundle$paths[["fig2_trajectory"]],
      note = paste(
        "Means and the existing Welch P are summarized without recomputing",
        "the test."
      )
    )
  }

  # Figure 2C: retain both within-group positive-association evidence and the
  # purely descriptive HDR-versus-LDR ordering of rho.  The latter has no
  # formal between-correlation test and must not be described as inferential.
  for (current_group in groups) {
    record <- bundle$fig2_correlation[bundle$fig2_correlation$group == current_group, ]
    append_row(
      figure_panel = "Fig2C",
      item_id = sprintf("Fig2C_positive_baseline_endpoint__%s", current_group),
      definition_id = definition_id,
      comparison = sprintf("%s baseline versus endpoint Shannon", current_group),
      metric = "shannon",
      group = current_group,
      value_1_label = "Spearman_rho",
      value_1 = record$rho,
      value_2_label = "zero",
      value_2 = 0,
      effect_label = "rho",
      effect_value = record$rho,
      statistic_label = "Spearman_P",
      statistic_value = record$p,
      boundary_rule = "rho > 0 and Spearman_P < 0.05",
      boundary_pass = record$rho > 0 && record$p < 0.05,
      observed_direction = numeric_direction(record$rho, "positive", "negative"),
      evidence_scope = "inferential_existing_bottom_table",
      source_file_1 = bundle$paths[["fig2_correlation"]],
      note = "Within-group rho and P are copied from the finalized Figure 2C table."
    )
  }
  hdr_corr <- bundle$fig2_correlation$rho[bundle$fig2_correlation$group == "HDR"]
  ldr_corr <- bundle$fig2_correlation$rho[bundle$fig2_correlation$group == "LDR"]
  corr_difference <- hdr_corr - ldr_corr
  append_row(
    figure_panel = "Fig2C",
    item_id = "Fig2C_rho_point_estimate_order_HDR_vs_LDR",
    definition_id = definition_id,
    comparison = "descriptive ordering of within-group Spearman rho",
    metric = "shannon",
    value_1_label = "HDR_rho",
    value_1 = hdr_corr,
    value_2_label = "LDR_rho",
    value_2 = ldr_corr,
    effect_label = "HDR_rho_minus_LDR_rho",
    effect_value = corr_difference,
    statistic_label = NA_character_,
    statistic_value = NA_real_,
    boundary_rule = "descriptive only; no formal test of rho difference",
    boundary_pass = NA,
    observed_direction = numeric_direction(
      corr_difference, "HDR_rho>LDR_rho", "HDR_rho<LDR_rho"
    ),
    evidence_scope = "descriptive_point_estimate_only",
    source_file_1 = bundle$paths[["fig2_correlation"]],
    note = "No between-correlation inferential test is introduced by this audit."
  )

  # Figure 4B: include every non-A1 LDR cell so the disturbance-period and
  # endpoint statements are auditable without selecting only favorable cells.
  ldr_fig4b <- bundle$fig4b[
    bundle$fig4b$group == "LDR" & bundle$fig4b$canon != "A1",
  ]
  for (index in seq_len(nrow(ldr_fig4b))) {
    record <- ldr_fig4b[index, ]
    difference <- record$paired_mean_difference
    append_row(
      figure_panel = "Fig4B",
      item_id = sprintf("Fig4B_LDR_vs_A1__%s__%s", record$metric, record$canon),
      definition_id = definition_id,
      comparison = "LDR current time point versus paired A1",
      metric = record$metric,
      canon = record$canon,
      group = "LDR",
      value_1_label = "current_mean",
      value_1 = record$current_mean,
      value_2_label = "baseline_A1_mean",
      value_2 = record$baseline_A1_mean,
      effect_label = "paired_mean_difference_current_minus_A1",
      effect_value = difference,
      statistic_label = "BH_q_within_metric_group",
      statistic_value = record$q_BH_within_metric_group,
      boundary_rule = "BH_q_within_metric_group < 0.05",
      boundary_pass = record$q_BH_within_metric_group < 0.05,
      observed_direction = numeric_direction(
        difference, "current>A1", "current<A1"
      ),
      evidence_scope = "inferential_existing_bottom_table",
      source_file_1 = bundle$paths[["fig4b"]],
      note = "Paired difference and within-metric/group BH q are copied directly."
    )
  }

  # Figure 4C endpoint: compare complete caption-qualified eligible networks.
  endpoint <- bundle$fig4c[bundle$fig4c$stage == "End", ]
  for (network_measure in c("nodes", "edges")) {
    column <- if (network_measure == "nodes") {
      "eligible_nodes_caption_rule"
    } else {
      "eligible_edges_caption_rule"
    }
    hdr_value <- endpoint[[column]][endpoint$group == "HDR"]
    ldr_value <- endpoint[[column]][endpoint$group == "LDR"]
    difference <- hdr_value - ldr_value
    append_row(
      figure_panel = "Fig4C",
      item_id = sprintf("Fig4C_End_HDR_vs_LDR__%s", network_measure),
      definition_id = definition_id,
      comparison = "pooled C4-C6 endpoint HDR versus LDR network count",
      metric = network_measure,
      canon = "C4-C6 pooled End",
      value_1_label = sprintf("HDR_%s", network_measure),
      value_1 = hdr_value,
      value_2_label = sprintf("LDR_%s", network_measure),
      value_2 = ldr_value,
      effect_label = sprintf("HDR_minus_LDR_%s", network_measure),
      effect_value = difference,
      statistic_label = "HDR_to_LDR_count_ratio",
      statistic_value = hdr_value / ldr_value,
      boundary_rule = "descriptive counts after |rho|>0.6 and bootstrap raw P<0.05",
      boundary_pass = NA,
      observed_direction = numeric_direction(
        difference, "HDR_count>LDR_count", "HDR_count<LDR_count"
      ),
      evidence_scope = "descriptive_qualified_network_counts",
      source_file_1 = bundle$paths[["fig4c"]],
      note = paste(
        "Counts use the full caption-qualified network; this audit does not",
        "treat the count contrast as an inferential HDR-LDR test."
      )
    )
  }

  # Figure 4D: include the complete LDR Faecalibacterium curve for each metric;
  # the response can then transparently identify which time points reproduce.
  ldr_faec <- bundle$fig4d[
    bundle$fig4d$group == "LDR" &
      bundle$fig4d$genus == "Faecalibacterium",
  ]
  for (index in seq_len(nrow(ldr_faec))) {
    record <- ldr_faec[index, ]
    append_row(
      figure_panel = "Fig4D",
      item_id = sprintf(
        "Fig4D_LDR_Faecalibacterium__%s__%s", record$metric, record$canon
      ),
      definition_id = definition_id,
      comparison = "LDR Faecalibacterium abundance versus alpha diversity",
      metric = record$metric,
      canon = record$canon,
      group = "LDR",
      value_1_label = "Spearman_rho",
      value_1 = record$rho,
      value_2_label = "zero",
      value_2 = 0,
      effect_label = "rho",
      effect_value = record$rho,
      statistic_label = "BH_q_within_curve_sensitivity",
      statistic_value = record$q_BH_within_curve_sensitivity,
      boundary_rule = "rho > 0 and BH_q_within_curve_sensitivity < 0.05",
      boundary_pass = record$rho > 0 && record$q_BH_within_curve_sensitivity < 0.05,
      observed_direction = numeric_direction(record$rho, "positive", "negative"),
      evidence_scope = "inferential_existing_bottom_table",
      source_file_1 = bundle$paths[["fig4d"]],
      note = "Rho and within-curve BH sensitivity q are copied directly."
    )
  }
}

recurrence <- do.call(rbind, rows)
assert_unique_keys(
  recurrence, c("item_id", "definition_id"), "combined Figure 2/4 recurrence table"
)
assert_true(
  all(table(recurrence$item_id) == 7L),
  "Every recurrence item must contain exactly D0-D6"
)

# Attach the D0 direction/boundary and cross-definition claim boundaries.
recurrence$d0_observed_direction <- NA_character_
recurrence$direction_matches_D0 <- NA
recurrence$d0_boundary_pass <- NA
recurrence$boundary_status_matches_D0 <- NA
recurrence$all_D1_D6_direction_match_D0 <- NA
recurrence$all_D1_D6_boundary_pass <- NA
recurrence$safe_to_call_D0_reproduced_this_definition <- NA
recurrence$safe_to_call_D0_reproduced_across_all_D1_D6 <- NA
recurrence$claim_category <- NA_character_

for (item_id in unique(recurrence$item_id)) {
  selected <- recurrence$item_id == item_id
  d0_index <- which(selected & recurrence$definition_id == "D0")
  alt_indices <- which(selected & recurrence$definition_id != "D0")
  assert_true(length(d0_index) == 1L, "%s does not have exactly one D0 row", item_id)
  assert_true(length(alt_indices) == 6L, "%s does not have six alternative rows", item_id)

  d0_direction <- recurrence$observed_direction[d0_index]
  direction_matches <- recurrence$observed_direction[selected] == d0_direction
  recurrence$d0_observed_direction[selected] <- d0_direction
  recurrence$direction_matches_D0[selected] <- direction_matches

  d0_boundary <- recurrence$boundary_pass[d0_index]
  recurrence$d0_boundary_pass[selected] <- d0_boundary
  if (!is.na(d0_boundary)) {
    recurrence$boundary_status_matches_D0[selected] <-
      recurrence$boundary_pass[selected] == d0_boundary
  }

  all_direction <- all(recurrence$direction_matches_D0[alt_indices])
  boundary_values <- recurrence$boundary_pass[c(d0_index, alt_indices)]
  all_boundary <- if (all(is.na(boundary_values))) {
    NA
  } else {
    assert_true(
      !anyNA(boundary_values),
      "%s mixes applicable and non-applicable boundary flags", item_id
    )
    all(recurrence$boundary_pass[alt_indices])
  }
  recurrence$all_D1_D6_direction_match_D0[selected] <- all_direction
  recurrence$all_D1_D6_boundary_pass[selected] <- all_boundary

  if (all(is.na(boundary_values))) {
    safe_all <- all_direction
  } else {
    safe_all <- isTRUE(d0_boundary) && all_direction && isTRUE(all_boundary)
  }
  recurrence$safe_to_call_D0_reproduced_across_all_D1_D6[selected] <- safe_all

  for (index in c(d0_index, alt_indices)) {
    if (recurrence$definition_id[index] == "D0") {
      recurrence$claim_category[index] <- "D0_reference"
      next
    }
    if (all(is.na(boundary_values))) {
      safe_current <- isTRUE(recurrence$direction_matches_D0[index])
      recurrence$safe_to_call_D0_reproduced_this_definition[index] <- safe_current
      if (recurrence$evidence_scope[index] == "descriptive_point_estimate_only") {
        recurrence$claim_category[index] <- if (safe_current) {
          "descriptive_D0_order_reproduced"
        } else {
          "six_alternatives_consistent_but_D0_order_not_reproduced"
        }
      } else {
        recurrence$claim_category[index] <- if (safe_current) {
          "descriptive_D0_direction_reproduced"
        } else {
          "descriptive_D0_direction_not_reproduced"
        }
      }
    } else {
      safe_current <-
        isTRUE(d0_boundary) &&
        isTRUE(recurrence$direction_matches_D0[index]) &&
        isTRUE(recurrence$boundary_pass[index])
      recurrence$safe_to_call_D0_reproduced_this_definition[index] <- safe_current
      recurrence$claim_category[index] <- if (safe_current) {
        "D0_direction_and_boundary_reproduced"
      } else if (isTRUE(recurrence$direction_matches_D0[index])) {
        "D0_direction_reproduced_but_boundary_not_reproduced"
      } else {
        "D0_direction_not_reproduced"
      }
    }
  }
}

# Figure 2F direction/significance agreement, using a fixed sorted cell key.
sort_fig2f <- function(frame) {
  frame[order(match(frame$canon, time_order), frame$species), ]
}
fig2f_sorted <- lapply(bundles, function(bundle) sort_fig2f(bundle$fig2f))
d0_f2f <- fig2f_sorted$D0
d0_sign <- sign(d0_f2f$lda_score_LDR_over_HDR)
d0_significant <- d0_f2f$displayed_significant
assert_true(sum(d0_significant) == 352L, "D0 Figure 2F significant-cell count drifted")

alt_ids <- paste0("D", 1:6)
alt_sign_matrix <- do.call(
  cbind,
  lapply(fig2f_sorted[alt_ids], function(frame) sign(frame$lda_score_LDR_over_HDR))
)
colnames(alt_sign_matrix) <- alt_ids
alt_sig_matrix <- do.call(
  cbind,
  lapply(fig2f_sorted[alt_ids], function(frame) frame$displayed_significant)
)
colnames(alt_sig_matrix) <- alt_ids

universal_significant <- rowSums(alt_sig_matrix) == length(alt_ids)
unanimous_direction <- apply(alt_sign_matrix, 1L, function(values) length(unique(values)) == 1L)
unanimous_matches_d0 <- unanimous_direction & alt_sign_matrix[, 1L] == d0_sign
majority_sum <- rowSums(alt_sign_matrix)
majority_non_tie <- majority_sum != 0
majority_matches_d0 <- majority_non_tie & sign(majority_sum) == d0_sign

fig2f_rows <- list()
for (definition_id in names(fig2f_sorted)) {
  frame <- fig2f_sorted[[definition_id]]
  current_sign <- sign(frame$lda_score_LDR_over_HDR)
  current_sig <- frame$displayed_significant
  direction_agree <- current_sign == d0_sign
  joint_significant <- current_sig & d0_significant
  fig2f_rows[[length(fig2f_rows) + 1L]] <- data.frame(
    definition_id = definition_id,
    definition_label = unname(definition_labels[[definition_id]]),
    n_displayed_species = length(unique(frame$species)),
    n_canonical_timepoints = length(unique(frame$canon)),
    n_displayed_cells = nrow(frame),
    direction_definition = "sign(lda_score_LDR_over_HDR): positive=LDR-enriched; negative=HDR-enriched",
    significance_definition = paste(
      "exact Mann-Whitney P, BH-adjusted among 45 displayed species",
      "within each canonical timepoint; q<0.05"
    ),
    d0_n_significant_cells = sum(d0_significant),
    n_direction_agree_with_D0_all_cells = sum(direction_agree),
    proportion_direction_agree_with_D0_all_cells = mean(direction_agree),
    n_direction_disagree_with_D0_all_cells = sum(!direction_agree),
    n_direction_agree_with_D0_among_D0_significant = sum(direction_agree[d0_significant]),
    proportion_direction_agree_with_D0_among_D0_significant = mean(direction_agree[d0_significant]),
    n_definition_significant_cells = sum(current_sig),
    n_joint_D0_and_definition_significant_cells = sum(joint_significant),
    n_joint_significant_and_same_direction = sum(joint_significant & direction_agree),
    n_cells_significant_in_all_D1_D6 = sum(universal_significant),
    n_cells_with_unanimous_D1_D6_direction = sum(unanimous_direction),
    n_unanimous_D1_D6_direction_cells_matching_D0 = sum(unanimous_matches_d0),
    n_cells_with_non_tied_D1_D6_majority_direction = sum(majority_non_tie),
    n_D1_D6_majority_direction_cells_matching_D0 = sum(majority_matches_d0),
    safe_to_claim_universal_taxon_level_replication = FALSE,
    safe_claim = if (definition_id == "D0") {
      "D0 reference"
    } else {
      "broad directional agreement only; not universal BH-significant replication"
    },
    source_file = bundles[[definition_id]]$paths[["fig2f"]],
    stringsAsFactors = FALSE
  )
}
fig2f_agreement <- do.call(rbind, fig2f_rows)

# Regression tests for the exact audited conclusions.  These are deliberately
# claim-level checks, not new analyses.
required_reproduced <- c(
  "Fig2_direct_HDR_vs_LDR__richness__B3",
  "Fig2_direct_HDR_vs_LDR__shannon__B3",
  "Fig2_direct_HDR_vs_LDR__shannon__B4",
  "Fig4C_End_HDR_vs_LDR__nodes",
  "Fig4C_End_HDR_vs_LDR__edges",
  paste0("Fig4B_LDR_vs_A1__richness__", c("B3", "B4", "B5", "B6", "C2")),
  paste0("Fig4B_LDR_vs_A1__shannon__", c("B3", "B4", "B5", "B6")),
  paste0("Fig4D_LDR_Faecalibacterium__shannon__", c("B2", "B3")),
  "Fig4D_LDR_Faecalibacterium__simpson__B2"
)
for (item_id in required_reproduced) {
  values <- unique(recurrence$safe_to_call_D0_reproduced_across_all_D1_D6[
    recurrence$item_id == item_id
  ])
  assert_true(
    identical(values, TRUE),
    "Required reproduced claim failed its audit gate: %s", item_id
  )
}
rho_order <- unique(recurrence$safe_to_call_D0_reproduced_across_all_D1_D6[
  recurrence$item_id == "Fig2C_rho_point_estimate_order_HDR_vs_LDR"
])
assert_true(
  identical(rho_order, FALSE),
  "Figure 2C relative rho ordering must not be labeled D0-reproduced"
)
required_not_universally_reproduced <- c(
  paste0("Fig4B_LDR_vs_A1__simpson__", c("B3", "B4")),
  paste0("Fig4B_LDR_vs_A1__shannon__", c("C4", "C5", "C6")),
  paste0("Fig4D_LDR_Faecalibacterium__richness__", c("B2", "B3"))
)
for (item_id in required_not_universally_reproduced) {
  values <- unique(recurrence$safe_to_call_D0_reproduced_across_all_D1_D6[
    recurrence$item_id == item_id
  ])
  assert_true(
    identical(values, FALSE),
    "Non-universal claim was incorrectly labeled reproduced: %s", item_id
  )
}
expected_all_agreement <- c(D0 = 720L, D1 = 520L, D2 = 533L, D3 = 654L,
                            D4 = 427L, D5 = 556L, D6 = 659L)
observed_all_agreement <- stats::setNames(
  fig2f_agreement$n_direction_agree_with_D0_all_cells,
  fig2f_agreement$definition_id
)
assert_true(
  identical(observed_all_agreement[names(expected_all_agreement)], expected_all_agreement),
  "Figure 2F all-cell D0 direction-agreement counts drifted"
)
assert_true(
  all(fig2f_agreement$n_cells_significant_in_all_D1_D6 == 0L),
  "Figure 2F unexpectedly has a cell significant in all D1-D6"
)
assert_true(
  all(fig2f_agreement$n_cells_with_unanimous_D1_D6_direction == 326L) &&
    all(fig2f_agreement$n_unanimous_D1_D6_direction_cells_matching_D0 == 312L),
  "Figure 2F unanimous-direction audit counts drifted"
)

# Stable ordering makes diffs and downstream parsing deterministic.
panel_order <- c("Fig2D/E", "Fig2C", "Fig4B", "Fig4C", "Fig4D")
recurrence <- recurrence[
  order(
    match(recurrence$figure_panel, panel_order),
    recurrence$item_id,
    match(recurrence$definition_id, names(definition_slugs))
  ),
]
row.names(recurrence) <- NULL
fig2f_agreement <- fig2f_agreement[
  match(names(definition_slugs), fig2f_agreement$definition_id),
]
row.names(fig2f_agreement) <- NULL

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
recurrence_path <- file.path(output_root, "Fig24_D0_vs_D1_D6_recurrence.csv")
fig2f_path <- file.path(output_root, "Fig2F_D0_direction_agreement.csv")
utils::write.csv(recurrence, recurrence_path, row.names = FALSE, na = "")
utils::write.csv(fig2f_agreement, fig2f_path, row.names = FALSE, na = "")

cat(sprintf(
  paste0(
    "PASS: validated D0-D6 bottom-table keys and wrote %d recurrence rows to %s\n",
    "PASS: wrote %d Figure 2F agreement rows to %s\n"
  ),
  nrow(recurrence), recurrence_path,
  nrow(fig2f_agreement), fig2f_path
))
