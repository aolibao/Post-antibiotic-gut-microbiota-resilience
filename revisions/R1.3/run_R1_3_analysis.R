#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(base_dir, "tables")
figure_dir <- file.path(base_dir, "figures")
log_dir <- file.path(base_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

z <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

hc3_vcov <- function(fit) {
  x <- model.matrix(fit)
  h <- hatvalues(fit)
  if (any(h >= 1 - 1e-10)) stop("HC3 undefined: leverage too close to 1")
  bread <- solve(crossprod(x))
  meat <- crossprod(x, x * (residuals(fit)^2 / (1 - h)^2))
  out <- bread %*% meat %*% bread
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

hc3_tidy <- function(fit, model_id) {
  v <- hc3_vcov(fit)
  b <- coef(fit)
  se <- sqrt(diag(v))
  df <- df.residual(fit)
  crit <- qt(0.975, df)
  tibble(
    model_id = model_id,
    term = names(b),
    estimate = as.numeric(b),
    HC3_SE = se,
    statistic = as.numeric(b / se),
    df_residual = df,
    p_value = 2 * pt(abs(b / se), df = df, lower.tail = FALSE),
    CI95_low = as.numeric(b - crit * se),
    CI95_high = as.numeric(b + crit * se),
    n_subjects = nobs(fit),
    R2 = summary(fit)$r.squared,
    adjusted_R2 = summary(fit)$adj.r.squared
  )
}

robust_block_test <- function(fit, pattern, block_label) {
  terms <- grep(pattern, names(coef(fit)), value = TRUE)
  if (length(terms) == 0L) stop("No terms found for block: ", block_label)
  b <- coef(fit)[terms]
  v <- hc3_vcov(fit)[terms, terms, drop = FALSE]
  q <- length(terms)
  f_value <- as.numeric(t(b) %*% solve(v, b) / q)
  tibble(
    model_id = deparse(formula(fit)),
    block = block_label,
    terms_tested = paste(terms, collapse = "; "),
    numerator_df = q,
    denominator_df = df.residual(fit),
    HC3_Wald_F = f_value,
    p_value = pf(f_value, q, df.residual(fit), lower.tail = FALSE)
  )
}

input_path <- file.path(root_dir, "reviewer_response", "03_tables", "core_subject_recovery_metrics.csv")
d <- read_csv(input_path, show_col_types = FALSE) |>
  filter(dataset_branch == "rarefied_primary") |>
  arrange(subject) |>
  mutate(
    recovery_score_D5 = -recovery_auc_bray_actual_day,
    z_recovery = z(recovery_score_D5),
    z_shannon = z(baseline_shannon),
    z_bp = z(baseline_bacteroides_plus_phocaeicola_fraction),
    bp_dominance = factor(
      if_else(bacteroides_plus_phocaeicola_dominant,
              "Bacteroides/Phocaeicola dominant", "Other genus dominant"),
      levels = c("Other genus dominant", "Bacteroides/Phocaeicola dominant")
    ),
    enterotype_k3 = factor(exploratory_enterotype_k3),
    endpoint_group = factor(original_endpoint_group, levels = c("LDR", "Middle", "HDR"))
  )

stopifnot(nrow(d) == 60L, n_distinct(d$subject) == 60L)

fits <- list(
  shannon_only = lm(z_recovery ~ z_shannon, data = d),
  shannon_and_bp_continuous = lm(z_recovery ~ z_shannon + z_bp, data = d),
  shannon_and_bp_dominance = lm(z_recovery ~ z_shannon + bp_dominance, data = d),
  shannon_and_k3_cluster = lm(z_recovery ~ z_shannon + enterotype_k3, data = d)
)

coef_table <- bind_rows(lapply(names(fits), function(nm) hc3_tidy(fits[[nm]], nm))) |>
  mutate(
    outcome = "standardized negative D5 Bray-Curtis AUC (higher = better recovery)",
    inference = "OLS with HC3 covariance; participant is the sampling unit"
  )
write_csv(coef_table, file.path(table_dir, "R1_3_conditional_models_HC3.csv"))

block_tests <- bind_rows(
  robust_block_test(fits$shannon_and_bp_dominance, "^bp_dominance", "baseline Bacteroides/Phocaeicola dominance | Shannon"),
  robust_block_test(fits$shannon_and_k3_cluster, "^enterotype_k3", "exploratory k=3 CLR cluster | Shannon")
)
write_csv(block_tests, file.path(table_dir, "R1_3_categorical_block_tests.csv"))

classification_counts <- bind_rows(
  d |> count(classification = "Exact dominant genus", category = dominant_genus, name = "subjects_n"),
  d |> count(classification = "Bacteroides/Phocaeicola dominance", category = as.character(bp_dominance), name = "subjects_n"),
  d |> count(classification = "Exploratory k=3 CLR cluster", category = as.character(enterotype_k3), name = "subjects_n")
)
write_csv(classification_counts, file.path(table_dir, "R1_3_baseline_classification_counts.csv"))

extremes <- d |> filter(endpoint_group %in% c("LDR", "HDR")) |> droplevels()
cross_bp <- with(extremes, table(endpoint_group, bp_dominance))
cross_k3 <- with(extremes, table(endpoint_group, enterotype_k3))
cross_table <- bind_rows(
  as.data.frame(cross_bp) |> transmute(classification = "Bacteroides/Phocaeicola dominance", endpoint_group = as.character(endpoint_group), category = as.character(bp_dominance), subjects_n = Freq),
  as.data.frame(cross_k3) |> transmute(classification = "Exploratory k=3 CLR cluster", endpoint_group = as.character(endpoint_group), category = as.character(enterotype_k3), subjects_n = Freq)
)
write_csv(cross_table, file.path(table_dir, "R1_3_original_HDR_LDR_crossclassification.csv"))

fisher_table <- tibble(
  classification = c("Bacteroides/Phocaeicola dominance", "Exploratory k=3 CLR cluster"),
  test = "Fisher exact test",
  n_subjects = 30L,
  p_value = c(fisher.test(cross_bp)$p.value, fisher.test(cross_k3, simulate.p.value = FALSE)$p.value),
  note = c(
    "Original endpoint-defined HDR/LDR only; descriptive association audit",
    "Original endpoint-defined HDR/LDR only; exploratory cluster, not a validated enterotype"
  )
)
write_csv(fisher_table, file.path(table_dir, "R1_3_original_HDR_LDR_classification_tests.csv"))

distance_audit_path <- file.path(root_dir, "reviewer_response", "03_tables", "core_distance_PERMANOVA_PERMDISP.csv")
baseline_composition_audit <- read_csv(distance_audit_path, show_col_types = FALSE) |>
  filter(dataset_branch == "rarefied_primary", stage == "baseline", comparison == "original_HDR_vs_LDR")
write_csv(baseline_composition_audit, file.path(table_dir, "R1_3_original_HDR_LDR_baseline_composition_audit.csv"))

# Figure source data are intentionally explicit.
figure_source <- d |>
  select(subject, baseline_shannon, recovery_score_D5, z_shannon, z_recovery,
         baseline_bacteroides_plus_phocaeicola_fraction, bp_dominance,
         enterotype_k3, dominant_genus, endpoint_group)
write_csv(figure_source, file.path(table_dir, "R1_3_figure_source_data.csv"))

blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"; grey <- "#666666"

# Publication typography: retain only essential labels, enlarge all type, and
# force every text element to pure black for consistent PDF/TIFF rendering.
theme_pub <- theme_classic(base_size = 11.5, base_family = "Liberation Sans") +
  theme(
    text = element_text(colour = "#000000"),
    axis.title = element_text(colour = "#000000", size = 11.5),
    axis.text = element_text(colour = "#000000", size = 10),
    axis.line = element_line(colour = "#000000", linewidth = 0.55),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.45),
    legend.text = element_text(colour = "#000000", size = 9.5),
    legend.title = element_text(colour = "#000000", size = 10),
    plot.tag = element_text(colour = "#000000", face = "bold", size = 15),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )
theme_set(theme_pub)

p_a <- ggplot(d, aes(baseline_shannon, recovery_score_D5, color = bp_dominance, shape = bp_dominance)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(data = d, aes(baseline_shannon, recovery_score_D5), method = "lm", se = TRUE,
              color = "#222222", fill = "#BDBDBD", linewidth = 0.7, inherit.aes = FALSE) +
  scale_color_manual(
    values = c("Other genus dominant" = blue, "Bacteroides/Phocaeicola dominant" = orange),
    breaks = c("Other genus dominant", "Bacteroides/Phocaeicola dominant")
  ) +
  scale_shape_manual(
    values = c("Other genus dominant" = 16, "Bacteroides/Phocaeicola dominant" = 17),
    breaks = c("Other genus dominant", "Bacteroides/Phocaeicola dominant")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(x = "Baseline Shannon (A1-A4 mean)", y = "D5 recovery score", color = NULL, shape = NULL) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 3),
    legend.key.width = grid::unit(0.7, "lines"),
    legend.spacing.x = grid::unit(0.25, "lines")
  )

forest <- bind_rows(
  coef_table |> filter(model_id == "shannon_and_bp_continuous", term == "z_shannon") |> mutate(label = "Shannon | Bacteroides + Phocaeicola abundance"),
  coef_table |> filter(model_id == "shannon_and_bp_continuous", term == "z_bp") |> mutate(label = "Bacteroides + Phocaeicola abundance | Shannon"),
  coef_table |> filter(model_id == "shannon_and_bp_dominance", term == "z_shannon") |> mutate(label = "Shannon | dominance class"),
  coef_table |> filter(model_id == "shannon_and_bp_dominance", grepl("^bp_dominance", term)) |> mutate(label = "Bacteroides/Phocaeicola dominant | Shannon"),
  coef_table |> filter(model_id == "shannon_and_k3_cluster", term == "z_shannon") |> mutate(label = "Shannon | exploratory k=3 cluster")
) |>
  mutate(label = factor(label, levels = rev(unique(label))), sig = p_value < 0.05)

p_b <- ggplot(forest, aes(estimate, label)) +
  geom_vline(xintercept = 0, color = "#888888", linewidth = 0.45, linetype = 2) +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.16, linewidth = 0.65, color = grey) +
  geom_point(aes(color = sig), size = 2.4) +
  scale_color_manual(values = c(`FALSE` = blue, `TRUE` = orange), guide = "none") +
  scale_y_discrete(labels = function(x) {
    x <- gsub("Bacteroides \\+ Phocaeicola abundance", "Bacteroides +\nPhocaeicola abundance", x)
    x <- gsub("Bacteroides/Phocaeicola dominant", "Bacteroides/Phocaeicola\ndominant", x, fixed = TRUE)
    x
  }) +
  labs(x = "Conditional coefficient\n(HC3 95% CI)", y = NULL) +
  theme(
    axis.text.y = element_text(size = 9.5, lineheight = 0.95),
    plot.margin = margin(6, 6, 6, 10)
  )

p_c <- ggplot(d, aes(bp_dominance, recovery_score_D5, color = bp_dominance, shape = bp_dominance)) +
  geom_boxplot(aes(group = bp_dominance), width = 0.46, outlier.shape = NA, color = "#555555", fill = "white", linewidth = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c("Other genus dominant" = blue, "Bacteroides/Phocaeicola dominant" = orange), guide = "none") +
  scale_shape_manual(values = c("Other genus dominant" = 16, "Bacteroides/Phocaeicola dominant" = 17), guide = "none") +
  scale_x_discrete(labels = c(
    "Other genus dominant" = "Other genus\ndominant",
    "Bacteroides/Phocaeicola dominant" = "Bacteroides/Phocaeicola\ndominant"
  )) +
  labs(x = NULL, y = "D5 recovery score") +
  theme(axis.text.x = element_text(size = 10, lineheight = 0.95))

cross_plot <- cross_table |>
  filter(classification == "Bacteroides/Phocaeicola dominance") |>
  group_by(endpoint_group) |>
  mutate(proportion = subjects_n / sum(subjects_n)) |>
  ungroup()
p_d <- ggplot(cross_plot, aes(endpoint_group, proportion, fill = category)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.4) +
  geom_text(aes(label = subjects_n), position = position_stack(vjust = 0.5),
            size = 3.8, colour = "#000000", family = "Liberation Sans") +
  annotate(
    "text", x = 1.5, y = 1.075,
    label = paste0("P = ", formatC(fisher.test(cross_bp)$p.value, digits = 3, format = "f")),
    colour = "#000000", family = "Liberation Sans", size = 4.0
  ) +
  scale_fill_manual(
    values = c("Other genus dominant" = blue, "Bacteroides/Phocaeicola dominant" = orange),
    breaks = c("Other genus dominant", "Bacteroides/Phocaeicola dominant")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.12),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Endpoint-Shannon group", y = "Proportion", fill = NULL) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 3),
    legend.key.width = grid::unit(0.7, "lines"),
    legend.spacing.x = grid::unit(0.25, "lines")
  )

p_a <- p_a + labs(tag = "A")
p_b <- p_b + labs(tag = "B")
p_c <- p_c + labs(tag = "C")
p_d <- p_d + labs(tag = "D")

# Draw with base grid so this script remains reproducible in the project's
# pinned R environment without adding another package dependency.
draw_combined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 2, ncol = 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(1, 1), "null")
  )))
  print(p_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1), newpage = FALSE)
  print(p_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2), newpage = FALSE)
  print(p_c, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
  print(p_d, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2), newpage = FALSE)
  grid::popViewport()
}

pdf_path <- file.path(figure_dir, "R1_3_enterotype_conditional_analysis.pdf")
tif_path <- file.path(figure_dir, "R1_3_enterotype_conditional_analysis.tif")
png_path <- file.path(figure_dir, "R1_3_enterotype_conditional_analysis_preview.png")

grDevices::cairo_pdf(pdf_path, width = 210 / 25.4, height = 185 / 25.4,
                     family = "Liberation Sans", bg = "white")
draw_combined()
grDevices::dev.off()

grDevices::tiff(tif_path, width = 210, height = 185, units = "mm", res = 600,
                compression = "lzw", bg = "white", type = "cairo",
                family = "Liberation Sans")
draw_combined()
grDevices::dev.off()

grDevices::png(png_path, width = 210, height = 185, units = "mm", res = 160,
               bg = "white", type = "cairo", family = "Liberation Sans")
draw_combined()
grDevices::dev.off()

capture.output(sessionInfo(), file = file.path(log_dir, "R1_3_R_sessionInfo.txt"))

key <- coef_table |>
  filter(model_id == "shannon_and_bp_continuous", term %in% c("z_shannon", "z_bp")) |>
  select(term, estimate, CI95_low, CI95_high, p_value)
write_lines(c(
  "R1.3 completed",
  paste0("Participants: ", nrow(d)),
  paste0("Exploratory k=3 silhouette inherited from the audited core pipeline: 0.139"),
  paste(capture.output(print(key)), collapse = "\n"),
  paste(capture.output(print(block_tests)), collapse = "\n"),
  paste(capture.output(print(fisher_table)), collapse = "\n")
), file.path(log_dir, "R1_3_run_summary.txt"))

