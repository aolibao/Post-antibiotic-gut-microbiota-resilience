#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))] |> sub("--file=", "", x = _)), "..", "..", ".."))
src_tables <- file.path(root, "reviewer_response", "03_tables")
out <- file.path(root, "reviewer_response", "R1.6")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)

needed <- c(
  "core_sample_longitudinal_metrics.csv",
  "core_repeated_GEE_coefficients.csv",
  "core_repeated_GEE_omnibus.csv",
  "core_repeated_MixedLM_sensitivity.csv",
  "core_distance_PERMANOVA_PERMDISP.csv",
  "lefse_multiplicity_audit.csv"
)
stopifnot(all(file.exists(file.path(src_tables, needed))))

long <- read_csv(file.path(src_tables, needed[1]), show_col_types = FALSE)
gee_coef <- read_csv(file.path(src_tables, needed[2]), show_col_types = FALSE)
gee_omni <- read_csv(file.path(src_tables, needed[3]), show_col_types = FALSE)
mixed <- read_csv(file.path(src_tables, needed[4]), show_col_types = FALSE)
dist <- read_csv(file.path(src_tables, needed[5]), show_col_types = FALSE)
lefse <- read_csv(file.path(src_tables, needed[6]), show_col_types = FALSE)
phase_clr_file <- file.path(out, "tables", "R1_6_participant_phase_CLR_global_FDR.csv")
phase_clr_summary_file <- file.path(out, "tables", "R1_6_participant_phase_CLR_summary.csv")
stopifnot(file.exists(phase_clr_file), file.exists(phase_clr_summary_file))
phase_clr <- read_csv(phase_clr_file, show_col_types = FALSE)
phase_clr_summary <- read_csv(phase_clr_summary_file, show_col_types = FALSE)

write_csv(gee_coef, file.path(out, "tables", "R1_6_GEE_coefficients.csv"))
write_csv(gee_omni, file.path(out, "tables", "R1_6_GEE_omnibus.csv"))
write_csv(mixed, file.path(out, "tables", "R1_6_MixedLM_sensitivity.csv"))
write_csv(dist, file.path(out, "tables", "R1_6_PERMANOVA_PERMDISP.csv"))
write_csv(lefse, file.path(out, "tables", "R1_6_LEfSe_multiplicity_audit.csv"))

time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
phase_bg <- tibble(
  xmin = c(0.5, 4.5, 10.5), xmax = c(4.5, 10.5, 16.5),
  ymin = -Inf, ymax = Inf,
  phase = c("Baseline", "Antibiotic exposure", "Post-exposure follow-up")
)

set.seed(20261012)
traj0 <- long |>
  filter(dataset_branch == "rarefied_primary", original_endpoint_group %in% c("HDR", "LDR")) |>
  mutate(canon = factor(canon, levels = time_order))

boot_ci <- function(v, B = 2000L) {
  means <- replicate(B, mean(sample(v, length(v), replace = TRUE)))
  c(mean = mean(v), lo = unname(quantile(means, 0.025)), hi = unname(quantile(means, 0.975)))
}

traj <- traj0 |>
  group_by(original_endpoint_group, canon) |>
  summarise(n_subjects = n_distinct(subject), vals = list(shannon), .groups = "drop") |>
  rowwise() |>
  mutate(tmp = list(boot_ci(unlist(vals))), mean = tmp[1], ci_low = tmp[2], ci_high = tmp[3]) |>
  ungroup() |>
  select(group = original_endpoint_group, canon, n_subjects, mean, ci_low, ci_high)
write_csv(traj, file.path(out, "tables", "R1_6_Shannon_trajectory_summary.csv"))

summary_tbl <- bind_rows(
  gee_omni |> filter(dataset_branch == "rarefied_primary", outcome == "shannon", model == "original_HDR_LDR") |>
    transmute(section = "Repeated measures", result = "GEE HDR/LDR group x time", estimate = wald_chi2, df, p_value,
              note = interpretation),
  dist |> filter(dataset_branch == "rarefied_primary", distance == "Bray-Curtis", comparison == "all60_within_subject_phase") |>
    transmute(section = "Distance inference", result = "Paired baseline vs endpoint PERMANOVA", estimate = PERMANOVA_R2,
              df = NA_real_, p_value = PERMANOVA_p, note = paste0("PERMDISP P=", PERMDISP_p)),
  dist |> filter(dataset_branch == "rarefied_primary", distance == "binary Jaccard", comparison == "all60_within_subject_phase") |>
    transmute(section = "Distance inference", result = "Paired baseline vs endpoint PERMANOVA", estimate = PERMANOVA_R2,
              df = NA_real_, p_value = PERMANOVA_p, note = paste0("PERMDISP P=", PERMDISP_p))
)
write_csv(summary_tbl, file.path(out, "tables", "R1_6_key_results.csv"))

pal <- c(HDR = "#0072B2", LDR = "#D55E00")
base_theme <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, colour = "grey25"),
    legend.position = "top", legend.title = element_blank(),
    axis.text = element_text(colour = "black"),
    plot.margin = margin(5, 7, 5, 5)
  )

pA <- ggplot(traj, aes(as.numeric(canon), mean, colour = group, shape = group)) +
  geom_rect(data = phase_bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = c("#F3F4F6", "#FFF2D5", "#EAF3F8"), alpha = 0.55) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = group), colour = NA, alpha = 0.12) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.8, stroke = 0.35) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_shape_manual(values = c(HDR = 16, LDR = 17)) +
  scale_x_continuous(breaks = seq_along(time_order), labels = time_order, expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "A  Repeated Shannon measurements", subtitle = "Mean and participant-level bootstrap 95% CI; reference groups only",
       x = "Canonical sampling point", y = "Shannon diversity") + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

omni_plot <- gee_omni |>
  filter(dataset_branch == "rarefied_primary", model %in% c("all60_time", "original_HDR_LDR")) |>
  mutate(
    model_label = recode(model, all60_time = "All 60: time", original_HDR_LDR = "HDR/LDR: group x time"),
    outcome_label = recode(outcome,
      shannon = "Shannon",
      signed_shannon_deviation_from_baseline = "Signed Shannon deviation",
      bray_to_personal_baseline_centroid = "Bray-Curtis to own baseline",
      jaccard_to_personal_baseline_consensus = "Jaccard to own baseline"
    ),
    score = pmin(-log10(p_value), 30),
    selected = model == "original_HDR_LDR"
  )
pB <- ggplot(omni_plot, aes(score, reorder(outcome_label, score), colour = model_label, shape = model_label)) +
  geom_vline(xintercept = -log10(0.05), linetype = 3, colour = "grey55", linewidth = 0.45) +
  geom_point(size = 2.5, stroke = 0.5) +
  scale_colour_manual(values = c("All 60: time" = "#009E73", "HDR/LDR: group x time" = "#CC79A7")) +
  scale_shape_manual(values = c("All 60: time" = 16, "HDR/LDR: group x time" = 17)) +
  labs(title = "B  GEE omnibus tests",
       subtitle = "Participant-clustered; values >30 capped\nHDR/LDR interaction is selection-coupled",
       x = expression(-log[10](italic(P))), y = NULL) + base_theme

dist_plot <- dist |>
  filter(dataset_branch == "rarefied_primary",
         comparison %in% c("original_HDR_vs_LDR", "all60_within_subject_phase")) |>
  mutate(
    stage_label = case_when(
      comparison == "all60_within_subject_phase" ~ "All 60: baseline vs endpoint",
      stage == "baseline" ~ "Reference groups: baseline",
      stage == "antibiotic_end" ~ "Reference groups: antibiotic end",
      stage == "endpoint" ~ "Reference groups: endpoint"
    ),
    dispersion = if_else(PERMDISP_p < 0.05, "PERMDISP P < 0.05", "PERMDISP P >= 0.05"),
    distance = factor(distance, levels = c("Bray-Curtis", "binary Jaccard"))
  )
pC <- ggplot(dist_plot, aes(PERMANOVA_R2, reorder(stage_label, PERMANOVA_R2), colour = distance, shape = dispersion)) +
  geom_point(size = 2.7, stroke = 0.7) +
  scale_colour_manual(values = c("Bray-Curtis" = "#0072B2", "binary Jaccard" = "#D55E00")) +
  scale_shape_manual(values = c("PERMDISP P < 0.05" = 17, "PERMDISP P >= 0.05" = 16)) +
  scale_x_continuous(labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "C  Restricted PERMANOVA paired with PERMDISP",
       subtitle = "Triangle marks comparisons where dispersion also differs",
       x = expression(PERMANOVA~R^2), y = NULL) + base_theme

lefse_long <- bind_rows(
  lefse |> transmute(family, correction_scope, metric = "q < 0.05", count = n_q_lt_0_05),
  lefse |> transmute(family, correction_scope, metric = "q < 0.05 and |LDA| >= 2", count = n_q_lt_0_05_and_abs_LDA_ge_2)
) |>
  mutate(family_label = recode(family,
    eligible_268_species_within_each_time = "268 eligible x 16; within-time BH",
    eligible_268_species_global_all_times = "268 eligible x 16; global BH",
    displayed_45_species_within_each_time = "45 displayed x 16; within-time BH",
    displayed_45_species_global_all_times = "45 displayed x 16; global BH"
  ))
phase_labels <- c(
  baseline_HDR_minus_LDR = "Baseline HDR-LDR",
  antibiotic_change_HDR_minus_LDR = "Antibiotic change HDR-LDR",
  late_recovery_change_HDR_minus_LDR = "Late recovery change HDR-LDR"
)
phase_plot <- phase_clr |>
  group_by(contrast) |>
  slice_min(order_by = q_global, n = 5, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    contrast_label = recode(contrast, !!!phase_labels),
    species_label = if_else(nchar(species) > 31, paste0(substr(species, 1, 29), "..."), species),
    label = species_label,
    sig = q_global < 0.05
  )
pD <- ggplot(phase_plot, aes(HDR_minus_LDR, reorder(label, HDR_minus_LDR), colour = contrast_label)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 3, linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.45) +
  geom_point(aes(shape = sig), size = 2.2, stroke = 0.5) +
  scale_colour_manual(values = c(
    "Baseline HDR-LDR" = "#7B3294", "Antibiotic change HDR-LDR" = "#E66101",
    "Late recovery change HDR-LDR" = "#008837"
  )) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("Global q >= 0.05", "Global q < 0.05")) +
  labs(title = "D  Phase CLR contrasts",
       subtitle = "Top 5/contrast; global BH across 804 tests",
       x = "Mean CLR contrast (HDR minus LDR)", y = NULL) + base_theme +
  theme(axis.text.y = element_text(size = 6.3), legend.position = "bottom")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = "R1.6 | Dependence-aware longitudinal and distance inference",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, margin = margin(b = 5))))

pdf_file <- file.path(out, "figures", "R1_6_response_figure.pdf")
tif_file <- file.path(out, "figures", "R1_6_response_figure.tif")
png_file <- file.path(out, "figures", "R1_6_response_figure_preview.png")
ggsave(pdf_file, fig, width = 210, height = 175, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(tif_file, fig, width = 210, height = 175, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", background = "white")
ggsave(png_file, fig, width = 210, height = 175, units = "mm", dpi = 180,
       device = ragg::agg_png, background = "white")

qa <- tibble(
  check = c("60 exposed participants in primary branch", "16 canonical time points", "GEE exchangeable correlation recorded",
            "PERMANOVA uses 9999 restricted permutations", "PERMANOVA paired with PERMDISP", "LEfSe audit has four families",
            "Phase CLR uses 804-test global BH family", "Phase CLR uses 15 participants per group",
            "PDF exists", "TIFF exists", "Preview exists"),
  passed = c(
    n_distinct(long$subject[long$dataset_branch == "rarefied_primary"]) == 60,
    n_distinct(long$canon[long$dataset_branch == "rarefied_primary"]) == 16,
    all(gee_coef$working_correlation == "exchangeable"),
    all(dist$permutations == 9999),
    all(!is.na(dist$PERMDISP_p)), nrow(lefse) == 4,
    all(phase_clr_summary$global_family_size == 804), all(phase_clr_summary$n_HDR == 15 & phase_clr_summary$n_LDR == 15),
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "R1_6_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R1_6_R_sessionInfo.txt"))
if (!all(qa$passed)) stop("R1.6 QA failed")
cat(sprintf("R1.6 complete: %d/%d checks passed\n", sum(qa$passed), nrow(qa)))
