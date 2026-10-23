#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ragg)
  library(systemfonts)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.6_multiplicity")
for (d in c("figures", "tables", "logs", "QA"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

audit_file <- file.path(root, "reviewer_response", "03_tables", "lefse_multiplicity_audit.csv")
result_file <- file.path(root, "reviewer_response", "R1.6", "tables",
                         "R1_6_participant_phase_CLR_global_FDR.csv")
summary_file <- file.path(root, "reviewer_response", "R1.6", "tables",
                          "R1_6_participant_phase_CLR_summary.csv")
stopifnot(file.exists(audit_file), file.exists(result_file), file.exists(summary_file))

audit <- read_csv(audit_file, show_col_types = FALSE)
results <- read_csv(result_file, show_col_types = FALSE)
summary <- read_csv(summary_file, show_col_types = FALSE)

phase_names <- c(
  baseline_HDR_minus_LDR = "Baseline",
  antibiotic_change_HDR_minus_LDR = "Antibiotic change",
  late_recovery_change_HDR_minus_LDR = "Late recovery change"
)

count_data <- summary |>
  transmute(
    phase = factor(unname(phase_names[contrast]), levels = unname(phase_names)),
    `Nominal P < 0.05` = nominal_p_lt_0_05,
    `Global q < 0.05` = global_q_lt_0_05
  ) |>
  tidyr::pivot_longer(-phase, names_to = "criterion", values_to = "count") |>
  mutate(criterion = factor(criterion, levels = c("Nominal P < 0.05", "Global q < 0.05")))

top <- results |>
  mutate(phase = factor(unname(phase_names[contrast]), levels = unname(phase_names))) |>
  group_by(phase) |>
  slice_min(q_global, n = 3, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    direction = if_else(HDR_minus_LDR >= 0, "Higher in HDR", "Higher in LDR"),
    fdr = if_else(q_global < 0.05, "Global q < 0.05", "Global q >= 0.05"),
    q_label = paste0("q = ", format.pval(q_global, digits = 2, eps = 1e-3)),
    species_key = paste(species, phase, sep = "___")
  ) |>
  group_by(phase) |>
  arrange(HDR_minus_LDR, .by_group = TRUE) |>
  mutate(species_key = factor(species_key, levels = unique(species_key))) |>
  ungroup()

write_csv(audit, file.path(out, "tables", "R1_6_LEfSe_multiplicity_audit.csv"))
write_csv(summary, file.path(out, "tables", "R1_6_phase_CLR_global_FDR_summary.csv"))
write_csv(top, file.path(out, "tables", "R1_6_phase_CLR_top_effects.csv"))

base_theme <- theme_classic(base_size = 9.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8.3, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 8.2, colour = "black"),
    plot.margin = margin(5, 8, 5, 5)
  )

pA <- ggplot(count_data, aes(phase, count, fill = criterion)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.58, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = count), position = position_dodge(width = 0.68), vjust = -0.35,
            family = "Arial", size = 2.8, colour = "black") +
  scale_fill_manual(values = c("Nominal P < 0.05" = "#B7B7B7", "Global q < 0.05" = "#2F6FA3")) +
  scale_y_continuous(limits = c(0, 91), breaks = seq(0, 80, 20), expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = NULL, y = "Number of species") +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

direction_pal <- c("Higher in HDR" = "#C44E52", "Higher in LDR" = "#2F6FA3")
pB <- ggplot(top, aes(HDR_minus_LDR, species_key, colour = direction, shape = fdr)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.45, colour = "black") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.55) +
  geom_point(aes(fill = direction), size = 2.7, stroke = 0.8) +
  geom_text(aes(x = if_else(HDR_minus_LDR >= 0, ci_high + 0.10, ci_low - 0.10),
                label = q_label, hjust = if_else(HDR_minus_LDR >= 0, 0, 1)),
            colour = "black", size = 2.45, family = "Arial", show.legend = FALSE) +
  facet_grid(phase ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = direction_pal) +
  scale_fill_manual(values = direction_pal) +
  scale_shape_manual(values = c("Global q < 0.05" = 16, "Global q >= 0.05" = 21)) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_x_continuous(limits = c(-4.8, 4.2), breaks = seq(-4, 4, 2), expand = expansion(mult = c(0, 0))) +
  labs(title = "B", x = "Participant-level CLR contrast (HDR minus LDR)", y = NULL) +
  base_theme +
  theme(
    strip.background = element_rect(fill = "#F2F2F2", colour = "black", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8.2, colour = "black"),
    axis.text.y = element_text(size = 7.2, face = "italic"),
    legend.position = "top"
  ) +
  guides(fill = "none", colour = guide_legend(order = 1), shape = "none")

draw <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 2, widths = grid::unit(c(0.40, 0.60), "null"),
                           heights = grid::unit(c(0.92, 0.08), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(
    "268 species × 3 prespecified contrasts = 804 tests in one global BH family (15 HDR, 15 LDR); filled points denote global q < 0.05.",
    x = 0.5, y = 0.5, gp = grid::gpar(fontfamily = "Arial", fontsize = 8.0, col = "black"),
    vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
  grid::popViewport()
}

pdf_file <- file.path(out, "figures", "R1_6_multiplicity_phase_CLR.pdf")
tif_file <- file.path(out, "figures", "R1_6_multiplicity_phase_CLR.tif")
png_file <- file.path(out, "figures", "R1_6_multiplicity_phase_CLR_preview.png")
grDevices::cairo_pdf(pdf_file, width = 183 / 25.4, height = 125 / 25.4, family = "Arial")
draw(); grDevices::dev.off()
ragg::agg_tiff(tif_file, width = 183, height = 125, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw(); grDevices::dev.off()
ragg::agg_png(png_file, width = 183, height = 125, units = "mm", res = 180,
              background = "white")
draw(); grDevices::dev.off()

qa <- tibble(
  check = c("Original eligible family has 4288 tests", "Replacement global family has 804 tests",
            "Three prespecified contrasts", "One globally significant feature",
            "30 independent participants", "PDF", "TIFF", "Preview"),
  passed = c(
    audit$n_tests[audit$family == "eligible_268_species_global_all_times"] == 4288,
    unique(summary$global_family_size) == 804, nrow(summary) == 3,
    sum(results$q_global < 0.05) == 1,
    all(summary$n_HDR == 15 & summary$n_LDR == 15),
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "R1_6_multiplicity_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)
