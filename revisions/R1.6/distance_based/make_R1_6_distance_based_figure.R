#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ragg)
  library(systemfonts)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grep("--file=", args)]
script_file <- normalizePath(sub("--file=", "", script_arg))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.6_distance_based")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "QA"), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file)) {
  register_font(name = "Arial", plain = font_file)
}

source_file <- file.path(root, "reviewer_response", "R1.6_revised", "tables", "R1_6_PERMANOVA_PERMDISP.csv")
stopifnot(file.exists(source_file))
all_results <- read_csv(source_file, show_col_types = FALSE)

plot_data <- all_results |>
  filter(
    dataset_branch == "rarefied_primary",
    comparison == "original_HDR_vs_LDR",
    stage %in% c("baseline", "antibiotic_end", "endpoint")
  ) |>
  mutate(
    stage = factor(stage, levels = c("baseline", "antibiotic_end", "endpoint"),
                   labels = c("Baseline", "Antibiotic end", "Study endpoint")),
    distance = factor(distance, levels = c("Bray-Curtis", "binary Jaccard"),
                      labels = c("Bray-Curtis", "Jaccard")),
    permanova_label = paste0("P = ", format.pval(PERMANOVA_p, digits = 2, eps = 1e-4)),
    permdisp_label = paste0("P = ", format.pval(PERMDISP_p, digits = 2, eps = 1e-4))
  )

stopifnot(nrow(plot_data) == 6L, all(plot_data$permutations == 9999L))
write_csv(plot_data, file.path(out, "tables", "R1_6_distance_based_results.csv"))

pal <- c("Bray-Curtis" = "#2F6FA3", "Jaccard" = "#C44E52")
pd <- position_dodge(width = 0.52)
base_theme <- theme_classic(base_size = 9.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5, colour = "black"),
    legend.key.width = grid::unit(10, "pt"),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 12, 5, 6)
  )

pA <- ggplot(plot_data, aes(stage, PERMANOVA_R2, colour = distance, group = distance)) +
  geom_line(position = pd, linewidth = 0.75) +
  geom_point(position = pd, size = 2.7, stroke = 0.45) +
  geom_text(aes(label = permanova_label), position = pd, vjust = -1.05,
            size = 2.7, family = "Arial", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(limits = c(0, 0.285), breaks = seq(0, 0.25, 0.05),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = NULL, y = expression(PERMANOVA~R^2)) +
  base_theme

pB <- ggplot(plot_data, aes(stage, PERMDISP_group_mean_distance_range,
                            colour = distance, group = distance)) +
  geom_line(position = pd, linewidth = 0.75) +
  geom_point(position = pd, size = 2.7, stroke = 0.45) +
  geom_text(aes(label = permdisp_label), position = pd, vjust = -1.05,
            size = 2.7, family = "Arial", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(limits = c(0, 0.235), breaks = seq(0, 0.20, 0.05),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B", x = NULL, y = "Difference in mean distance\nto group centroid") +
  base_theme + theme(legend.position = "none")

draw_figure <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 2, ncol = 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(0.93, 0.07), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(
    "P values were obtained from 9,999 participant-label permutations; PERMDISP was evaluated in parallel.",
    x = 0.5, y = 0.48, just = "centre",
    gp = grid::gpar(fontfamily = "Arial", fontsize = 8.2, col = "black"),
    vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
  grid::popViewport()
}

pdf_file <- file.path(out, "figures", "R1_6_distance_based_permutation_analysis.pdf")
tif_file <- file.path(out, "figures", "R1_6_distance_based_permutation_analysis.tif")
png_file <- file.path(out, "figures", "R1_6_distance_based_permutation_analysis_preview.png")

grDevices::cairo_pdf(pdf_file, width = 183 / 25.4, height = 79 / 25.4, family = "Arial")
draw_figure()
grDevices::dev.off()

ragg::agg_tiff(tif_file, width = 183, height = 79, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw_figure()
grDevices::dev.off()

ragg::agg_png(png_file, width = 183, height = 79, units = "mm", res = 180,
              background = "white")
draw_figure()
grDevices::dev.off()

qa <- tibble(
  check = c(
    "Six prespecified primary comparisons retained",
    "All tests used 9999 permutations",
    "All PERMANOVA and PERMDISP P values present",
    "PDF exists", "TIFF exists", "Preview exists"
  ),
  passed = c(
    nrow(plot_data) == 6L,
    all(plot_data$permutations == 9999L),
    all(is.finite(plot_data$PERMANOVA_p)) && all(is.finite(plot_data$PERMDISP_p)),
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "R1_6_distance_based_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)
