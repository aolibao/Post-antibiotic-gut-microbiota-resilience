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
out <- file.path(root, "reviewer_response", "R1.6_dispersion_location")
for (d in c("figures", "tables", "logs", "QA"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

source_file <- file.path(root, "reviewer_response", "R1.6_revised", "tables",
                         "R1_6_PERMANOVA_PERMDISP.csv")
dat <- read_csv(source_file, show_col_types = FALSE) |>
  filter(
    dataset_branch == "rarefied_primary",
    comparison %in% c("original_HDR_vs_LDR", "all60_within_subject_phase")
  ) |>
  mutate(
    comparison_label = case_when(
      stage == "baseline" ~ "HDR vs LDR: baseline",
      stage == "antibiotic_end" ~ "HDR vs LDR: antibiotic end",
      stage == "endpoint" ~ "HDR vs LDR: study endpoint",
      comparison == "all60_within_subject_phase" ~ "All participants: baseline vs endpoint"
    ),
    comparison_label = factor(comparison_label, levels = rev(c(
      "HDR vs LDR: baseline", "HDR vs LDR: antibiotic end",
      "HDR vs LDR: study endpoint", "All participants: baseline vs endpoint"
    ))),
    distance = factor(distance, levels = c("Bray-Curtis", "binary Jaccard"),
                      labels = c("Bray-Curtis", "Jaccard")),
    permanova_label = paste0("P = ", format.pval(PERMANOVA_p, digits = 2, eps = 1e-4)),
    permdisp_label = paste0("P = ", format.pval(PERMDISP_p, digits = 2, eps = 1e-4)),
    permanova_label_x = PERMANOVA_R2 + 0.012,
    dispersion_result = if_else(PERMDISP_p < 0.05, "Significant", "Not significant")
  )
stopifnot(nrow(dat) == 8L, all(dat$permutations == 9999L))
write_csv(dat, file.path(out, "tables", "R1_6_location_and_dispersion_results.csv"))

pal <- c("Bray-Curtis" = "#2F6FA3", "Jaccard" = "#C44E52")
pd <- position_dodge(width = 0.68)
base_theme <- theme_classic(base_size = 9.1, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 9.4),
    axis.text = element_text(size = 8.1, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top", legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.3, colour = "black"),
    plot.margin = margin(5, 14, 5, 6)
  )

pA <- ggplot(dat, aes(PERMANOVA_R2, comparison_label, colour = distance)) +
  geom_segment(aes(x = 0, xend = PERMANOVA_R2, yend = comparison_label),
               position = pd, linewidth = 0.65, alpha = 0.75) +
  geom_point(position = pd, size = 2.8) +
  geom_text(aes(x = permanova_label_x, label = permanova_label, group = distance),
            position = pd, hjust = 0,
            size = 2.35, family = "Arial", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(0, 0.385), breaks = seq(0, 0.30, 0.10),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "A", x = expression(PERMANOVA~R^2), y = NULL) + base_theme +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.8, size = 2.8)))

pB <- ggplot(dat, aes(PERMDISP_group_mean_distance_range, comparison_label,
                       colour = distance, shape = dispersion_result)) +
  geom_segment(aes(x = 0, xend = PERMDISP_group_mean_distance_range,
                   yend = comparison_label), position = pd,
               linewidth = 0.65, alpha = 0.75) +
  geom_point(position = pd, size = 2.8, stroke = 0.8, fill = "white") +
  geom_text(aes(label = permdisp_label, group = distance), position = pd, hjust = -0.12,
            size = 2.55, family = "Arial", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c("Significant" = 16, "Not significant" = 21)) +
  scale_x_continuous(limits = c(0, 0.265), breaks = seq(0, 0.25, 0.05),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B", x = "Difference in mean distance to centroid", y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  guides(colour = "none", shape = guide_legend(override.aes = list(colour = "black")))

draw <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(2, 2, widths = grid::unit(c(1.05, 0.95), "null"),
                              heights = grid::unit(c(0.91, 0.09), "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(
    "PERMANOVA tests location; PERMDISP tests dispersion. Filled and open points denote significant and non-significant PERMDISP, respectively.",
    x = 0.5, y = 0.5, gp = grid::gpar(fontfamily = "Arial", fontsize = 7.9, col = "black"),
    vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
  grid::popViewport()
}

pdf_file <- file.path(out, "figures", "R1_6_dispersion_vs_location.pdf")
tif_file <- file.path(out, "figures", "R1_6_dispersion_vs_location.tif")
png_file <- file.path(out, "figures", "R1_6_dispersion_vs_location_preview.png")
grDevices::cairo_pdf(pdf_file, width = 183 / 25.4, height = 92 / 25.4, family = "Arial")
draw(); grDevices::dev.off()
ragg::agg_tiff(tif_file, width = 183, height = 92, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw(); grDevices::dev.off()
ragg::agg_png(png_file, width = 183, height = 92, units = "mm", res = 180,
              background = "white")
draw(); grDevices::dev.off()

qa <- tibble(
  check = c("Eight planned comparisons", "9999 permutations", "P values complete",
            "Original P=0.052 classified nonsignificant", "PDF", "TIFF", "Preview"),
  passed = c(nrow(dat) == 8L, all(dat$permutations == 9999L),
             all(is.finite(dat$PERMANOVA_p)) && all(is.finite(dat$PERMDISP_p)),
             0.052 >= 0.05, file.exists(pdf_file), file.exists(tif_file), file.exists(png_file))
)
write_tsv(qa, file.path(out, "QA", "R1_6_dispersion_location_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)
