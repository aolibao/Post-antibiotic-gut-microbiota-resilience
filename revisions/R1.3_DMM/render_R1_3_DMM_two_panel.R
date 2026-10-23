#!/home/qqt/anaconda3/bin/Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(ragg)
  library(systemfonts)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
table_dir <- file.path(base_dir, "tables")
figure_dir <- file.path(base_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# The host does not contain Microsoft's proprietary Arial files. Register the
# metric-compatible Liberation Sans files under the requested Arial family for
# ragg devices; cairo_pdf resolves the same Arial request through fontconfig.
font_dir <- "/usr/share/fonts/liberation-sans"
font_files <- c(
  plain = file.path(font_dir, "LiberationSans-Regular.ttf"),
  bold = file.path(font_dir, "LiberationSans-Bold.ttf"),
  italic = file.path(font_dir, "LiberationSans-Italic.ttf"),
  bolditalic = file.path(font_dir, "LiberationSans-BoldItalic.ttf")
)
stopifnot(all(file.exists(font_files)))
systemfonts::register_font(
  name = "Arial",
  plain = font_files[["plain"]],
  bold = font_files[["bold"]],
  italic = font_files[["italic"]],
  bolditalic = font_files[["bolditalic"]]
)

cross <- read.csv(file.path(table_dir, "R1_3_DMM_HDR_LDR_crossclassification.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
assoc <- read.csv(file.path(table_dir, "R1_3_DMM_HDR_LDR_association_test.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
models <- read.csv(file.path(table_dir, "R1_3_DMM_conditional_models_HC3.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(nrow(assoc) == 1L)
cross$group <- factor(cross$Var1, levels = c("LDR", "HDR"))
cross$type <- factor(cross$Var2,
                     levels = c("DMM1:Bacteroides", "DMM2:Segatella"))
cross <- cross[order(cross$group, cross$type), ]
cross$proportion <- ave(cross$Freq, cross$group, FUN = function(x) x / sum(x))

joint <- models[models$model_id == "Shannon + DMM" & models$term != "(Intercept)", ]
stopifnot(nrow(joint) == 2L)
joint$parameter <- ifelse(joint$term == "z_shannon", "Baseline Shannon", "DMM type")
joint$y <- ifelse(joint$parameter == "Baseline Shannon", 2, 1)
joint$shape <- ifelse(joint$parameter == "Baseline Shannon", 21, 22)

source_a <- data.frame(
  panel = "A",
  left_to_right_group = as.character(cross$group),
  fill_code = as.character(cross$type),
  count = cross$Freq,
  proportion = cross$proportion,
  test = assoc$test,
  p_value = assoc$p_value,
  stringsAsFactors = FALSE
)
source_b <- data.frame(
  panel = "B",
  top_to_bottom_parameter = joint$parameter,
  estimate = joint$estimate,
  CI95_low = joint$CI95_low,
  CI95_high = joint$CI95_high,
  p_value = joint$p_value,
  inference = "HC3",
  stringsAsFactors = FALSE
)
write.csv(rbind(
  transform(source_a, top_to_bottom_parameter = NA_character_, estimate = NA_real_,
            CI95_low = NA_real_, CI95_high = NA_real_, inference = NA_character_),
  transform(source_b, left_to_right_group = NA_character_, fill_code = NA_character_,
            count = NA_real_, proportion = NA_real_, test = NA_character_)
), file.path(table_dir, "R1_3_DMM_two_panel_source_data.csv"), row.names = FALSE)

fmt_p <- function(x) paste0("P = ", formatC(x, digits = 3, format = "f"))

# Manuscript-wide HDR/LDR palette, sampled from the restored main figures.
ldr_blue <- "#0072B5"
ldr_blue_light <- "#B9DCEC"
hdr_red <- "#BC3C29"
hdr_red_light <- "#E9BCB5"
zero_grey <- "#A9A9A9"

base_theme <- theme_void(base_family = "Arial", base_size = 8.2) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    plot.margin = margin(5, 6, 5, 6, unit = "pt"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Panel A: each dot is one participant.  The row colour preserves the
# manuscript's recovery-group vocabulary; fill state carries DMM type.
dot_a <- do.call(rbind, lapply(c("LDR", "HDR"), function(g) {
  n_dmm1 <- cross$Freq[cross$group == g & cross$type == "DMM1:Bacteroides"]
  n_total <- sum(cross$Freq[cross$group == g])
  data.frame(
    group = g,
    x = seq_len(n_total),
    y = if (g == "LDR") 1.55 else 0.75,
    dmm = c(rep("DMM1", n_dmm1), rep("DMM2", n_total - n_dmm1)),
    outline = if (g == "LDR") ldr_blue else hdr_red,
    fill = c(rep(if (g == "LDR") ldr_blue else hdr_red, n_dmm1),
             rep("white", n_total - n_dmm1)),
    stringsAsFactors = FALSE
  )
}))

p_a <- ggplot(dot_a, aes(x, y)) +
  geom_point(aes(colour = outline, fill = fill), shape = 21, size = 2.55,
             stroke = 0.58) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0.25, 15.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.32, 2.28), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  annotate("text", x = 0.28, y = 2.225, label = "A", family = "Arial",
           fontface = "bold", colour = "black", size = 3.18, hjust = 0) +
  annotate("text", x = 8, y = 2.03,
           label = paste0("italic(P) == ", formatC(assoc$p_value, digits = 3, format = "f")),
           parse = TRUE, family = "Arial", colour = "black", size = 2.66, hjust = 0.5) +
  base_theme +
  theme(plot.margin = margin(5, 11, 5, 5, unit = "pt"))

# Panel B: top=covariate-adjusted baseline Shannon; bottom=DMM2 versus DMM1.
# Circle and square provide redundant row encoding; manuscript red/blue hues
# improve visual continuity while the caption carries all decoding.
x_limits <- c(-1.42, 1.20)
joint$plot_colour <- ifelse(joint$parameter == "Baseline Shannon", ldr_blue, hdr_red)
joint$y <- ifelse(joint$parameter == "Baseline Shannon", 1.55, 0.75)
joint$p_plot <- paste0("italic(P) == ", formatC(joint$p_value, digits = 3, format = "f"))
p_b <- ggplot(joint, aes(estimate, y)) +
  geom_segment(aes(x = 0, xend = 0, y = 0.34, yend = 1.92),
               inherit.aes = FALSE, linetype = "22", linewidth = 0.42,
               colour = zero_grey) +
  geom_segment(aes(x = CI95_low, xend = CI95_high, yend = y),
               linewidth = 0.58, colour = joint$plot_colour, lineend = "round") +
  geom_segment(aes(x = CI95_low, xend = CI95_low, y = y - 0.10, yend = y + 0.10),
               linewidth = 0.52, colour = joint$plot_colour, lineend = "round") +
  geom_segment(aes(x = CI95_high, xend = CI95_high, y = y - 0.10, yend = y + 0.10),
               linewidth = 0.52, colour = joint$plot_colour, lineend = "round") +
  geom_point(aes(shape = factor(shape)), size = 3.05, stroke = 0.48,
             colour = "white", fill = joint$plot_colour, show.legend = FALSE) +
  scale_shape_manual(values = c("21" = 21, "22" = 22), guide = "none") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.32, 2.28), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  annotate("text", x = x_limits[1], y = 2.225, label = "B", family = "Arial",
           fontface = "bold", colour = "black", size = 3.18, hjust = 0) +
  geom_text(aes(x = 1.18, label = p_plot), parse = TRUE, family = "Arial",
            colour = "black", size = 2.66, hjust = 1) +
  base_theme +
  theme(plot.margin = margin(5, 5, 5, 11, unit = "pt"))

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.90, 1.35), "null"))))
  print(p_a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), newpage = FALSE)
  print(p_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2), newpage = FALSE)
}

out_base <- file.path(figure_dir, "R1_3_DMM_two_panel_supplement")
width_mm <- 178
height_mm <- 68

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = width_mm / 25.4,
                     height = height_mm / 25.4, family = "Arial",
                     bg = "white", onefile = TRUE)
draw_figure()
dev.off()

ragg::agg_tiff(paste0(out_base, ".tif"), width = width_mm, height = height_mm,
               units = "mm", res = 600, compression = "lzw",
               background = "white", scaling = 1)
draw_figure()
dev.off()

ragg::agg_png(paste0(out_base, "_word.png"), width = width_mm, height = height_mm,
              units = "mm", res = 240, background = "white", scaling = 1)
draw_figure()
dev.off()

writeLines(c(
  "R-only two-panel rendering complete.",
  paste0("Panel A Fisher P: ", format(assoc$p_value, digits = 8)),
  paste0("Panel B baseline Shannon HC3 P: ", format(joint$p_value[joint$parameter == "Baseline Shannon"], digits = 8)),
  paste0("Panel B DMM type HC3 P: ", format(joint$p_value[joint$parameter == "DMM type"], digits = 8)),
  "Figure text restricted to A, B and P values.",
  paste0("Manuscript palette: LDR blue ", ldr_blue, "; HDR red ", hdr_red, "."),
  "Requested family: Arial; host-compatible registered file: Liberation Sans.",
  paste0("Final size: ", width_mm, " x ", height_mm, " mm; TIFF 600 dpi.")
), file.path(base_dir, "logs", "R1_3_DMM_two_panel_render_summary.txt"))

message("Wrote: ", out_base, ".pdf/.tif/_word.png")
