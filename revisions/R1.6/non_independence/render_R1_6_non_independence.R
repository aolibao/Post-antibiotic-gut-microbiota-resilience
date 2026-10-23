#!/home/qqt/anaconda3/bin/Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(grid)
  library(ragg)
  library(systemfonts)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
figure_dir <- file.path(base_dir, "figures")
table_dir <- file.path(base_dir, "tables")
log_dir <- file.path(base_dir, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# The host lacks Microsoft's proprietary Arial files. Register the metrically
# compatible Liberation Sans family under Arial for ragg output.
font_dir <- "/usr/share/fonts/liberation-sans"
font_files <- c(
  plain = file.path(font_dir, "LiberationSans-Regular.ttf"),
  bold = file.path(font_dir, "LiberationSans-Bold.ttf"),
  italic = file.path(font_dir, "LiberationSans-Italic.ttf"),
  bolditalic = file.path(font_dir, "LiberationSans-BoldItalic.ttf")
)
stopifnot(all(file.exists(font_files)))
systemfonts::register_font(
  name = "Arial", plain = font_files[["plain"]], bold = font_files[["bold"]],
  italic = font_files[["italic"]], bolditalic = font_files[["bolditalic"]]
)

time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
group_cols <- c(LDR = "#0072B5", HDR = "#BC3C29")
individual_file <- file.path(
  project_root, "reviewer_response", "R1.1", "full_original_layout",
  "D1_endpoint_abs_shannon", "tables", "sample_metrics_with_alternative_groups.csv"
)
summary_file <- file.path(
  project_root, "reviewer_response", "R1.6_revised", "tables", "R1_6_trajectory_summary.csv"
)
omnibus_file <- file.path(
  project_root, "reviewer_response", "R1.6_revised", "tables", "R1_6_repeated_measure_omnibus.csv"
)
stopifnot(file.exists(individual_file), file.exists(summary_file), file.exists(omnibus_file))

individual <- read_csv(individual_file, show_col_types = FALSE) |>
  filter(recovery_group %in% c("LDR", "HDR"), canon %in% time_order) |>
  mutate(
    recovery_group = factor(recovery_group, levels = c("LDR", "HDR")),
    canon = factor(canon, levels = time_order),
    time_index = match(as.character(canon), time_order)
  )
stopifnot(n_distinct(individual$subject) == 30L, nrow(individual) == 480L)
stopifnot(all(table(individual$recovery_group) == 240L))

individual_long <- bind_rows(
  transmute(individual, subject, recovery_group, canon, time_index,
            outcome = "Shannon diversity", value = shannon),
  transmute(individual, subject, recovery_group, canon, time_index,
            outcome = "Bray-Curtis distance", value = bray_to_baseline_centroid),
  transmute(individual, subject, recovery_group, canon, time_index,
            outcome = "Jaccard distance", value = mean_jaccard_to_baseline)
)

summary_check <- individual_long |>
  group_by(recovery_group, canon, time_index, outcome) |>
  summarise(
    n_subjects = n_distinct(subject), mean = mean(value),
    se = sd(value) / sqrt(n()), .groups = "drop"
  )
stopifnot(all(summary_check$n_subjects == 15L))

# Verify that the independently reconstructed summaries match the established
# audited R1.6 trajectory table before drawing.
audited_summary <- read_csv(summary_file, show_col_types = FALSE) |>
  mutate(canon = factor(canon, levels = time_order), time_index = match(as.character(canon), time_order)) |>
  select(recovery_group, canon, time_index,
         shannon_mean, shannon_se, bray_mean, bray_se, jaccard_mean, jaccard_se)
check_wide <- summary_check |>
  tidyr::pivot_wider(names_from = outcome, values_from = c(mean, se)) |>
  mutate(recovery_group = as.character(recovery_group))
check_join <- left_join(audited_summary, check_wide,
                        by = c("recovery_group", "canon", "time_index"))
stopifnot(
  max(abs(check_join$shannon_mean - check_join$`mean_Shannon diversity`)) < 1e-10,
  max(abs(check_join$shannon_se - check_join$`se_Shannon diversity`)) < 1e-10,
  max(abs(check_join$bray_mean - check_join$`mean_Bray-Curtis distance`)) < 1e-10,
  max(abs(check_join$bray_se - check_join$`se_Bray-Curtis distance`)) < 1e-10,
  max(abs(check_join$jaccard_mean - check_join$`mean_Jaccard distance`)) < 1e-10,
  max(abs(check_join$jaccard_se - check_join$`se_Jaccard distance`)) < 1e-10
)

omnibus <- read_csv(omnibus_file, show_col_types = FALSE) |>
  filter(
    dataset_branch == "rarefied_primary", model == "original_HDR_LDR",
    outcome %in% c("shannon", "bray_to_personal_baseline_centroid",
                   "jaccard_to_personal_baseline_consensus")
  ) |>
  mutate(
    outcome_label = recode(
      outcome,
      shannon = "Shannon diversity",
      bray_to_personal_baseline_centroid = "Bray-Curtis distance",
      jaccard_to_personal_baseline_consensus = "Jaccard distance"
    )
  )
stopifnot(nrow(omnibus) == 3L, all(omnibus$n_subjects == 30L),
          all(omnibus$n_observations == 480L), all(omnibus$df == 15L))

write_csv(individual_long, file.path(table_dir, "R1_6_non_independence_source_values.csv"))
write_csv(summary_check, file.path(table_dir, "R1_6_non_independence_trajectory_summary.csv"))
write_csv(omnibus, file.path(table_dir, "R1_6_non_independence_GEE_tests.csv"))

fmt_p <- function(x) {
  exponent <- floor(log10(x))
  mantissa <- x / (10^exponent)
  superscript <- function(value) {
    chars <- strsplit(as.character(value), "", fixed = TRUE)[[1]]
    map <- c("-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³",
             "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
    paste0(unname(map[chars]), collapse = "")
  }
  paste0(formatC(mantissa, digits = 2, format = "f"), " × 10", superscript(exponent))
}

panel_spec <- data.frame(
  outcome = c("Shannon diversity", "Bray-Curtis distance", "Jaccard distance"),
  tag = c("A", "B", "C"),
  y_label = c("Shannon diversity", "Bray-Curtis distance\nto personal baseline",
              "Jaccard distance\nto personal baseline"),
  stringsAsFactors = FALSE
)

theme_panel <- theme_classic(base_family = "Arial", base_size = 7.2) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks.length = unit(1.5, "pt"),
    axis.text = element_text(size = 5.6, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 6.6, margin = margin(r = 4)),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0,
                              margin = margin(b = 1.5)),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(4, 5, 3, 4, unit = "pt")
  )

make_panel <- function(outcome_name, tag, y_label) {
  raw <- individual_long |> filter(outcome == outcome_name)
  sm <- summary_check |> filter(outcome == outcome_name) |>
    mutate(ci_low = mean - 1.96 * se, ci_high = mean + 1.96 * se)
  test <- omnibus |> filter(outcome_label == outcome_name)
  stopifnot(nrow(test) == 1L)
  p_label <- paste0("Group × time: P = ", fmt_p(test$p_value))

  ggplot() +
    annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#F1F7FB", colour = NA) +
    annotate("rect", xmin = 4.5, xmax = 10.5, ymin = -Inf, ymax = Inf,
             fill = "#FCF1EE", colour = NA) +
    annotate("rect", xmin = 10.5, xmax = 16.5, ymin = -Inf, ymax = Inf,
             fill = "#F1F8F3", colour = NA) +
    geom_vline(xintercept = c(4.5, 10.5), colour = "#C9C9C9",
               linetype = "22", linewidth = 0.30) +
    geom_line(data = raw,
              aes(time_index, value, group = subject, colour = recovery_group),
              linewidth = 0.23, alpha = 0.13) +
    geom_ribbon(data = sm,
                aes(time_index, ymin = ci_low, ymax = ci_high,
                    group = recovery_group, fill = recovery_group),
                inherit.aes = FALSE, alpha = 0.14, colour = NA) +
    geom_line(data = sm,
              aes(time_index, mean, group = recovery_group, colour = recovery_group),
              linewidth = 0.72, lineend = "round") +
    geom_point(data = sm,
               aes(time_index, mean, shape = recovery_group,
                   fill = recovery_group, colour = recovery_group),
               size = 1.55, stroke = 0.30) +
    scale_colour_manual(values = group_cols, breaks = c("LDR", "HDR")) +
    scale_fill_manual(values = group_cols, breaks = c("LDR", "HDR")) +
    scale_shape_manual(values = c(LDR = 21, HDR = 22)) +
    scale_x_continuous(breaks = seq_along(time_order), labels = time_order,
                       expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.17))) +
    coord_cartesian(xlim = c(0.7, 16.3), clip = "on") +
    annotate("text", x = 8.5, y = Inf, label = p_label,
             family = "Arial", colour = "black", size = 2.35,
             hjust = 0.5, vjust = 1.25) +
    labs(title = tag, y = y_label) +
    theme_panel
}

panels <- lapply(seq_len(nrow(panel_spec)), function(i) {
  make_panel(panel_spec$outcome[i], panel_spec$tag[i], panel_spec$y_label[i])
})

draw_legend <- function() {
  grid.lines(x = unit(c(0.405, 0.445), "npc"), y = unit(c(0.5, 0.5), "npc"),
             gp = gpar(col = group_cols[["LDR"]], lwd = 1.7, lineend = "round"))
  grid.points(x = unit(0.425, "npc"), y = unit(0.5, "npc"), pch = 21,
              size = unit(2.3, "mm"), gp = gpar(col = group_cols[["LDR"]], fill = group_cols[["LDR"]], lwd = 0.5))
  grid.text("LDR", x = unit(0.466, "npc"), y = unit(0.5, "npc"),
            just = "left", gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
  grid.lines(x = unit(c(0.535, 0.575), "npc"), y = unit(c(0.5, 0.5), "npc"),
             gp = gpar(col = group_cols[["HDR"]], lwd = 1.7, lineend = "round"))
  grid.points(x = unit(0.555, "npc"), y = unit(0.5, "npc"), pch = 22,
              size = unit(2.3, "mm"), gp = gpar(col = group_cols[["HDR"]], fill = group_cols[["HDR"]], lwd = 0.5))
  grid.text("HDR", x = unit(0.596, "npc"), y = unit(0.5, "npc"),
            just = "left", gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2, 3, heights = unit(c(0.10, 0.90), "null"), widths = unit(c(1, 1, 1), "null")
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  draw_legend()
  popViewport()
  for (i in seq_along(panels)) {
    print(panels[[i]], vp = viewport(layout.pos.row = 2, layout.pos.col = i), newpage = FALSE)
  }
  popViewport()
}

out_base <- file.path(figure_dir, "R1_6_non_independence_GEE_trajectories")
width_mm <- 183
height_mm <- 78

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

ragg::agg_png(paste0(out_base, "_preview.png"), width = width_mm, height = height_mm,
              units = "mm", res = 240, background = "white", scaling = 1)
draw_figure()
dev.off()

writeLines(c(
  "R1.6 non-independence response figure rendered entirely in R.",
  "Model: participant-clustered GEE; exchangeable working correlation; categorical time; group-by-time interaction; robust inference.",
  "Analytic sample: 30 participants (15 LDR, 15 HDR), 16 visits, 480 observations per outcome.",
  paste0("Shannon: Wald chi-square(15) = ", formatC(omnibus$wald_chi2[omnibus$outcome_label == "Shannon diversity"], digits = 2, format = "f"),
         ", P = ", format(omnibus$p_value[omnibus$outcome_label == "Shannon diversity"], scientific = TRUE, digits = 4)),
  paste0("Bray-Curtis: Wald chi-square(15) = ", formatC(omnibus$wald_chi2[omnibus$outcome_label == "Bray-Curtis distance"], digits = 2, format = "f"),
         ", P = ", format(omnibus$p_value[omnibus$outcome_label == "Bray-Curtis distance"], scientific = TRUE, digits = 4)),
  paste0("Jaccard: Wald chi-square(15) = ", formatC(omnibus$wald_chi2[omnibus$outcome_label == "Jaccard distance"], digits = 2, format = "f"),
         ", P = ", format(omnibus$p_value[omnibus$outcome_label == "Jaccard distance"], scientific = TRUE, digits = 4)),
  "Thin lines: participant trajectories. Thick lines and ribbons: group means and descriptive 95% confidence intervals.",
  paste0("Output size: ", width_mm, " x ", height_mm, " mm; TIFF 600 dpi LZW."),
  "Requested family: Arial; host-compatible registered file: Liberation Sans."
), file.path(log_dir, "R1_6_non_independence_render_summary.txt"))

qa <- data.frame(
  check = c(
    "30 participants", "15 participants per group", "16 repeated visits",
    "480 observations per outcome", "audited summaries reproduced",
    "three GEE interaction tests", "PDF exists", "TIFF exists"
  ),
  passed = c(
    n_distinct(individual$subject) == 30L,
    all(table(droplevels(individual$recovery_group)) / length(time_order) == 15L),
    n_distinct(individual$canon) == 16L,
    all(table(individual_long$outcome) == 480L),
    TRUE, nrow(omnibus) == 3L,
    file.exists(paste0(out_base, ".pdf")), file.exists(paste0(out_base, ".tif"))
  )
)
write_tsv(qa, file.path(log_dir, "R1_6_non_independence_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "R_sessionInfo.txt"))
stopifnot(all(qa$passed))

message("Wrote: ", out_base, ".pdf/.tif/_preview.png")
