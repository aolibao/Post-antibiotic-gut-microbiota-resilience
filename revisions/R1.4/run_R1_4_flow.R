#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(base_dir, "tables"); figure_dir <- file.path(base_dir, "figures"); log_dir <- file.path(base_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(root_dir, "raw", "meta.csv"), show_col_types = FALSE)
required <- c("ids", "subject", "abx", "canon")
if (!all(required %in% names(meta))) stop("Missing required metadata columns")
if (anyDuplicated(meta$ids)) stop("Sample IDs are not unique")

per_subject <- meta |>
  group_by(subject) |>
  summarise(samples_n = n(), timepoints_n = n_distinct(canon), exposure = paste(sort(unique(abx)), collapse = ","), .groups = "drop")
if (any(per_subject$samples_n != 16L) || any(per_subject$timepoints_n != 16L)) stop("Incomplete 16-time-point grid")

n_all <- nrow(per_subject); samples_all <- nrow(meta)
n_exp <- sum(per_subject$exposure == "Y"); samples_exp <- sum(meta$abx == "Y")
n_ctrl <- sum(per_subject$exposure == "N"); samples_ctrl <- sum(meta$abx == "N")
stopifnot(n_all == 65L, samples_all == 1040L, n_exp == 60L, samples_exp == 960L, n_ctrl == 5L, samples_ctrl == 80L)

flow <- tibble(
  stage = c("Available longitudinal dataset", "Ciprofloxacin-exposed analysis cohort", "Unexposed descriptive reference", "Original LDR illustration", "Middle participants", "Original HDR illustration"),
  subjects_n = c(65, 60, 5, 15, 30, 15),
  samples_n = c(1040, 960, 80, 240, 480, 240),
  role = c("All supplied participants", "Primary full-cohort analysis", "Descriptive only", "Bottom 15 by C4-C6 Shannon", "Retained in full-cohort analysis", "Top 15 by C4-C6 Shannon"),
  exclusion_or_note = c(
    "Complete 16-sample grid; no local participant/sample exclusion",
    "Selected by exposure status; all 16 samples retained",
    "Not used for confirmatory inference about antibiotic specificity",
    "Secondary extreme-group illustration; not an exclusion from the 60-person cohort",
    "Excluded only from the original 15+15 illustration",
    "Secondary extreme-group illustration; not an exclusion from the 60-person cohort"
  )
)
write_csv(flow, file.path(table_dir, "R1_4_cohort_flow.csv"))
write_csv(per_subject, file.path(table_dir, "R1_4_subject_sample_audit.csv"))

sample_check <- meta |>
  count(abx, canon, name = "samples_n") |>
  arrange(abx, canon)
write_csv(sample_check, file.path(table_dir, "R1_4_timepoint_sample_counts.csv"))

box_plot <- function(title, line1, line2, border) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "white", color = border, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.70, label = title, fontface = "bold", color = border, size = 4.0) +
    annotate("text", x = 0.5, y = 0.43, label = line1, size = 3.5) +
    annotate("text", x = 0.5, y = 0.20, label = line2, size = 3.0, color = "#444444") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + theme_void()
}

all_box <- box_plot("Available longitudinal dataset", "65 participants · 1,040 samples", "16 samples per participant", "#294C60")
exp_box <- box_plot("Ciprofloxacin-exposed cohort", "60 participants · 960 samples", "Primary full-cohort analysis", "#0072B2")
ctrl_box <- box_plot("Unexposed reference", "5 participants · 80 samples", "Descriptive only", "#D55E00")
ldr_box <- box_plot("LDR illustration", "15 participants · 240 samples", "Bottom 15 by endpoint Shannon", "#0072B2")
mid_box <- box_plot("Middle participants", "30 participants · 480 samples", "Retained in full-cohort analysis", "#666666")
hdr_box <- box_plot("HDR illustration", "15 participants · 240 samples", "Top 15 by endpoint Shannon", "#D55E00")

top <- all_box + plot_spacer() + plot_layout(widths = c(2, 0.01))
middle <- exp_box | ctrl_box
bottom <- ldr_box | mid_box | hdr_box
note <- ggplot() +
  annotate("text", x = 0.5, y = 0.66,
           label = "The 15 + 15 groups are subsets of the 60 exposed participants, not cohort exclusions.", size = 3.1) +
  annotate("text", x = 0.5, y = 0.28,
           label = "Both the source report and supplied data identify 5 unexposed participants.",
           size = 2.8, color = "#444444") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

flow_fig <- top / middle / bottom / note +
  plot_layout(heights = c(1.0, 1.0, 1.0, 0.65)) +
  plot_annotation(title = "Human cohort and analysis flow") &
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

ggsave(file.path(figure_dir, "R1_4_cohort_flow.pdf"), flow_fig, width = 190, height = 170, units = "mm", bg = "white", useDingbats = FALSE)
ggsave(file.path(figure_dir, "R1_4_cohort_flow.tif"), flow_fig, width = 190, height = 170, units = "mm", dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(figure_dir, "R1_4_cohort_flow_preview.png"), flow_fig, width = 190, height = 170, units = "mm", dpi = 160, bg = "white")

capture.output(sessionInfo(), file = file.path(log_dir, "R1_4_R_sessionInfo.txt"))
write_lines(c(
  "R1.4 cohort audit passed",
  "65 participants / 1,040 samples / complete 16-sample grid",
  "60 exposed / 960 samples; 5 unexposed / 80 samples",
  "15 LDR + 30 middle + 15 HDR = 60 exposed participants",
  "Author-confirmed: source report and supplied data both identify n=5 unexposed controls"
), file.path(log_dir, "R1_4_run_summary.txt"))

