#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
audit_file <- file.path(root, "reviewer_response/R1.9/tables/R1_9_data_availability_audit.csv")
out_pdf <- file.path(root, "reviewer_response/R1.9/figures/R1_9_evidence_availability_audit.pdf")
out_tif <- file.path(root, "reviewer_response/R1.9/figures/R1_9_evidence_availability_audit.tif")

d <- read.csv(audit_file, check.names = FALSE)
d$item <- factor(d$requested_item, levels = rev(d$requested_item))
d$status <- factor(d$status, levels = c("missing", "reported_only"))

p <- ggplot(d, aes(x = status, y = item, colour = status)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c(missing = "#B44A3A", reported_only = "#C28A2C"), guide = "none") +
  scale_x_discrete(labels = c(missing = "Required data absent", reported_only = "Reported in manuscript only")) +
  labs(
    title = "R1.9 evidence availability audit",
    subtitle = "Internal audit only; animal-level reanalysis still requires source data",
    x = NULL,
    y = NULL,
    caption = "Absence was assessed within the supplied local workspace."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_line(colour = "#E6E8EB"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "#2F343B"),
    axis.text.x = element_text(colour = "#2F343B")
  )

ggsave(out_pdf, p, width = 8.0, height = 5.2, units = "in", device = "pdf")
ggsave(out_tif, p, width = 8.0, height = 5.2, units = "in", dpi = 600, compression = "lzw")

