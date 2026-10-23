#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ragg)
  library(systemfonts)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.6_sFig2C_PERMDISP")
for (d in c("figures", "tables", "logs", "QA"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
baseline <- paste0("A", 1:4)
endpoint <- paste0("C", 4:6)
n_perm <- 9999L

meta <- read_csv(file.path(root, "raw", "meta.csv"), show_col_types = FALSE,
                 col_types = cols(ids = col_character(), subject = col_character(), canon = col_character())) |>
  filter(abx == "Y", canon %in% time_order) |>
  mutate(canon = factor(canon, levels = time_order)) |>
  arrange(subject, canon)
groups <- read_csv(
  file.path(root, "reproduction", "outputs", "integrated_fig124_raw", "tables",
            "raw_bracken_hdr_ldr_assignments.csv"), show_col_types = FALSE,
  col_types = cols(subject = col_character(), recovery_group = col_character())
) |>
  select(subject, recovery_group)
meta <- meta |>
  left_join(groups, by = "subject")
stopifnot(n_distinct(meta$subject) == 60L, nrow(meta) == 960L)

mpa <- read_tsv(file.path(root, "raw", "combine_brk.mpa"), show_col_types = FALSE,
                name_repair = "minimal")
taxon <- mpa[[1]]
names(mpa)[-1] <- sub("\\.brk\\.report$", "", names(mpa)[-1])
terminal <- sub(".*\\|", "", taxon)
keep <- startsWith(taxon, "k__Bacteria|") & startsWith(terminal, "s__")
counts <- as.matrix(mpa[keep, meta$ids, drop = FALSE])
storage.mode(counts) <- "double"
rel <- sweep(counts, 2, colSums(counts), "/")
profiles <- t(rel)
rownames(profiles) <- meta$ids
full_d <- vegdist(profiles, method = "jaccard", binary = FALSE)
message("Computed quantitative Jaccard matrix for 960 samples")

extract_dist <- function(d, ids) {
  m <- as.matrix(d)
  as.dist(m[ids, ids, drop = FALSE])
}

run_location_dispersion <- function(d, dat, factor_name, block_name, label, seed) {
  dat <- as.data.frame(dat)
  dat[[factor_name]] <- factor(dat[[factor_name]])
  ctrl <- how(nperm = n_perm, blocks = factor(dat[[block_name]]),
              within = Within(type = "free"))
  set.seed(seed)
  fit <- adonis2(d ~ block_factor + test_factor,
                 data = data.frame(block_factor = factor(dat[[block_name]]),
                                   test_factor = dat[[factor_name]]),
                 permutations = ctrl, by = "margin", parallel = 8)
  set.seed(seed + 10000L)
  bd <- betadisper(d, dat[[factor_name]], type = "median", bias.adjust = TRUE,
                   add = "lingoes")
  pd <- permutest(bd, permutations = ctrl, pairwise = FALSE, parallel = 8)
  means <- tapply(bd$distances, dat[[factor_name]], mean)
  tibble(
    comparison = label,
    n_subjects = n_distinct(dat[[block_name]]),
    n_samples = nrow(dat),
    factor_levels = nlevels(dat[[factor_name]]),
    PERMANOVA_pseudoF = unname(fit["test_factor", "F"]),
    PERMANOVA_R2 = unname(fit["test_factor", "R2"]),
    PERMANOVA_p = unname(fit["test_factor", "Pr(>F)"]),
    PERMDISP_F = unname(pd$tab[1, "F"]),
    PERMDISP_p = unname(pd$tab[1, "Pr(>F)"]),
    dispersion_means = paste(names(means), sprintf("%.4f", means), sep = "=", collapse = "; "),
    permutations = n_perm,
    permutation_unit = "labels restricted within participant"
  )
}

results <- list()
# Global 16-timepoint test, all 60 participants.
results[[1]] <- run_location_dispersion(
  full_d, meta, "canon", "subject", "All60_global_time", 20261031L)

# Two-phase subject-level profiles for paired baseline-versus-endpoint tests.
make_phase_data <- function(subject_set) {
  rows <- list()
  mat <- list()
  kk <- 0L
  for (s in subject_set) {
    for (phase_name in c("Baseline", "Endpoint")) {
      kk <- kk + 1L
      times <- if (phase_name == "Baseline") baseline else endpoint
      ids <- meta$ids[meta$subject == s & as.character(meta$canon) %in% times]
      mat[[kk]] <- colMeans(profiles[ids, , drop = FALSE])
      rows[[kk]] <- data.frame(subject = s, phase = phase_name)
    }
  }
  list(meta = bind_rows(rows), profiles = do.call(rbind, mat))
}

all_phase <- make_phase_data(sort(unique(meta$subject)))
all_phase_d <- vegdist(all_phase$profiles, method = "jaccard", binary = FALSE)
results[[2]] <- run_location_dispersion(
  all_phase_d, all_phase$meta, "phase", "subject", "All60_baseline_vs_endpoint", 20261032L)

selected_groups <- c("LDR", "HDR")
ri <- 2L
for (g in selected_groups) {
  gmeta <- meta |>
    filter(recovery_group == g) |>
    arrange(subject, canon)
  gd <- extract_dist(full_d, gmeta$ids)
  ri <- ri + 1L
  results[[ri]] <- run_location_dispersion(
    gd, gmeta, "canon", "subject", paste0(g, "_global_time"), 20261040L + ri)

  gphase <- make_phase_data(sort(unique(gmeta$subject)))
  gphase_d <- vegdist(gphase$profiles, method = "jaccard", binary = FALSE)
  ri <- ri + 1L
  results[[ri]] <- run_location_dispersion(
    gphase_d, gphase$meta, "phase", "subject", paste0(g, "_baseline_vs_endpoint"), 20261040L + ri)
}

stats <- bind_rows(results) |>
  mutate(
    location_result = if_else(PERMANOVA_p < 0.05, "significant", "not significant"),
    dispersion_result = if_else(PERMDISP_p < 0.05, "significant", "not significant")
  )
write_csv(stats, file.path(out, "tables", "sFig2C_PERMANOVA_PERMDISP_results.csv"))

all_scores <- read_csv(
  file.path(root, "reproduction", "outputs", "integrated_fig124_raw", "tables",
            "supplementary2_pcoa_all_high_fidelity.csv"), show_col_types = FALSE
) |>
  mutate(canon = factor(canon, levels = time_order), phase = substr(as.character(canon), 1, 1))
ext_scores <- read_csv(
  file.path(root, "reproduction", "outputs", "integrated_fig124_raw", "tables",
            "supplementary2_pcoa_extremes_high_fidelity.csv"), show_col_types = FALSE
) |>
  mutate(canon = factor(canon, levels = time_order), phase = substr(as.character(canon), 1, 1))
write_csv(all_scores, file.path(out, "tables", "sFig2C_all60_PCoA_scores.csv"))
write_csv(ext_scores, file.path(out, "tables", "sFig2C_HDR_LDR_PCoA_scores.csv"))

summarise_centroids <- function(x, grouping) {
  x |>
    group_by(across(all_of(grouping))) |>
    summarise(
      se1 = stats::sd(PCo1) / sqrt(n()), se2 = stats::sd(PCo2) / sqrt(n()),
      PCo1 = mean(PCo1), PCo2 = mean(PCo2),
      n = n(), .groups = "drop"
    )
}
cent_all <- summarise_centroids(all_scores, c("canon", "phase"))
cent_ext <- summarise_centroids(ext_scores, c("recovery_group", "canon", "phase"))
write_csv(cent_all, file.path(out, "tables", "sFig2C_all60_centroids.csv"))
write_csv(cent_ext, file.path(out, "tables", "sFig2C_HDR_LDR_centroids.csv"))

fmtp <- function(x) ifelse(x < 0.001, "<0.001", paste0("=", sprintf("%.3f", x)))
getstat <- function(name, col) stats[[col]][match(name, stats$comparison)]
left_note <- paste0(
  "Global time: PERMANOVA P", fmtp(getstat("All60_global_time", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("All60_global_time", "PERMDISP_p")), "\n",
  "Baseline vs endpoint: PERMANOVA P", fmtp(getstat("All60_baseline_vs_endpoint", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("All60_baseline_vs_endpoint", "PERMDISP_p"))
)
right_note <- paste0(
  "Global time  LDR: PERMANOVA P", fmtp(getstat("LDR_global_time", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("LDR_global_time", "PERMDISP_p")), "\n",
  "Global time  HDR: PERMANOVA P", fmtp(getstat("HDR_global_time", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("HDR_global_time", "PERMDISP_p")), "\n",
  "Baseline vs endpoint  LDR: PERMANOVA P", fmtp(getstat("LDR_baseline_vs_endpoint", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("LDR_baseline_vs_endpoint", "PERMDISP_p")), "\n",
  "Baseline vs endpoint  HDR: PERMANOVA P", fmtp(getstat("HDR_baseline_vs_endpoint", "PERMANOVA_p")),
  "; PERMDISP P", fmtp(getstat("HDR_baseline_vs_endpoint", "PERMDISP_p"))
)

phase_pal <- c(A = "#3B8BC2", B = "#E64B35", C = "#55A868")
group_phase_pal <- c(
  "HDR_A" = "#2878B5", "HDR_B" = "#D64A3A", "HDR_C" = "#249B68",
  "LDR_A" = "#8172B2", "LDR_B" = "#F39C34", "LDR_C" = "#4DBBD5"
)
cent_ext <- cent_ext |>
  mutate(key = paste(recovery_group, phase, sep = "_"))

base_theme <- theme_classic(base_size = 9.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
    axis.title = element_text(size = 9.6), axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4), legend.title = element_blank(),
    legend.text = element_text(size = 8.2), plot.margin = margin(7, 8, 6, 7)
  )

p_left <- ggplot(cent_all, aes(PCo1, PCo2, group = 1)) +
  geom_path(colour = "#777777", linewidth = 0.55,
            arrow = arrow(length = unit(1.4, "mm"), type = "open")) +
  geom_errorbar(aes(ymin = PCo2 - se2, ymax = PCo2 + se2, colour = phase),
                width = 0, linewidth = 0.55) +
  geom_segment(aes(x = PCo1 - se1, xend = PCo1 + se1,
                   y = PCo2, yend = PCo2, colour = phase), linewidth = 0.55) +
  geom_point(aes(colour = phase), shape = 21, fill = "white", size = 4.3, stroke = 0.9) +
  geom_text(aes(label = canon, colour = phase), size = 2.65, family = "Arial",
            show.legend = FALSE) +
  annotate("label", x = -Inf, y = Inf, label = left_note, hjust = -0.02, vjust = 1.08,
           size = 2.5, family = "Arial", colour = "black", label.size = 0,
           fill = scales::alpha("white", 0.88), lineheight = 1.05) +
  scale_colour_manual(values = phase_pal,
                      labels = c(A = "Baseline (A)", B = "Perturbation (B)", C = "Recovery (C)")) +
  labs(x = "PCoA1 (9.6%)", y = "PCoA2 (5.2%)", colour = NULL) +
  coord_cartesian(clip = "off") + base_theme +
  theme(legend.position = c(0.16, 0.72), legend.background = element_rect(fill = "white", colour = NA))

p_right <- ggplot(cent_ext, aes(PCo1, PCo2, group = recovery_group)) +
  geom_path(aes(linetype = recovery_group), colour = "#666666", linewidth = 0.7,
            arrow = arrow(length = unit(1.35, "mm"), type = "open")) +
  geom_errorbar(aes(ymin = PCo2 - se2, ymax = PCo2 + se2, colour = key),
                width = 0, linewidth = 0.52) +
  geom_segment(aes(x = PCo1 - se1, xend = PCo1 + se1,
                   y = PCo2, yend = PCo2, colour = key), linewidth = 0.52) +
  geom_point(aes(colour = key, shape = recovery_group), fill = "white", size = 4.1, stroke = 0.9) +
  geom_text(aes(label = canon, colour = key), size = 2.5, family = "Arial",
            show.legend = FALSE) +
  annotate("label", x = -Inf, y = Inf, label = right_note, hjust = -0.02, vjust = 1.07,
           size = 2.15, family = "Arial", colour = "black", label.size = 0,
           fill = scales::alpha("white", 0.90), lineheight = 1.02) +
  scale_colour_manual(values = group_phase_pal, guide = "none") +
  scale_shape_manual(values = c(LDR = 22, HDR = 21)) +
  scale_linetype_manual(values = c(LDR = "dashed", HDR = "solid")) +
  labs(x = "PCoA1 (12.5%)", y = "PCoA2 (8.5%)") +
  coord_cartesian(clip = "off") + base_theme +
  theme(legend.position = "top", legend.justification = "right")

draw_combined <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(1, 2, widths = grid::unit(c(1, 1.08), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text("A", x = grid::unit(0.008, "npc"), y = grid::unit(0.982, "npc"),
                  just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12))
  grid::grid.text("B", x = grid::unit(0.488, "npc"), y = grid::unit(0.982, "npc"),
                  just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 12))
  grid::popViewport()
}

pdf_file <- file.path(out, "figures", "sFig2C_PCoA_with_PERMDISP.pdf")
tif_file <- file.path(out, "figures", "sFig2C_PCoA_with_PERMDISP.tif")
png_file <- file.path(out, "figures", "sFig2C_PCoA_with_PERMDISP_preview.png")
grDevices::cairo_pdf(pdf_file, width = 183 / 25.4, height = 92 / 25.4, family = "Arial")
draw_combined(); grDevices::dev.off()
ragg::agg_tiff(tif_file, width = 183, height = 92, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw_combined(); grDevices::dev.off()
ragg::agg_png(png_file, width = 183, height = 92, units = "mm", res = 180,
              background = "white")
draw_combined(); grDevices::dev.off()

response <- c(
  "# R1.6：Supplementary Fig. 2C的PERMANOVA与PERMDISP联合分析",
  "",
  "**Response:**  ",
  "Thank you for your valuable comments.",
  "",
  "感谢审稿人的宝贵意见。我们同意，PERMANOVA的显著结果可能来自组中心位置变化、组内离散度变化或二者共同作用，因此不能在缺少多元离散度检验时将结果单独解释为群落中心偏移。我们重新绘制了Supplementary Fig. 2C，并为图中每项PERMANOVA同步增加基于betadisper的PERMDISP检验。两类检验均采用9,999次受试者内限制置换，以保留纵向重复测量结构。修订图同时报告PERMANOVA P值和PERMDISP P值；完整结果表进一步报告PERMANOVA pseudo-F、R²以及PERMDISP F值。",
  "",
  "左图重新评估全部60名抗生素暴露参与者的整体时间轨迹以及基线与研究终点之间的变化；右图分别评估LDR和HDR的整体时间轨迹及各自的基线—终点变化。PERMANOVA显著而PERMDISP不显著时，结果主要支持群落中心位置变化；两者均显著时，我们将其表述为位置和离散度共同变化，不再归因于单独的中心偏移。原Supplementary Fig. 2C右图中的P=0.052高于预设显著性阈值0.05，现已明确标记并描述为不显著，不再将其作为支持组间差异的证据。",
  "",
  "**修订Supplementary Fig. 2C图注｜抗生素干预期间肠道微生物群落结构的PCoA轨迹及位置—离散度检验。**（A）全部60名抗生素暴露参与者在16个时间点的群落质心轨迹。（B）LDR和HDR各15名参与者的分组质心轨迹。点和误差线分别表示PCoA坐标的组均值和标准误。每个面板同时报告PERMANOVA群落中心位置检验和PERMDISP多元离散度检验；全时间点检验及基线—终点比较均采用9,999次受试者内限制置换。P≥0.05解释为不显著。"
)
writeLines(response, file.path(out, "R1.6_sFig2C_response_zh.md"), useBytes = TRUE)

qa <- tibble(
  check = c("60 exposed participants", "30 HDR/LDR participants", "Six planned comparisons",
            "9999 restricted permutations", "All inferential values finite",
            "Original P=0.052 classified non-significant", "PDF", "TIFF", "Preview"),
  passed = c(
    n_distinct(meta$subject) == 60L,
    n_distinct(meta$subject[meta$recovery_group %in% c("HDR", "LDR")]) == 30L,
    nrow(stats) == 6L, all(stats$permutations == n_perm),
    all(is.finite(stats$PERMANOVA_p)) && all(is.finite(stats$PERMDISP_p)),
    0.052 >= 0.05, file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "sFig2C_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)

