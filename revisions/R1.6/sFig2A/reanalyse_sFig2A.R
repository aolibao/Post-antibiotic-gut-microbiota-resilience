#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ragg)
  library(systemfonts)
  library(vegan)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.6_sFig2A_subject_permutation")
for (d in c("figures", "tables", "logs", "QA"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
threshold_percent <- 0.0001
n_perm <- 9999L
n_boot <- 5000L

meta <- read_csv(file.path(root, "raw", "meta.csv"), show_col_types = FALSE,
                 col_types = cols(ids = col_character(), subject = col_character(), canon = col_character())) |>
  filter(abx == "Y", canon %in% time_order) |>
  select(-group)
groups <- read_csv(
  file.path(root, "reproduction", "outputs", "integrated_fig124_raw", "tables",
            "raw_bracken_hdr_ldr_assignments.csv"),
  show_col_types = FALSE,
  col_types = cols(subject = col_character(), recovery_group = col_character())
) |>
  filter(recovery_group %in% c("HDR", "LDR")) |>
  select(subject, group = recovery_group)
meta <- meta |>
  inner_join(groups, by = "subject") |>
  mutate(canon = factor(canon, levels = time_order)) |>
  arrange(group, subject, canon)
stopifnot(n_distinct(meta$subject[meta$group == "HDR"]) == 15L,
          n_distinct(meta$subject[meta$group == "LDR"]) == 15L,
          nrow(meta) == 480L)

mpa <- read_tsv(file.path(root, "raw", "combine_brk.mpa"), show_col_types = FALSE,
                name_repair = "minimal")
taxon <- mpa[[1]]
names(mpa)[-1] <- sub("\\.brk\\.report$", "", names(mpa)[-1])
terminal <- sub(".*\\|", "", taxon)
keep <- startsWith(taxon, "k__Bacteria|") & startsWith(terminal, "g__")
genus_name <- sub("^g__", "", terminal[keep])
counts <- as.matrix(mpa[keep, meta$ids, drop = FALSE])
storage.mode(counts) <- "double"
# Match the source reconstruction: collapse terminal genus names, normalize
# genus-leaf counts to percent, then apply the displayed >0.0001% rule.
counts <- rowsum(counts, group = genus_name, reorder = FALSE)
percent <- sweep(counts, 2, colSums(counts), "/") * 100
presence <- t(percent > threshold_percent)
rownames(presence) <- meta$ids
message("Loaded and thresholded genus profiles: ", nrow(presence), " samples x ", ncol(presence), " genera")

jaccard_cross <- function(x, y) {
  x <- x * 1
  y <- y * 1
  inter <- tcrossprod(x, y)
  union <- rowSums(x) + rep(rowSums(y), each = nrow(x)) - inter
  out <- 1 - inter / union
  out[union == 0] <- 0
  out
}

subjects_ldr <- sort(unique(groups$subject[groups$group == "LDR"]))
subjects_hdr <- sort(unique(groups$subject[groups$group == "HDR"]))
subjects <- c(subjects_ldr, subjects_hdr)
subject_group <- c(rep("LDR", 15), rep("HDR", 15))

get_ids <- function(time_label) {
  z <- meta |>
    filter(as.character(.data$canon) == .env$time_label) |>
    select(subject, ids)
  setNames(z$ids, z$subject)[subjects]
}

a1_ids <- get_ids("A1")
pair_rows <- list()
stat_rows <- list()

for (ii in seq_along(time_order)) {
  target <- time_order[ii]
  target_ids <- get_ids(target)
  target_dist <- jaccard_cross(
    presence[target_ids[subject_group == "LDR"], , drop = FALSE],
    presence[target_ids[subject_group == "HDR"], , drop = FALSE]
  )
  pair_rows[[ii]] <- expand.grid(
    ldr_subject = subjects_ldr,
    hdr_subject = subjects_hdr,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    mutate(canon = target, jaccard = as.vector(target_dist)) |>
    select(canon, ldr_subject, hdr_subject, jaccard)

  # PERMANOVA uses the complete 30-subject distance matrix and permutes the
  # HDR/LDR labels among subjects. The 225 cross-group distances above are
  # descriptive only and never enter the inferential test as replicates.
  all_dist <- jaccard_cross(presence[target_ids, , drop = FALSE],
                            presence[target_ids, , drop = FALSE])
  test_meta <- data.frame(group = factor(subject_group, levels = c("LDR", "HDR")))
  set.seed(20261020L + ii)
  fit <- adonis2(as.dist(all_dist) ~ group, data = test_meta,
                 permutations = n_perm, by = "terms")
  fit_row <- fit["group", , drop = FALSE]

  # Subject bootstrap CI for the descriptive mean cross-group distance.
  set.seed(20262020L + ii)
  boot_mean <- numeric(n_boot)
  for (bb in seq_len(n_boot)) {
    l <- sample(seq_len(15L), 15L, replace = TRUE)
    h <- sample(16:30, 15L, replace = TRUE)
    boot_mean[bb] <- mean(all_dist[l, h, drop = FALSE])
  }
  ci <- unname(quantile(boot_mean, c(0.025, 0.975), type = 6))
  stat_rows[[ii]] <- tibble(
    canon = target, n_subjects = 30L, n_descriptive_pairs = 225L,
    mean_jaccard = mean(target_dist), mean_ci95_low = ci[1], mean_ci95_high = ci[2],
    pseudo_F = unname(fit_row$F), R2 = unname(fit_row$R2),
    permutation_p = unname(fit_row$`Pr(>F)`),
    permutations = n_perm
  )
  message("Completed ", target, ": p=", format(unname(fit_row$`Pr(>F)`), digits = 4))
}

pair_data <- bind_rows(pair_rows) |>
  mutate(canon = factor(canon, levels = time_order))
stats <- bind_rows(stat_rows) |>
  mutate(
    q_BH = p.adjust(permutation_p, method = "BH"),
    canon = factor(canon, levels = time_order),
    label = case_when(
      q_BH < 0.001 ~ "q<0.001",
      TRUE ~ paste0("q=", sprintf("%.3f", q_BH))
    )
  )

write_csv(pair_data, file.path(out, "tables", "sFig2A_intergroup_Jaccard_values.csv"))
write_csv(stats, file.path(out, "tables", "sFig2A_subject_permutation_results.csv"))

plot_data <- pair_data |>
  mutate(phase = case_when(
    canon %in% paste0("A", 1:4) ~ "Baseline",
    canon %in% paste0("B", 1:5) ~ "Antibiotic exposure",
    TRUE ~ "Post-exposure"
  ))
means <- plot_data |>
  group_by(canon) |>
  summarise(mean = mean(jaccard), .groups = "drop")
ann <- stats |>
  mutate(x = as.numeric(canon), y = 0.91)
phase_pal <- c(Baseline = "#4C8DBB", `Antibiotic exposure` = "#C44E52", `Post-exposure` = "#55A868")

p <- ggplot(plot_data, aes(canon, jaccard, fill = phase)) +
  geom_boxplot(width = 0.70, outlier.shape = NA, linewidth = 0.48,
               colour = "black", alpha = 0.92) +
  geom_line(data = means, aes(canon, mean, group = 1), inherit.aes = FALSE,
            colour = "#B22222", linewidth = 0.65, linetype = "dashed") +
  geom_point(data = means, aes(canon, mean), inherit.aes = FALSE,
             colour = "#B22222", size = 1.2) +
  geom_text(data = ann, aes(x, y, label = label), inherit.aes = FALSE,
            angle = 90, hjust = 1, family = "Arial", size = 2.55, colour = "black") +
  annotate("segment", x = 5, xend = 9, y = -0.052, yend = -0.052,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("segment", x = 5, xend = 5, y = -0.052, yend = -0.032,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("segment", x = 9, xend = 9, y = -0.052, yend = -0.032,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("text", x = 7, y = -0.082, label = "Antibiotic exposure",
           colour = "#C44E52", family = "Arial", size = 3.0) +
  scale_fill_manual(values = phase_pal) +
  scale_y_continuous(limits = c(-0.10, 0.94), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Pairwise Jaccard distance (LDR vs HDR)") +
  theme_classic(base_size = 9.5, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8.6, colour = "black"),
    axis.line = element_line(linewidth = 0.48, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 20, 8)
  )

pdf_file <- file.path(out, "figures", "sFig2A_Jaccard_subject_permutation.pdf")
tif_file <- file.path(out, "figures", "sFig2A_Jaccard_subject_permutation.tif")
png_file <- file.path(out, "figures", "sFig2A_Jaccard_subject_permutation_preview.png")
ggsave(pdf_file, p, width = 183, height = 105, units = "mm", device = cairo_pdf,
       bg = "white", family = "Arial")
ggsave(tif_file, p, width = 183, height = 105, units = "mm", dpi = 600,
       device = agg_tiff, compression = "lzw", background = "white")
ggsave(png_file, p, width = 183, height = 105, units = "mm", dpi = 180,
       device = agg_png, background = "white")

response <- c(
  "# R1.6第二点补充：Supplementary Fig. 2A距离检验修订",
  "",
  "**Response:**  ",
  "Thank you for your valuable comments.",
  "",
  "感谢审稿人的宝贵意见。我们同意，Supplementary Fig. 2A中的225个LDR–HDR组间Jaccard距离由同一组受试者重复形成，不能视为相互独立的观测；原先对固定距离对采用Wilcoxon符号秩检验也与实际研究设计不匹配。因此，我们保留箱线图作为组间距离分布的描述性展示，但完全替换了显著性检验。",
  "",
  "修订分析以受试者而非距离对作为置换单位。对于A1至C6的每个时间点，我们基于30名受试者的完整Jaccard距离矩阵进行PERMANOVA，以HDR/LDR为分组变量，并在受试者层面置换组标签9,999次。16个时间点的置换P值进一步采用Benjamini–Hochberg方法校正。图中报告校正后的q值，结果表同时报告PERMANOVA的pseudo-F、R²和置换P值；箱线图及按受试者重抽样获得的95% bootstrap置信区间仅用于描述组间距离。该分析不再把共享同一受试者端点的225个距离当作独立生物学重复。",
  "",
  "**修订Supplementary Fig. 2A图注｜HDR与LDR之间的纵向组间Jaccard距离。** 箱线图描述每个时间点15名LDR与15名HDR受试者形成的225个组间Jaccard距离，中线和箱体分别表示中位数和四分位距，红色虚线连接各时间点均值；这些距离仅用于描述分布。每个时间点的组间差异采用完整30受试者Jaccard距离矩阵的PERMANOVA检验，HDR/LDR组标签在受试者层面置换9,999次，并对16项检验进行BH校正。图中报告校正后的q值。"
)
writeLines(response, file.path(out, "R1.6_sFig2A_response_zh.md"), useBytes = TRUE)

qa <- tibble(
  check = c("15 HDR participants", "15 LDR participants", "16 time points",
            "225 descriptive distances per time", "16 subject-level PERMANOVA tests",
            "Time-specific distance distributions vary", "9999 permutations per test",
            "Finite inferential results", "PDF", "TIFF", "Preview"),
  passed = c(
    length(subjects_hdr) == 15L, length(subjects_ldr) == 15L,
    n_distinct(pair_data$canon) == 16L,
    all(table(pair_data$canon) == 225L), nrow(stats) == 16L,
    sd(stats$mean_jaccard) > 1e-6,
    all(stats$permutations == n_perm),
    all(is.finite(stats$permutation_p)) && all(is.finite(stats$q_BH)),
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "sFig2A_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)

