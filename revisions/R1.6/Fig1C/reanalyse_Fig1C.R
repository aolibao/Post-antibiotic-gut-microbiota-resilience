#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(vegan)
  library(permute)
  library(ragg)
  library(systemfonts)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.6_Fig1C_permutation")
for (d in c("figures", "tables", "logs", "QA"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
meta <- read_csv(file.path(root, "raw", "meta.csv"), show_col_types = FALSE,
                 col_types = cols(ids = col_character(), subject = col_character(), canon = col_character())) |>
  filter(abx == "Y", canon %in% time_order) |>
  mutate(canon = factor(canon, levels = time_order)) |>
  arrange(subject, canon)
stopifnot(n_distinct(meta$subject) == 60L, nrow(meta) == 960L)

mpa <- read_tsv(file.path(root, "raw", "combine_brk.mpa"), show_col_types = FALSE,
                name_repair = "minimal")
message("Loaded Bracken MPA table")
taxon <- mpa[[1]]
sample_names <- sub("\\.brk\\.report$", "", names(mpa)[-1])
names(mpa)[-1] <- sample_names
terminal <- sub(".*\\|", "", taxon)
keep <- startsWith(taxon, "k__Bacteria|") & startsWith(terminal, "s__")
counts <- as.matrix(mpa[keep, meta$ids, drop = FALSE])
storage.mode(counts) <- "double"
presence <- t(counts > 0)
rownames(presence) <- meta$ids

# Full Jaccard matrix is computed once; each time-specific paired analysis
# extracts the corresponding A1/time-point block.
full_jaccard <- as.matrix(vegdist(presence, method = "jaccard", binary = TRUE))
message("Computed full binary Jaccard matrix")

a1_id <- meta |>
  filter(canon == "A1") |>
  select(subject, a1_id = ids)
plot_data <- meta |>
  left_join(a1_id, by = "subject") |>
  rowwise() |>
  mutate(jaccard_to_A1 = full_jaccard[ids, a1_id]) |>
  ungroup() |>
  select(ids, subject, canon, jaccard_to_A1)

run_paired_permanova <- function(target, seed) {
  pair_meta <- meta |>
    filter(canon %in% c("A1", target)) |>
    mutate(subject = factor(subject), timepoint = factor(as.character(canon), levels = c("A1", target))) |>
    arrange(subject, timepoint)
  stopifnot(nrow(pair_meta) == 120L, all(table(pair_meta$subject) == 2L))
  d <- as.dist(full_jaccard[pair_meta$ids, pair_meta$ids])
  ctrl <- how(nperm = 9999, blocks = pair_meta$subject,
              within = Within(type = "free"))
  set.seed(seed)
  # A paired/restricted PERMANOVA: participant is the exchangeability block,
  # so A1/target labels can only be swapped within the same participant.
  fit <- adonis2(d ~ timepoint, data = pair_meta, permutations = ctrl,
                 by = "terms")
  row <- fit["timepoint", , drop = FALSE]
  tibble(
    canon = target,
    n_subjects = n_distinct(pair_meta$subject),
    pseudo_F = unname(row$F),
    R2 = unname(row$R2),
    permutation_p = unname(row$`Pr(>F)`),
    permutations = 9999L,
    permutation_scheme = "A1/target labels exchanged only within participant",
    model = "paired PERMANOVA: Jaccard ~ timepoint; participant as permutation block"
  )
}

targets <- time_order[-1]
n_workers <- min(8L, length(targets), parallel::detectCores(logical = FALSE))
message("Running 15 paired PERMANOVA tests with ", n_workers, " forked workers")
tests <- bind_rows(parallel::mclapply(seq_along(targets), function(i)
  run_paired_permanova(targets[i], 20261020L + i), mc.cores = n_workers)) |>
  mutate(
    q_BH = p.adjust(permutation_p, method = "BH"),
    canon = factor(canon, levels = time_order),
    label = if_else(q_BH < 0.001, "q<0.001", paste0("q=", sprintf("%.3f", q_BH)))
  )

plot_data <- plot_data |>
  mutate(
    canon = factor(canon, levels = time_order),
    phase = case_when(
      canon %in% paste0("A", 1:4) ~ "Baseline",
      canon %in% paste0("B", 1:5) ~ "Antibiotic exposure",
      TRUE ~ "Post-exposure"
    )
  )
write_csv(plot_data, file.path(out, "tables", "Fig1C_Jaccard_to_A1_values.csv"))
write_csv(tests, file.path(out, "tables", "Fig1C_paired_PERMANOVA_results.csv"))

ann <- tests |>
  mutate(x = as.numeric(canon), y = 0.88)
ref_ann <- tibble(x = 1, y = 0.88, label = "reference")
means <- plot_data |>
  group_by(canon) |>
  summarise(mean = mean(jaccard_to_A1), .groups = "drop")

phase_pal <- c(Baseline = "#4C8DBB", `Antibiotic exposure` = "#C44E52", `Post-exposure` = "#55A868")
p <- ggplot(plot_data, aes(canon, jaccard_to_A1, fill = phase)) +
  geom_boxplot(width = 0.72, outlier.shape = NA, linewidth = 0.48,
               colour = "black", alpha = 0.92) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20261020),
             shape = 16, size = 0.65, alpha = 0.24, colour = "black") +
  geom_line(data = means, aes(x = canon, y = mean, group = 1),
            inherit.aes = FALSE, colour = "#B22222", linewidth = 0.65,
            linetype = "dashed") +
  geom_point(data = means, aes(x = canon, y = mean), inherit.aes = FALSE,
             colour = "#B22222", size = 1.15) +
  geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            angle = 90, hjust = 1, vjust = 0.5, family = "Arial", size = 2.45, colour = "black") +
  geom_text(data = ref_ann, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            angle = 90, hjust = 1, vjust = 0.5, family = "Arial", size = 2.45, colour = "black") +
  annotate("segment", x = 5, xend = 9, y = -0.055, yend = -0.055,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("segment", x = 5, xend = 5, y = -0.055, yend = -0.035,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("segment", x = 9, xend = 9, y = -0.055, yend = -0.035,
           linewidth = 0.55, colour = "#C44E52") +
  annotate("text", x = 7, y = -0.085, label = "Antibiotic exposure",
           colour = "#C44E52", family = "Arial", size = 3.0) +
  scale_fill_manual(values = phase_pal) +
  scale_y_continuous(limits = c(-0.105, 0.91), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Jaccard distance to the participant's A1 sample") +
  theme_classic(base_size = 9.5, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8.6, colour = "black"),
    axis.text.y = element_text(size = 8.6, colour = "black"),
    axis.line = element_line(linewidth = 0.48, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 20, 8)
  )

pdf_file <- file.path(out, "figures", "Fig1C_Jaccard_paired_PERMANOVA.pdf")
tif_file <- file.path(out, "figures", "Fig1C_Jaccard_paired_PERMANOVA.tif")
png_file <- file.path(out, "figures", "Fig1C_Jaccard_paired_PERMANOVA_preview.png")
ggsave(pdf_file, p, width = 183, height = 105, units = "mm", device = cairo_pdf,
       bg = "white", family = "Arial")
ggsave(tif_file, p, width = 183, height = 105, units = "mm", dpi = 600,
       device = agg_tiff, compression = "lzw", background = "white")
ggsave(png_file, p, width = 183, height = 105, units = "mm", dpi = 180,
       device = agg_png, background = "white")

response <- c(
  "# R1.6第二点：Fig. 1C距离检验修订",
  "",
  "**Response:**  ",
  "Thank you for your valuable comments.",
  "",
  "感谢审稿人的宝贵意见。我们同意，原Fig. 1C中各时间点相对于A1的Jaccard距离共享同一受试者和同一A1参照，将这些距离直接用于普通t检验不能正确反映数据的配对结构。修订后的Fig. 1C仍以箱线图展示全部60名抗生素暴露参与者在各时间点相对于本人A1样本的Jaccard距离，但箱线图中的距离值仅用于描述效应分布，不作为独立观测进行参数检验。",
  "",
  "显著性检验改为受试者内限制置换PERMANOVA。对于A2至C6的每个时间点，我们使用该时间点与A1的二元物种组成矩阵，以受试者作为置换区组，并仅允许A1/目标时间点标签在同一受试者内交换，共进行9,999次置换。15个时间点的置换P值进一步采用Benjamini–Hochberg方法校正。图中标注的是校正后的q值，A1标记为reference，不进行A1与自身的统计检验。该方法以受试者为配对和置换限制单位，避免将共享A1样本的两两Jaccard距离错误地视为相互独立。",
  "",
  "**图注。修订Fig. 1C｜各时间点相对于个人A1样本的Jaccard距离。** 箱线图显示60名抗生素暴露参与者的中位数、四分位距和1.5×四分位距范围；散点为受试者，红色虚线连接各时间点均值。显著性来自A1与各目标时间点二元群落组成的配对PERMANOVA，A1/目标时间标签仅在同一受试者内交换（9,999次置换），并对15项检验进行BH校正。图中报告校正后的q值。"
)
writeLines(response, file.path(out, "R1.6_second_point_response_zh.md"), useBytes = TRUE)

qa <- tibble(
  check = c("60 exposed participants", "16 time points", "960 plotted observations",
            "15 paired PERMANOVA tests", "9999 restricted permutations",
            "All P and q values finite", "A1 distances are zero", "PDF", "TIFF", "Preview"),
  passed = c(
    n_distinct(plot_data$subject) == 60L, n_distinct(plot_data$canon) == 16L,
    nrow(plot_data) == 960L, nrow(tests) == 15L,
    all(tests$permutations == 9999L),
    all(is.finite(tests$permutation_p)) && all(is.finite(tests$q_BH)),
    all(abs(plot_data$jaccard_to_A1[plot_data$canon == "A1"]) < 1e-12),
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "Fig1C_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)
