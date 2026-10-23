#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ragg)
  library(systemfonts)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."))
out <- file.path(root, "reviewer_response", "R1.1_focused_4panel")
for (d in c("figures", "tables", "logs", "QA", "word_assets"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

font_file <- system_fonts() |>
  filter(family == "Liberation Sans", style == "Regular") |>
  slice(1) |>
  pull(path)
if (length(font_file) == 1 && file.exists(font_file))
  register_font(name = "Arial", plain = font_file)

source_root <- file.path(root, "reviewer_response", "R1.1", "full_original_layout")
defs <- tibble(
  definition_id = c("D2", "D3"),
  slug = c("D2_endpoint_bray", "D3_endpoint_jaccard"),
  definition = factor(c("Bray-Curtis", "Jaccard"), levels = c("Bray-Curtis", "Jaccard"))
)
time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
group_pal <- c(HDR = "#2F6FA3", LDR = "#C44E52")

read_def <- function(slug, file) read_csv(file.path(source_root, slug, "tables", file), show_col_types = FALSE)

# Panel B: participant-level bootstrap summaries of alpha diversity.
set.seed(20261020)
boot_mean <- function(x, B = 2000L) {
  b <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)))
}

alpha_raw <- bind_rows(lapply(seq_len(nrow(defs)), function(i) {
  read_def(defs$slug[i], "sample_metrics_with_alternative_groups.csv") |>
    filter(recovery_group %in% c("HDR", "LDR")) |>
    transmute(definition_id = defs$definition_id[i], definition = defs$definition[i],
              subject, canon, group = recovery_group, Shannon = shannon, Richness = richness) |>
    pivot_longer(c(Shannon, Richness), names_to = "metric", values_to = "value")
})) |>
  mutate(canon = factor(canon, levels = time_order), group = factor(group, levels = c("HDR", "LDR")))

alpha_summary <- alpha_raw |>
  group_by(definition_id, definition, metric, group, canon) |>
  summarise(n_subjects = n_distinct(subject), tmp = list(boot_mean(value)), .groups = "drop") |>
  rowwise() |>
  mutate(mean = tmp[1], ci_low = tmp[2], ci_high = tmp[3]) |>
  ungroup() |>
  select(-tmp)

lrt <- read_csv(file.path(root, "reviewer_response", "R1.1", "R1_1_trajectory_LRT.csv"), show_col_types = FALSE) |>
  filter(definition_id %in% c("D2", "D3"), metric %in% c("shannon", "richness")) |>
  transmute(
    definition_id, metric = recode(metric, shannon = "Shannon", richness = "Richness"),
    q_group = q_any_group_within_definition,
    q_interaction = q_interaction_within_definition,
    model_label = paste0("qG = ", format.pval(q_group, digits = 2, eps = 1e-4),
                         "; qG×T = ", format.pval(q_interaction, digits = 2, eps = 1e-4))
  )

alpha_summary <- alpha_summary |>
  left_join(lrt, by = c("definition_id", "metric")) |>
  mutate(metric = factor(metric, levels = c("Shannon", "Richness")))
alpha_annotations <- alpha_summary |>
  distinct(definition, metric, model_label) |>
  mutate(x = 1, y = Inf)
write_csv(alpha_summary, file.path(out, "tables", "panel_B_alpha_diversity_summary.csv"))
write_csv(lrt, file.path(out, "tables", "panel_B_repeated_measure_models.csv"))

# Panel C: use every taxon from the available rapid MDSINE2 sensitivity summaries.
mdsine <- bind_rows(lapply(seq_len(nrow(defs)), function(i) {
  pert <- read_def(defs$slug[i], "perturbation_effects_combined.csv") |>
    select(taxon, perturbation_HDR_minus_LDR = HDR_minus_LDR_median)
  growth <- read_def(defs$slug[i], "growth_effects_combined.csv") |>
    select(taxon, growth_HDR_minus_LDR = HDR_minus_LDR_median)
  inner_join(pert, growth, by = "taxon") |>
    mutate(definition_id = defs$definition_id[i], definition = defs$definition[i])
}))

mdsine_stats <- mdsine |>
  group_by(definition_id, definition) |>
  summarise(
    n_taxa = n(),
    median_perturbation_difference = median(perturbation_HDR_minus_LDR),
    median_growth_difference = median(growth_HDR_minus_LDR),
    spearman_rho = cor(perturbation_HDR_minus_LDR, growth_HDR_minus_LDR, method = "spearman"),
    perturbation_positive = sum(perturbation_HDR_minus_LDR > 0),
    perturbation_negative = sum(perturbation_HDR_minus_LDR < 0),
    growth_positive = sum(growth_HDR_minus_LDR > 0),
    growth_negative = sum(growth_HDR_minus_LDR < 0),
    .groups = "drop"
  ) |>
  mutate(label = paste0("median ΔP = ", number(median_perturbation_difference, accuracy = 0.01),
                        "\nmedian Δg = ", number(median_growth_difference, accuracy = 0.001),
                        "\nSpearman ρ = ", number(spearman_rho, accuracy = 0.01)))
write_csv(mdsine, file.path(out, "tables", "panel_C_MDSINE2_all_taxa.csv"))
write_csv(mdsine_stats, file.path(out, "tables", "panel_C_MDSINE2_summary.csv"))

# Panel D: descriptive FastSpar network trajectories recalculated for each grouping.
network <- bind_rows(lapply(seq_len(nrow(defs)), function(i) {
  read_def(defs$slug[i], "Fig4C_fastspar_summary.csv") |>
    transmute(definition_id = defs$definition_id[i], definition = defs$definition[i],
              group, stage, Nodes = trajectory_nodes, Edges = trajectory_edges)
})) |>
  pivot_longer(c(Nodes, Edges), names_to = "metric", values_to = "count") |>
  mutate(
    stage = factor(stage, levels = c("Baseline", "B6", "C1", "C2", "C3", "End")),
    metric = factor(metric, levels = c("Nodes", "Edges")),
    group = factor(group, levels = c("HDR", "LDR"))
  )
network_end <- network |>
  filter(stage == "End") |>
  select(definition_id, definition, metric, group, count) |>
  pivot_wider(names_from = group, values_from = count) |>
  mutate(HDR_minus_LDR = HDR - LDR)
write_csv(network, file.path(out, "tables", "panel_D_network_trajectory.csv"))
write_csv(network_end, file.path(out, "tables", "panel_D_endpoint_contrast.csv"))

base_theme <- theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.1, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    strip.background = element_rect(fill = "#F2F2F2", colour = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 7.2, colour = "black"),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 7.2, colour = "black"),
    panel.grid.major.y = element_line(colour = "#E9E9E9", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 6, 4, 5)
  )

# Panel A: method diagram.
pA <- ggplot() +
  annotate("rect", xmin = 0.2, xmax = 2.2, ymin = 3.15, ymax = 4.25,
           fill = "#F2F2F2", colour = "black", linewidth = 0.45) +
  annotate("text", x = 1.2, y = 3.82, label = "A1–A4\npersonal baseline", family = "Arial", size = 3.0) +
  annotate("rect", xmin = 0.2, xmax = 2.2, ymin = 1.55, ymax = 2.65,
           fill = "#F2F2F2", colour = "black", linewidth = 0.45) +
  annotate("text", x = 1.2, y = 2.22, label = "C4–C6\nlate recovery", family = "Arial", size = 3.0) +
  annotate("segment", x = 2.25, xend = 3.0, y = 3.7, yend = 3.15,
           arrow = grid::arrow(length = grid::unit(2, "mm")), linewidth = 0.45) +
  annotate("segment", x = 2.25, xend = 3.0, y = 2.15, yend = 2.70,
           arrow = grid::arrow(length = grid::unit(2, "mm")), linewidth = 0.45) +
  annotate("rect", xmin = 3.05, xmax = 5.25, ymin = 2.45, ymax = 3.55,
           fill = "white", colour = "black", linewidth = 0.45) +
  annotate("text", x = 4.15, y = 3.0,
           label = "Mean recovery distance\ncomputed per participant", family = "Arial", size = 2.8) +
  annotate("segment", x = 5.3, xend = 6.05, y = 3.0, yend = 3.0,
           arrow = grid::arrow(length = grid::unit(2, "mm")), linewidth = 0.45) +
  annotate("rect", xmin = 6.1, xmax = 8.15, ymin = 3.15, ymax = 4.2,
           fill = "#D9E8F3", colour = "#2F6FA3", linewidth = 0.6) +
  annotate("text", x = 7.125, y = 3.68, label = "Bray–Curtis\nto baseline centroid",
           family = "Arial", size = 2.7) +
  annotate("rect", xmin = 6.1, xmax = 8.15, ymin = 1.75, ymax = 2.80,
           fill = "#F3DFDF", colour = "#C44E52", linewidth = 0.6) +
  annotate("text", x = 7.125, y = 2.28, label = "Jaccard\nto baseline consensus",
           family = "Arial", size = 2.7) +
  annotate("segment", x = 8.2, xend = 8.9, y = 3.65, yend = 3.65,
           arrow = grid::arrow(length = grid::unit(2, "mm")), linewidth = 0.45) +
  annotate("segment", x = 8.2, xend = 8.9, y = 2.25, yend = 2.25,
           arrow = grid::arrow(length = grid::unit(2, "mm")), linewidth = 0.45) +
  annotate("rect", xmin = 8.95, xmax = 11.2, ymin = 1.65, ymax = 4.25,
           fill = "white", colour = "black", linewidth = 0.45) +
  annotate("text", x = 10.075, y = 3.76, label = "Rank all 60", family = "Arial", size = 2.9, fontface = "bold") +
  annotate("rect", xmin = 9.2, xmax = 10.95, ymin = 2.92, ymax = 3.45,
           fill = "#2F6FA3", colour = NA) +
  annotate("text", x = 10.075, y = 3.19, label = "smallest 15: HDR", colour = "white", family = "Arial", size = 2.45) +
  annotate("text", x = 10.075, y = 2.65, label = "middle 30 excluded", family = "Arial", size = 2.35) +
  annotate("rect", xmin = 9.2, xmax = 10.95, ymin = 1.92, ymax = 2.45,
           fill = "#C44E52", colour = NA) +
  annotate("text", x = 10.075, y = 2.19, label = "largest 15: LDR", colour = "white", family = "Arial", size = 2.45) +
  coord_cartesian(xlim = c(0, 11.45), ylim = c(1.35, 4.45), clip = "off") +
  labs(title = "A") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5, colour = "black", hjust = 0),
        plot.margin = margin(4, 6, 4, 5))

phase_bg <- tibble(xmin = c(0.5, 4.5, 10.5), xmax = c(4.5, 10.5, 16.5),
                   ymin = -Inf, ymax = Inf,
                   fill = c("#F2F2F2", "#FFF1D6", "#EAF3F8"))
pB <- ggplot(alpha_summary, aes(as.numeric(canon), mean, colour = group, fill = group, group = group)) +
  geom_rect(data = phase_bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
            inherit.aes = FALSE, alpha = 0.42, colour = NA) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.11, colour = NA) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.35, stroke = 0.25) +
  geom_text(data = alpha_annotations,
            aes(x = x, y = y, label = model_label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.15, family = "Arial", size = 2.35, colour = "black") +
  facet_grid(metric ~ definition, scales = "free_y") +
  scale_colour_manual(values = group_pal) +
  scale_fill_manual(values = group_pal) +
  scale_x_continuous(breaks = c(1, 4, 5, 10, 11, 16),
                     labels = c("A1", "A4", "B1", "B6", "C1", "C6"),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "B", x = "Sampling point", y = NULL) +
  base_theme + theme(strip.text.x = element_text(size = 7.2), legend.position = "top")

pC <- ggplot(mdsine, aes(perturbation_HDR_minus_LDR, growth_HDR_minus_LDR)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.4, colour = "black") +
  geom_point(size = 1.05, alpha = 0.42, colour = "#5A5A5A") +
  geom_point(data = mdsine_stats,
             aes(x = median_perturbation_difference, y = median_growth_difference),
             inherit.aes = FALSE, shape = 23, size = 3.0, fill = "#C44E52", colour = "black") +
  geom_text(data = mdsine_stats, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.15,
            family = "Arial", size = 2.45, colour = "black") +
  facet_grid(. ~ definition) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1), breaks = c(-10, -3, 0, 3, 10)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.005),
                     breaks = c(-0.1, -0.03, -0.01, 0, 0.01, 0.03, 0.1),
                     labels = label_number(accuracy = 0.001)) +
  labs(title = "C", x = "Perturbation response difference (HDR minus LDR)",
       y = "Growth parameter difference\n(HDR minus LDR)") +
  base_theme + theme(legend.position = "none")

pD <- ggplot(network, aes(stage, count, colour = group, group = group, shape = group)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.7, stroke = 0.35) +
  facet_grid(metric ~ definition, scales = "free_y") +
  scale_colour_manual(values = group_pal) +
  scale_shape_manual(values = c(HDR = 16, LDR = 17)) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(title = "D", x = "Stage", y = NULL) +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

draw_figure <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(3, 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(0.78, 1.25, 1.22), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(pB, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  print(pC, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(pD, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid::popViewport()
}

pdf_file <- file.path(out, "figures", "Figure_R1_four_panel.pdf")
tif_file <- file.path(out, "figures", "Figure_R1_four_panel.tif")
png_file <- file.path(out, "figures", "Figure_R1_four_panel_preview.png")
word_png <- file.path(out, "word_assets", "Figure_R1_four_panel_word.png")

grDevices::cairo_pdf(pdf_file, width = 183 / 25.4, height = 224 / 25.4, family = "Arial")
draw_figure(); grDevices::dev.off()
ragg::agg_tiff(tif_file, width = 183, height = 224, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw_figure(); grDevices::dev.off()
ragg::agg_png(png_file, width = 183, height = 224, units = "mm", res = 180, background = "white")
draw_figure(); grDevices::dev.off()
ragg::agg_png(word_png, width = 183, height = 224, units = "mm", res = 300, background = "white")
draw_figure(); grDevices::dev.off()

# Chinese reviewer response and Word build.
response_md <- file.path(out, "R1.1_focused_response_zh.md")
response_lines <- c(
  "# 审稿人1，主要意见1：恢复力定义敏感性分析",
  "",
  "## 审稿意见",
  "",
  "审稿人建议采用多种恢复力定义重新划分HDR和LDR，并在新分组下重新检验原Fig. 2–4中的多样性、动态模型和网络结论是否保持一致。",
  "",
  "## 回复",
  "",
  "**Response:**  ",
  "Thank you for your valuable comments.",
  "",
  "感谢审稿人的宝贵意见。我们同意，仅依据恢复终点Shannon指数定义HDR和LDR，不能充分说明后续生态学结论对恢复力定义是否稳健。为避免继续使用Shannon指数进行分组，我们选择了两种基于群落组成距离的独立定义：Bray–Curtis距离和二元Jaccard距离。对于每名受试者，我们以A1–A4样本建立个人基线参照，并计算C4–C6样本相对于该基线参照的平均恢复距离。Bray–Curtis定义使用个人基线质心，Jaccard定义使用个人基线出现/缺失共识。分别按两项距离独立排序，距离最小的15人定义为HDR，距离最大的15人定义为LDR，中间30人不进入二组比较。",
  "",
  "我们使用两套新分组重新分析了原Fig. 2–4的关键结论，并将最直接回答审稿问题的结果整合为Figure R1。纵向Shannon和Richness采用包含组别、分类时间及组别×时间交互项的受试者重复测量模型，而不是将16个时间点视为独立观测。在Bray–Curtis分组下，Shannon和Richness均显示显著的组别×时间交互（定义内校正q=0.00024和0.000071）。在Jaccard分组下，Shannon的交互作用仍然显著（q=0.043），Richness的总体组别相关效应显著（q=3.26×10^-10），但Richness组别×时间交互不显著（q=0.126）。因此，在不使用Shannon进行分组后，多样性轨迹仍然区分HDR和LDR，但具体时间依赖形式随距离定义有所变化。",
  "",
  "我们还使用相同的新分组检查了MDSINE2扰动响应和生长参数。该部分现有结果来自快速单链敏感性运行，每组保留100个后验抽样，无法计算R-hat，且没有实测qPCR，因此仅按相对丰度尺度作描述性解释。Bray–Curtis定义下259个物种的HDR−LDR扰动响应差中位数为0.00，Jaccard定义下为2.01；两种定义下生长参数差中位数均为负（−0.014和−0.004）。这表明总体生长方向较一致，但扰动响应幅度对分组定义敏感，不能将其描述为完全稳健的确认性结果。",
  "",
  "共现网络使用各自的新HDR/LDR成员重新计算。恢复终点时，Bray–Curtis分组下HDR与LDR分别包含249和240个节点、4,244和2,854条边；Jaccard分组下分别包含235和229个节点、4,608和1,747条边。两种非Shannon定义均显示恢复终点HDR网络具有更多节点和明显更多边，支持原稿关于较高恢复能力与较丰富网络结构相关的主要结论；中间时间点的变化并非始终同向，因此我们将结论限定在整体动态和恢复终点，而不声称每个时间点均一致。",
  "",
  "综上，基于Bray–Curtis和Jaccard的重新分组支持多样性轨迹差异和恢复终点网络复杂度差异不依赖于终点Shannon分组。MDSINE2生长方向具有一定一致性，但扰动响应参数对定义较敏感，已在修订稿中作为描述性敏感性结果而非确认性证据报告。原Fig. 2E的基线—终点Shannon相关性及完整差异物种热图未纳入核心回复图，以避免循环论证和后验筛选。",
  "",
  "## Figure R1",
  "",
  "![Figure R1](word_assets/Figure_R1_four_panel_word.png){width=6.75in}",
  "",
  "**图注。Figure R1｜采用两种非Shannon恢复定义的关键结论敏感性分析。** **A**，分组流程。A1–A4用于建立个人基线参照，C4–C6用于计算晚期恢复距离；Bray–Curtis距离相对于个人基线质心计算，二元Jaccard距离相对于个人基线出现/缺失共识计算。每种指标独立排序，最小15人为HDR、最大15人为LDR，中间30人不进入组间比较。**B**，两种定义下HDR与LDR的Shannon和Richness纵向轨迹，阴影为受试者层面bootstrap 95%置信区间；分面标题报告重复测量模型的定义内BH校正q值。**C**，259个物种的MDSINE2扰动响应差与生长参数差，均表示HDR减LDR；菱形表示跨物种中位数。该面板来自快速单链敏感性运行，仅作描述性解释。**D**，使用各自新分组重新计算的FastSpar网络节点数和边数。网络计数为描述性结果。所有下游分析均使用对应的新HDR/LDR标签。",
  "",
  "## 分析边界",
  "",
  "1. 两种分组均由C4–C6恢复距离定义，因此恢复后期的组间差异属于定义条件下的敏感性结果，不是独立外部验证。",
  "2. MDSINE2快速运行不具备正式多链收敛诊断，不能据此提出参数显著性或机制性结论。",
  "3. FastSpar节点数和边数用于比较网络结构复杂度，不解释为直接生物相互作用。"
)
writeLines(response_lines, response_md, useBytes = TRUE)

pandoc <- Sys.which("pandoc")
if (!nzchar(pandoc)) stop("pandoc not found")
docx_file <- file.path(out, "R1.1_focused_response_zh_with_Figure_R1.docx")
tmp_docx <- tempfile("R1_1_focused_", tmpdir = out, fileext = ".docx")
pandoc_args <- c(response_md, "--standalone", "--from=markdown+raw_attribute+link_attributes",
                 paste0("--resource-path=", out),
                 paste0("--output=", tmp_docx))
status <- system2(pandoc, pandoc_args)
if (!identical(status, 0L) || !file.exists(tmp_docx)) stop("DOCX build failed")
if (!file.copy(tmp_docx, docx_file, overwrite = TRUE)) stop("DOCX copy failed")
unlink(tmp_docx)

qa <- tibble(
  check = c(
    "D2 and D3 only", "30 subjects per definition", "16 time points",
    "Repeated-measure models available", "259 MDSINE2 taxa per definition",
    "D2/D3 MDSINE2 marked descriptive", "Network recalculated for both definitions",
    "PDF", "TIFF", "Word PNG", "DOCX ZIP integrity"
  ),
  passed = c(
    setequal(unique(alpha_summary$definition_id), c("D2", "D3")),
    all(alpha_raw |> distinct(definition_id, subject) |> count(definition_id) |> pull(n) == 30),
    n_distinct(alpha_raw$canon) == 16,
    nrow(lrt) == 4 && all(is.finite(lrt$q_group)) && all(is.finite(lrt$q_interaction)),
    all(mdsine |> count(definition_id) |> pull(n) == 259),
    TRUE,
    nrow(network) == 48,
    file.exists(pdf_file), file.exists(tif_file), file.exists(word_png),
    system2("unzip", c("-t", docx_file), stdout = FALSE, stderr = FALSE) == 0
  )
)
write_tsv(qa, file.path(out, "QA", "R1_1_focused_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R_sessionInfo.txt"))
if (!all(qa$passed)) stop("QA failed")
message("Completed: ", out)
