#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
src_tables <- file.path(root, "reviewer_response", "03_tables")
audit_tables <- file.path(root, "reproduction", "outputs", "integrated_fig124_raw", "tables", "sfig3_variant_audit")
out <- file.path(root, "reviewer_response", "R1.7")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)

cor_file <- file.path(src_tables, "core_compositional_CLR_focal_excluded_correlations.csv")
source_file <- file.path(audit_tables, "sfig3_source_upset_counts.csv")
score_file <- file.path(audit_tables, "sfig3_variant_global_scores.csv")
inventory_file <- file.path(audit_tables, "sfig3_variant_feature_inventory.csv")
stopifnot(all(file.exists(c(cor_file, source_file, score_file, inventory_file))))

cors <- read_csv(cor_file, show_col_types = FALSE)
source_upset <- read_csv(source_file, show_col_types = FALSE)
scores <- read_csv(score_file, show_col_types = FALSE)
inventory <- read_csv(inventory_file, show_col_types = FALSE)

write_csv(cors, file.path(out, "tables", "R1_7_CLR_focal_excluded_correlations.csv"))
write_csv(source_upset, file.path(out, "tables", "R1_7_original_UpSet_source_counts_audit.csv"))
write_csv(scores, file.path(out, "tables", "R1_7_UpSet_reconstruction_scores.csv"))
write_csv(inventory, file.path(out, "tables", "R1_7_feature_inventory.csv"))

primary <- cors |> filter(dataset_branch == "rarefied_primary") |>
  mutate(phase = substr(canon, 1, 1),
         corrected_positive_hit = rho_CLR_vs_focal_genus_excluded_shannon > 0 & q_BH_global_16_times < 0.05)

patterns <- primary |>
  group_by(genus, phase) |>
  summarise(present = any(corrected_positive_hit), .groups = "drop") |>
  complete(genus, phase = c("A", "B", "C"), fill = list(present = FALSE)) |>
  pivot_wider(names_from = phase, values_from = present) |>
  mutate(pattern = paste0(if_else(A, "A", ""), if_else(B, "B", ""), if_else(C, "C", "")),
         pattern = if_else(pattern == "", "None", pattern))
write_csv(patterns, file.path(out, "tables", "R1_7_corrected_phase_membership_by_genus.csv"))

pattern_levels <- c("A", "B", "C", "AB", "AC", "BC", "ABC", "None")
pattern_counts <- patterns |>
  count(pattern, name = "n_genera") |>
  complete(pattern = pattern_levels, fill = list(n_genera = 0)) |>
  mutate(pattern = factor(pattern, levels = pattern_levels),
         definition = "At least one positive CLR-vs-focal-excluded-Shannon association with global BH q < 0.05 in the phase")
write_csv(pattern_counts, file.path(out, "tables", "R1_7_corrected_phase_pattern_counts.csv"))

selected_taxa <- c("Bacteroides", "Phocaeicola", "Faecalibacterium", "Blautia", "Alistipes", "Akkermansia")
c6_selected <- primary |> filter(canon == "C6", genus %in% selected_taxa) |>
  select(genus, n_subjects, rho_relative_abundance_vs_unmodified_shannon,
         rho_CLR_vs_focal_genus_excluded_shannon, p_subject_label_permutation,
         q_BH_within_time, q_BH_global_16_times)
write_csv(c6_selected, file.path(out, "tables", "R1_7_C6_selected_taxa.csv"))

summary_tbl <- tibble(
  item = c("Rarefied genera tested at 16 times", "Corrected genus-time associations with global q < 0.05",
           "Original genus-rank database features", "Best reconstruction relative L1 error",
           "C6 Faecalibacterium corrected rho", "C6 Faecalibacterium global q",
           "C6 Bacteroides corrected rho", "C6 Bacteroides global q"),
  value = c(
    paste0(n_distinct(primary$genus), " x ", n_distinct(primary$canon)),
    as.character(sum(primary$q_BH_global_16_times < 0.05)),
    as.character(inventory |> filter(level == "genus", taxon_mode == "direct_path", abundance_scale == "relative_abundance") |> pull(n_features)),
    format(min(scores$upset_pattern_relative_L1), digits = 4),
    format(c6_selected |> filter(genus == "Faecalibacterium") |> pull(rho_CLR_vs_focal_genus_excluded_shannon), digits = 4),
    format(c6_selected |> filter(genus == "Faecalibacterium") |> pull(q_BH_global_16_times), digits = 4),
    format(c6_selected |> filter(genus == "Bacteroides") |> pull(rho_CLR_vs_focal_genus_excluded_shannon), digits = 4),
    format(c6_selected |> filter(genus == "Bacteroides") |> pull(q_BH_global_16_times), digits = 4)
  )
)
write_csv(summary_tbl, file.path(out, "tables", "R1_7_key_results.csv"))

base_theme <- theme_classic(base_size = 9, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        axis.text = element_text(colour = "black"),
        legend.position = "top", legend.title = element_blank(),
        plot.margin = margin(5, 7, 5, 5))

dumb <- c6_selected |>
  select(genus, conventional = rho_relative_abundance_vs_unmodified_shannon,
         corrected = rho_CLR_vs_focal_genus_excluded_shannon) |>
  pivot_longer(c(conventional, corrected), names_to = "analysis", values_to = "rho") |>
  mutate(analysis = recode(analysis,
    conventional = "Conventional same-vector",
    corrected = "CLR + focal exclusion"))
pA <- ggplot(dumb, aes(rho, reorder(genus, rho))) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_line(aes(group = genus), colour = "grey70", linewidth = 0.8) +
  geom_point(aes(colour = analysis, shape = analysis), size = 2.6, stroke = 0.5) +
  scale_colour_manual(values = c("Conventional same-vector" = "#999999",
                                 "CLR + focal exclusion" = "#0072B2")) +
  scale_shape_manual(values = c("Conventional same-vector" = 1,
                                "CLR + focal exclusion" = 16)) +
  scale_x_continuous(limits = c(-0.8, 0.65), breaks = seq(-0.8, 0.6, 0.2)) +
  labs(title = "A  Coupling sensitivity at C6",
       subtitle = "C6, n = 60 participants", x = expression(Spearman~rho), y = NULL) + base_theme

c6_all <- primary |> filter(canon == "C6") |>
  mutate(global_result = if_else(q_BH_global_16_times < 0.05, "Global q < 0.05", "Other"))
pB <- ggplot(c6_all, aes(rho_relative_abundance_vs_unmodified_shannon,
                         rho_CLR_vs_focal_genus_excluded_shannon,
                         colour = global_result, shape = global_result)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey55", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.35) +
  geom_point(alpha = 0.82, size = 1.8, stroke = 0.35) +
  scale_colour_manual(values = c("Global q < 0.05" = "#D55E00", "Other" = "#777777")) +
  scale_shape_manual(values = c("Global q < 0.05" = 17, "Other" = 16)) +
  coord_equal(xlim = c(-0.8, 0.9), ylim = c(-0.8, 0.9)) +
  labs(title = "B  All 124 genera at C6",
       subtitle = "Global BH across 1,984 tests",
       x = "Conventional correlation", y = "Corrected correlation") + base_theme

source_genus <- source_upset |> filter(level == "genus") |>
  mutate(pattern = factor(pattern, levels = pattern_levels[1:7]))
pC <- ggplot(source_genus, aes(pattern, source_n, fill = group)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = comma(source_n)), position = position_dodge(width = 0.76), vjust = -0.25, size = 2.4) +
  scale_fill_manual(values = c(HDR = "#0072B2", LDR = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "C  Original UpSet counts: clarified unit",
       subtitle = "Source plotting-table values; genus-rank Bracken features\nIncludes low-abundance and provisional labels",
       x = "Mutually exclusive phase-presence pattern", y = "Source figure count") + base_theme

pD <- ggplot(pattern_counts, aes(pattern, n_genera)) +
  geom_col(width = 0.68, fill = "#009E73", colour = "white", linewidth = 0.2) +
  geom_text(aes(label = n_genera), vjust = -0.3, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "D  Prevalence-filtered sensitivity",
       subtitle = "Full cohort; positive corrected associations\nGlobal BH q < 0.05",
       x = "Mutually exclusive phase-presence pattern", y = "Number of genera") + base_theme

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = "R1.7 | Compositionality-aware sensitivity and UpSet clarification",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, margin = margin(b = 5))))

pdf_file <- file.path(out, "figures", "R1_7_response_figure.pdf")
tif_file <- file.path(out, "figures", "R1_7_response_figure.tif")
png_file <- file.path(out, "figures", "R1_7_response_figure_preview.png")
ggsave(pdf_file, fig, width = 210, height = 175, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(tif_file, fig, width = 210, height = 175, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", background = "white")
ggsave(png_file, fig, width = 210, height = 175, units = "mm", dpi = 180,
       device = ragg::agg_png, background = "white")

qa <- tibble(
  check = c("Primary rarefied branch has 124 genera", "Primary branch has 16 time points",
            "Every corrected test uses 60 participants", "Permutation count is 999",
            "Global BH q is populated", "Source genus ABC counts are 862 and 1294",
            "Best reconstruction does not exactly match source", "Corrected pattern counts sum to 124",
            "PDF exists", "TIFF exists", "Preview exists"),
  passed = c(
    n_distinct(primary$genus) == 124, n_distinct(primary$canon) == 16,
    all(primary$n_subjects == 60), all(primary$permutations == 999),
    all(!is.na(primary$q_BH_global_16_times)),
    setequal(source_genus$source_n[source_genus$pattern == "ABC"], c(862, 1294)),
    min(scores$upset_pattern_relative_L1) > 0,
    sum(pattern_counts$n_genera) == 124,
    file.exists(pdf_file), file.exists(tif_file), file.exists(png_file)
  )
)
write_tsv(qa, file.path(out, "QA", "R1_7_QA.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "R1_7_R_sessionInfo.txt"))
if (!all(qa$passed)) stop("R1.7 QA failed")
cat(sprintf("R1.7 complete: %d/%d checks passed\n", sum(qa$passed), nrow(qa)))
