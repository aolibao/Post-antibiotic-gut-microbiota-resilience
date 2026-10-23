#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(nlme)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_path <- normalizePath(sub("^--file=", "", script_arg))
root <- dirname(dirname(script_path))

out6 <- file.path(root, "reviewer_response", "R1.6_revised")
out7 <- file.path(root, "reviewer_response", "R1.7_revised")
for (d in c(out6, out7)) {
  dir.create(file.path(d, "scripts"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "QA"), recursive = TRUE, showWarnings = FALSE)
}

set.seed(20261013)
time_order <- c(paste0("A", 1:4), paste0("B", 1:6), paste0("C", 1:6))
phase_map <- c(
  A1 = "Baseline", A2 = "Baseline", A3 = "Baseline", A4 = "Baseline",
  B1 = "Antibiotic", B2 = "Antibiotic", B3 = "Antibiotic", B4 = "Antibiotic", B5 = "Antibiotic",
  C4 = "Late recovery", C5 = "Late recovery", C6 = "Late recovery"
)
phase_order <- c("Baseline", "Antibiotic", "Late recovery")

message("Reading source count table...")
meta <- fread(file.path(root, "raw", "meta.csv"), colClasses = list(character = c("ids", "subject", "canon")))
mpa <- fread(file.path(root, "raw", "combine_brk.mpa"), check.names = FALSE)
tax_path <- mpa[[1]]
sample_columns <- names(mpa)[-1]
sample_ids <- sub("\\.brk\\.report$", "", sample_columns)
setnames(mpa, sample_columns, sample_ids)

terminal <- sub("^.*\\|", "", tax_path)
species_keep <- startsWith(tax_path, "k__Bacteria|") & startsWith(terminal, "s__")
species_path <- tax_path[species_keep]
species_label <- gsub("_", " ", sub("^s__", "", terminal[species_keep]))
genus_label <- vapply(strsplit(species_path, "\\|", fixed = FALSE), function(parts) {
  hit <- parts[startsWith(parts, "g__")]
  if (!length(hit)) return("Unclassified")
  ans <- gsub("_", " ", sub("^g__", "", hit[length(hit)]))
  ifelse(nchar(ans) == 0, "Unclassified", ans)
}, character(1))

exposed <- meta[abx == "Y" & canon %in% time_order]
stopifnot(uniqueN(exposed$subject) == 60L, nrow(exposed) == 960L)
missing_ids <- setdiff(exposed$ids, names(mpa))
if (length(missing_ids)) stop("Missing source count columns: ", paste(missing_ids, collapse = ", "))
exposed <- exposed[match(intersect(names(mpa), exposed$ids), ids)]
exposed_ids <- exposed$ids
count_mat <- as.matrix(mpa[species_keep, ..exposed_ids])
storage.mode(count_mat) <- "double"
rownames(count_mat) <- species_path
sample_totals <- colSums(count_mat)

# Collapse terminal duplicate species and genus labels reproducibly.
species_mat <- rowsum(count_mat, group = species_label, reorder = FALSE)
genus_mat <- rowsum(count_mat, group = genus_label, reorder = FALSE)
species_clogc <- count_mat
positive <- species_clogc > 0
species_clogc[positive] <- species_clogc[positive] * log(species_clogc[positive])
species_clogc[!positive] <- 0
total_clogc <- colSums(species_clogc)
genus_clogc <- rowsum(species_clogc, group = genus_label, reorder = FALSE)

meta_dt <- as.data.table(exposed)
meta_dt[, canon := factor(canon, levels = time_order)]
meta_dt[, time_index := match(as.character(canon), time_order)]
meta_dt[, log_depth := log10(sample_totals)]

message("R1.6: participant-level trajectory and gated species-profile analysis...")
assign_file <- file.path(root, "reviewer_response", "R1.1", "full_original_layout",
                         "D1_endpoint_abs_shannon", "tables", "subject_group_assignments_and_panelB_values.csv")
candidate_file <- file.path(root, "reviewer_response", "R1.1", "full_original_layout",
                            "D1_endpoint_abs_shannon", "tables", "Fig2F_LEfSe_candidate_universe.csv")
sample_metric_file <- file.path(root, "reviewer_response", "R1.1", "full_original_layout",
                               "D1_endpoint_abs_shannon", "tables", "sample_metrics_with_alternative_groups.csv")
assign <- fread(assign_file)[recovery_group %in% c("HDR", "LDR"), .(subject, recovery_group)]
stopifnot(nrow(assign) == 30L, all(table(assign$recovery_group) == 15L))
candidates <- fread(candidate_file)$species
stopifnot(length(candidates) == 268L)

subject_phase_meta <- merge(meta_dt, assign, by = "subject")
subject_phase_meta[, phase := unname(phase_map[as.character(canon)])]
subject_phase_meta <- subject_phase_meta[!is.na(phase)]
subject_phase_meta[, phase := factor(phase, levels = phase_order)]

available_candidates <- intersect(candidates, rownames(species_mat))
if (length(available_candidates) != length(candidates)) {
  warning(length(candidates) - length(available_candidates), " candidate species absent from collapsed source table")
}
cpm <- sweep(species_mat[available_candidates, , drop = FALSE], 2, sample_totals, "/") * 1e6
log_cpm <- log10(cpm + 1)

# Phase means are calculated within participant, so the independent unit remains the participant.
phase_rows <- lapply(seq_len(nrow(subject_phase_meta)), function(i) {
  data.frame(
    subject = subject_phase_meta$subject[i],
    group = subject_phase_meta$recovery_group[i],
    phase = as.character(subject_phase_meta$phase[i]),
    sample_id = subject_phase_meta$ids[i],
    stringsAsFactors = FALSE
  )
}) |> bind_rows()
phase_long <- as.data.table(phase_rows)

phase_array <- array(NA_real_, dim = c(length(available_candidates), 30L, 3L),
                     dimnames = list(available_candidates, sort(assign$subject), phase_order))
for (s in sort(assign$subject)) {
  for (ph in phase_order) {
    ids <- phase_long[subject == s & phase == ph, sample_id]
    if (length(ids)) phase_array[, s, ph] <- rowMeans(log_cpm[, ids, drop = FALSE])
  }
}

group_by_subject <- assign$recovery_group[match(dimnames(phase_array)[[2]], assign$subject)]
omnibus_stat <- function(labels) {
  h <- labels == "HDR"
  l <- labels == "LDR"
  out <- numeric(dim(phase_array)[1])
  for (j in seq_len(dim(phase_array)[3])) {
    xh <- phase_array[, h, j, drop = FALSE][, , 1]
    xl <- phase_array[, l, j, drop = FALSE][, , 1]
    mh <- rowMeans(xh)
    ml <- rowMeans(xl)
    vh <- apply(xh, 1, var)
    vl <- apply(xl, 1, var)
    denom <- vh / sum(h) + vl / sum(l)
    out <- out + ifelse(denom > 0, (mh - ml)^2 / denom, 0)
  }
  out
}

B <- 9999L
species_cache <- file.path(out6, "tables", "R1_6_species_profile_global_gate.csv")
if (file.exists(species_cache)) {
  species_profiles <- read_csv(species_cache, show_col_types = FALSE)
} else {
  observed_stat <- omnibus_stat(group_by_subject)
  perm_exceed <- integer(length(observed_stat))
  for (b in seq_len(B)) {
    perm_stat <- omnibus_stat(sample(group_by_subject, replace = FALSE))
    perm_exceed <- perm_exceed + as.integer(perm_stat >= observed_stat - 1e-12)
  }
  omnibus_p <- (perm_exceed + 1) / (B + 1)
  omnibus_q <- p.adjust(omnibus_p, method = "BH")

  phase_effect_rows <- list()
  for (i in seq_along(available_candidates)) {
    for (j in seq_along(phase_order)) {
      h <- phase_array[i, group_by_subject == "HDR", j]
      l <- phase_array[i, group_by_subject == "LDR", j]
      diff <- mean(h) - mean(l)
      se <- sqrt(var(h) / length(h) + var(l) / length(l))
      df <- (var(h) / length(h) + var(l) / length(l))^2 /
        ((var(h) / length(h))^2 / (length(h) - 1) + (var(l) / length(l))^2 / (length(l) - 1))
      pooled <- sqrt(((length(h) - 1) * var(h) + (length(l) - 1) * var(l)) / (length(h) + length(l) - 2))
      d <- ifelse(pooled > 0, diff / pooled, 0)
      correction <- 1 - 3 / (4 * (length(h) + length(l)) - 9)
      phase_effect_rows[[length(phase_effect_rows) + 1L]] <- data.frame(
        species = available_candidates[i], phase = phase_order[j],
        HDR_minus_LDR_log10_CPM = diff, ci_low = diff - qt(0.975, df) * se,
        ci_high = diff + qt(0.975, df) * se, hedges_g = d * correction,
        omnibus_statistic = observed_stat[i], omnibus_p = omnibus_p[i], omnibus_q = omnibus_q[i]
      )
    }
  }
  species_profiles <- bind_rows(phase_effect_rows) |>
    mutate(global_pass = omnibus_q < 0.05,
           max_abs_g = ave(abs(hedges_g), species, FUN = max))
  write_csv(species_profiles, species_cache)
}

sample_metrics <- fread(sample_metric_file)[recovery_group %in% c("HDR", "LDR")]
sample_metrics[, canon := factor(canon, levels = time_order)]
trajectory_summary <- sample_metrics[, .(
  n_subjects = uniqueN(subject),
  shannon_mean = mean(shannon), shannon_se = sd(shannon) / sqrt(.N),
  bray_mean = mean(bray_to_baseline_centroid), bray_se = sd(bray_to_baseline_centroid) / sqrt(.N),
  jaccard_mean = mean(mean_jaccard_to_baseline), jaccard_se = sd(mean_jaccard_to_baseline) / sqrt(.N)
), by = .(recovery_group, canon)]
write_csv(trajectory_summary, file.path(out6, "tables", "R1_6_trajectory_summary.csv"))

gee_omnibus <- read_csv(file.path(root, "reviewer_response", "R1.6", "tables", "R1_6_GEE_omnibus.csv"), show_col_types = FALSE)
distance_tests <- read_csv(file.path(root, "reviewer_response", "R1.6", "tables", "R1_6_PERMANOVA_PERMDISP.csv"), show_col_types = FALSE)
write_csv(gee_omnibus, file.path(out6, "tables", "R1_6_repeated_measure_omnibus.csv"))
write_csv(distance_tests, file.path(out6, "tables", "R1_6_PERMANOVA_PERMDISP.csv"))

message("R1.7: focal-exclusion, CLR, and longitudinal within/between decomposition...")
genus_totals <- colSums(genus_mat)
genus_rel <- sweep(genus_mat, 2, genus_totals, "/")
prevalence <- rowMeans(genus_mat > 0)
mean_rel <- rowMeans(genus_rel)
focal_genera <- rownames(genus_mat)[prevalence >= 0.10 & mean_rel >= 0.0001 & rownames(genus_mat) != "Unclassified"]
audited_genus_file <- file.path(root, "reviewer_response", "R1.7", "tables", "R1_7_CLR_focal_excluded_correlations.csv")
audited_genera <- fread(audited_genus_file)[dataset_branch == "rarefied_primary", unique(genus)]
focal_genera <- intersect(focal_genera, audited_genera)
genus_log <- log(genus_mat + 0.5)
genus_clr <- sweep(genus_log, 2, colMeans(genus_log), "-")

fit_rows <- list()
for (g in focal_genera) {
  remain <- sample_totals - genus_mat[g, ]
  y <- log(remain) - (total_clogc - genus_clogc[g, ]) / remain
  x <- genus_clr[g, ]
  dat <- data.frame(
    y = y, x = x, subject = factor(meta_dt$subject), canon = factor(meta_dt$canon, levels = time_order),
    log_depth = meta_dt$log_depth
  )
  dat$x_between <- ave(dat$x, dat$subject, FUN = mean)
  dat$x_within <- dat$x - dat$x_between
  dat$y_z <- as.numeric(scale(dat$y))
  dat$x_between_z <- as.numeric(scale(dat$x_between))
  dat$x_within_z <- as.numeric(scale(dat$x_within))
  dat$log_depth_z <- as.numeric(scale(dat$log_depth))
  fit <- try(lme(y_z ~ x_between_z + x_within_z + canon + log_depth_z,
                 random = ~1 | subject, data = dat, method = "REML",
                 control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)), silent = TRUE)
  if (inherits(fit, "try-error")) next
  tt <- summary(fit)$tTable
  for (term in c("x_between_z", "x_within_z")) {
    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      genus = g,
      component = ifelse(term == "x_between_z", "Between-participant", "Within-participant"),
      beta_standardized = tt[term, "Value"], se = tt[term, "Std.Error"],
      df = tt[term, "DF"], p_value = tt[term, "p-value"],
      ci_low = tt[term, "Value"] - qt(0.975, tt[term, "DF"]) * tt[term, "Std.Error"],
      ci_high = tt[term, "Value"] + qt(0.975, tt[term, "DF"]) * tt[term, "Std.Error"],
      n_subjects = uniqueN(dat$subject), n_observations = nrow(dat)
    )
  }
}
longitudinal_assoc <- bind_rows(fit_rows) |>
  group_by(component) |>
  mutate(q_BH_across_genera = p.adjust(p_value, method = "BH")) |>
  ungroup()
write_csv(longitudinal_assoc, file.path(out7, "tables", "R1_7_longitudinal_CLR_focal_excluded_models.csv"))

per_time <- read_csv(file.path(root, "reviewer_response", "R1.7", "tables", "R1_7_CLR_focal_excluded_correlations.csv"), show_col_types = FALSE) |>
  filter(dataset_branch == "rarefied_primary")
write_csv(per_time, file.path(out7, "tables", "R1_7_time_resolved_sensitivity.csv"))

pattern_counts <- read_csv(file.path(root, "reviewer_response", "R1.7", "tables", "R1_7_corrected_phase_pattern_counts.csv"), show_col_types = FALSE)
upset_source <- read_csv(file.path(root, "reviewer_response", "R1.7", "tables", "R1_7_original_UpSet_source_counts_audit.csv"), show_col_types = FALSE)
write_csv(pattern_counts, file.path(out7, "tables", "R1_7_filtered_genus_pattern_counts.csv"))
write_csv(upset_source, file.path(out7, "tables", "R1_7_original_feature_count_audit.csv"))

# ---- Figures ----
pal <- c(HDR = "#2B6CB0", LDR = "#C05621")
base_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.8, colour = "grey30"),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(), plot.margin = margin(5, 7, 5, 5))

traj <- as.data.frame(trajectory_summary)
p6a <- ggplot(traj, aes(as.numeric(canon), shannon_mean, colour = recovery_group, fill = recovery_group)) +
  geom_ribbon(aes(ymin = shannon_mean - 1.96 * shannon_se, ymax = shannon_mean + 1.96 * shannon_se), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.7) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = seq_along(time_order), labels = time_order) +
  labs(title = "A  Diversity trajectories", subtitle = "Participant means and 95% CI; GEE group x time P = 3.58e-15",
       x = "Sampling point", y = "Shannon diversity") + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p6b <- ggplot(traj, aes(as.numeric(canon), bray_mean, colour = recovery_group, fill = recovery_group)) +
  geom_ribbon(aes(ymin = bray_mean - 1.96 * bray_se, ymax = bray_mean + 1.96 * bray_se), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.7) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = seq_along(time_order), labels = time_order) +
  labs(title = "B  Return toward personal baseline", subtitle = "Bray-Curtis distance to each participant's baseline centroid",
       x = "Sampling point", y = "Distance to personal baseline") + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

species_global <- species_profiles |>
  distinct(species, omnibus_p, omnibus_q, omnibus_statistic, max_abs_g) |>
  arrange(omnibus_q)
top_labels <- species_global |> slice_head(n = 12)
p6c <- ggplot(species_global, aes(max_abs_g, -log10(pmax(omnibus_q, 1e-12)), colour = omnibus_q < 0.05)) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, colour = "grey55") +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_text(data = top_labels, aes(label = species), size = 2.2, check_overlap = TRUE, vjust = -0.5) +
  scale_colour_manual(values = c(`TRUE` = "#7B3294", `FALSE` = "grey65"),
                      labels = c(`TRUE` = "Global q < 0.05", `FALSE` = "Other")) +
  labs(title = "C  Global-gated species profiles", subtitle = "One participant-label permutation test/species; BH across 268 species",
       x = "Maximum |Hedges' g| across three phases", y = expression(-log[10](q))) + base_theme

top_species <- species_global |> slice_head(n = 15) |> pull(species)
heat <- species_profiles |> filter(species %in% top_species) |>
  mutate(species = factor(species, levels = rev(top_species)), phase = factor(phase, levels = phase_order))
p6d <- ggplot(heat, aes(phase, species, fill = hedges_g)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(global_pass, "*", "")), size = 3.0, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(title = "D  Phase-specific effect sizes", subtitle = "Top-ranked profiles; * global q < 0.05; HDR minus LDR",
       x = NULL, y = NULL, fill = "Hedges' g") + base_theme +
  theme(axis.text.y = element_text(size = 5.8), legend.position = "right")

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(title = "R1.6 revised | Longitudinal inference aligned with the study design",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(out6, "figures", "R1_6_revised_response_figure.pdf"), fig6, width = 210, height = 180, units = "mm", device = "pdf", useDingbats = FALSE)
ggsave(file.path(out6, "figures", "R1_6_revised_response_figure.tif"), fig6, width = 210, height = 180, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", background = "white")
ggsave(file.path(out6, "figures", "R1_6_revised_response_figure_preview.png"), fig6, width = 210, height = 180, units = "mm", dpi = 180,
       device = ragg::agg_png, background = "white")

focus_genera <- c("Bacteroides", "Phocaeicola", "Faecalibacterium", "Blautia", "Alistipes", "Akkermansia")
forest <- longitudinal_assoc |> filter(genus %in% focus_genera) |>
  mutate(genus = factor(genus, levels = rev(focus_genera)))
p7a <- ggplot(forest, aes(beta_standardized, genus, colour = component, shape = component)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.13, position = position_dodge(width = 0.45)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Between-participant" = "#7B3294", "Within-participant" = "#008837")) +
  labs(title = "A  Longitudinal models", subtitle = "CLR + focal exclusion; time/depth adjusted",
       x = "Standardized association with focal-excluded Shannon", y = NULL) + base_theme

focus_time <- per_time |> filter(genus %in% c("Bacteroides", "Faecalibacterium")) |>
  mutate(canon = factor(canon, levels = time_order), significant = q_BH_global_16_times < 0.05)
p7b <- ggplot(focus_time, aes(as.numeric(canon), rho_CLR_vs_focal_genus_excluded_shannon, colour = genus, group = genus)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_line(linewidth = 0.7) +
  geom_point(aes(shape = significant, size = significant), stroke = 0.5) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) + scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.8)) +
  scale_colour_manual(values = c(Bacteroides = "#D55E00", Faecalibacterium = "#0072B2")) +
  scale_x_continuous(breaks = seq_along(time_order), labels = time_order) +
  labs(title = "B  Time-resolved associations", subtitle = "Filled: global BH q < 0.05",
       x = "Sampling point", y = expression(Spearman~rho)) + base_theme +
  guides(shape = "none", size = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

c6 <- per_time |> filter(canon == "C6") |>
  mutate(significant = q_BH_global_16_times < 0.05)
p7c <- ggplot(c6, aes(rho_relative_abundance_vs_unmodified_shannon, rho_CLR_vs_focal_genus_excluded_shannon, colour = significant)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_text(data = c6 |> filter(genus %in% focus_genera), aes(label = genus), size = 2.2, check_overlap = TRUE, vjust = -0.5) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65")) +
  coord_equal(xlim = c(-0.85, 0.9), ylim = c(-0.85, 0.9)) +
  labs(title = "C  Effect of correction at C6", subtitle = "C6 focal-exclusion sensitivity",
       x = "Relative abundance vs original Shannon", y = "CLR abundance vs focal-excluded Shannon") + base_theme +
  theme(legend.position = "none")

pc <- pattern_counts |> mutate(pattern = factor(pattern, levels = c("A", "B", "C", "AB", "AC", "BC", "ABC", "None")))
p7d <- ggplot(pc, aes(pattern, n_genera)) +
  geom_col(fill = "#009E73", width = 0.68) + geom_text(aes(label = n_genera), vjust = -0.25, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "D  Phase-intersection sensitivity", subtitle = "124 filtered named genera",
       x = "Mutually exclusive phase pattern", y = "Number of genera") + base_theme

fig7 <- (p7a | p7b) / (p7c | p7d) +
  plot_layout(widths = c(1.05, 1)) +
  plot_annotation(title = "R1.7 revised | Associations beyond same-vector mathematical coupling",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(out7, "figures", "R1_7_revised_response_figure.pdf"), fig7, width = 240, height = 185, units = "mm", device = "pdf", useDingbats = FALSE)
ggsave(file.path(out7, "figures", "R1_7_revised_response_figure.tif"), fig7, width = 240, height = 185, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", background = "white")
ggsave(file.path(out7, "figures", "R1_7_revised_response_figure_preview.png"), fig7, width = 240, height = 185, units = "mm", dpi = 180,
       device = ragg::agg_png, background = "white")

qa6 <- data.frame(
  check = c("60 exposed participants", "30 reference-group participants", "268 candidate species", "9999 participant-label permutations",
            "one omnibus test per candidate species", "PDF exists", "TIFF exists"),
  passed = c(uniqueN(meta_dt$subject) == 60L, nrow(assign) == 30L, length(candidates) == 268L, B == 9999L,
             n_distinct(species_profiles$species) == length(available_candidates),
             file.exists(file.path(out6, "figures", "R1_6_revised_response_figure.pdf")),
             file.exists(file.path(out6, "figures", "R1_6_revised_response_figure.tif")))
)
qa7 <- data.frame(
  check = c("60 exposed participants", "16 repeated visits", "focal genus universe matches audit", "within/between models available",
            "time-resolved sensitivity available", "PDF exists", "TIFF exists"),
  passed = c(uniqueN(meta_dt$subject) == 60L, uniqueN(meta_dt$canon) == 16L, length(focal_genera) == 124L,
             nrow(longitudinal_assoc) == 2L * length(focal_genera), nrow(per_time) == 124L * 16L,
             file.exists(file.path(out7, "figures", "R1_7_revised_response_figure.pdf")),
             file.exists(file.path(out7, "figures", "R1_7_revised_response_figure.tif")))
)
write_tsv(qa6, file.path(out6, "QA", "R1_6_revised_QA.tsv"))
write_tsv(qa7, file.path(out7, "QA", "R1_7_revised_QA.tsv"))
if (!all(qa6$passed) || !all(qa7$passed)) stop("QA failed")

summary6 <- species_global |> summarise(
  n_species = n(), n_global_q05 = sum(omnibus_q < 0.05),
  minimum_q = min(omnibus_q), median_max_abs_g = median(max_abs_g)
)
write_csv(summary6, file.path(out6, "tables", "R1_6_revised_key_results.csv"))
summary7 <- longitudinal_assoc |> filter(genus %in% focus_genera) |> arrange(component, q_BH_across_genera)
write_csv(summary7, file.path(out7, "tables", "R1_7_revised_key_results.csv"))

writeLines(capture.output(sessionInfo()), file.path(out6, "logs", "R_sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), file.path(out7, "logs", "R_sessionInfo.txt"))
message("R1.6/R1.7 revised analyses complete")
