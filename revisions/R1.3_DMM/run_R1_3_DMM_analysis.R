#!/home/qqt/anaconda3/bin/Rscript

suppressPackageStartupMessages({
  library(DirichletMultinomial)
  library(ggplot2)
  library(grid)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(base_dir, "tables")
figure_dir <- file.path(base_dir, "figures")
log_dir <- file.path(base_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

seed0 <- 20261013L
set.seed(seed0)

z <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

hc3 <- function(fit) {
  x <- model.matrix(fit)
  h <- hatvalues(fit)
  bread <- solve(crossprod(x))
  meat <- crossprod(x, x * (residuals(fit)^2 / (1 - h)^2))
  bread %*% meat %*% bread
}

hc3_table <- function(fit, model_id) {
  v <- hc3(fit); b <- coef(fit); se <- sqrt(diag(v)); df <- df.residual(fit)
  crit <- qt(0.975, df)
  data.frame(
    model_id = model_id, term = names(b), estimate = unname(b), HC3_SE = unname(se),
    statistic = unname(b / se), df_residual = df,
    p_value = 2 * pt(abs(b / se), df = df, lower.tail = FALSE),
    CI95_low = unname(b - crit * se), CI95_high = unname(b + crit * se),
    n_subjects = nobs(fit), R2 = summary(fit)$r.squared,
    adjusted_R2 = summary(fit)$adj.r.squared, check.names = FALSE
  )
}

block_test <- function(fit, pattern, label) {
  terms <- grep(pattern, names(coef(fit)), value = TRUE)
  if (!length(terms)) return(data.frame(block = label, df1 = 0L, df2 = df.residual(fit), F = NA, p_value = NA))
  b <- coef(fit)[terms]; v <- hc3(fit)[terms, terms, drop = FALSE]; q <- length(terms)
  f <- as.numeric(t(b) %*% solve(v, b) / q)
  data.frame(block = label, terms = paste(terms, collapse = "; "), df1 = q,
             df2 = df.residual(fit), F = f, p_value = pf(f, q, df.residual(fit), lower.tail = FALSE))
}

meta <- read.csv(file.path(root_dir, "raw", "meta.csv"), stringsAsFactors = FALSE,
                 colClasses = c(ids = "character", subject = "character"))
base_meta <- meta[meta$abx == "Y" & meta$canon %in% paste0("A", 1:4), ]
stopifnot(length(unique(base_meta$subject)) == 60L, nrow(base_meta) == 240L)

message("Reading Bracken MPA table ...")
mpa <- read.delim(file.path(root_dir, "raw", "combine_brk.mpa"), check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
tax <- mpa[[1]]
terminal <- sub(".*\\|", "", tax)
keep <- startsWith(tax, "k__Bacteria|") & startsWith(terminal, "s__")
mpa <- mpa[keep, , drop = FALSE]
tax <- tax[keep]
sample_names <- sub("\\.brk\\.report$", "", names(mpa)[-1])
names(mpa)[-1] <- sample_names
missing_samples <- setdiff(base_meta$ids, sample_names)
if (length(missing_samples)) stop("Baseline samples missing from Bracken table: ", paste(missing_samples, collapse = ", "))

genus <- vapply(strsplit(tax, "\\|", fixed = FALSE), function(parts) {
  hit <- parts[startsWith(parts, "g__")]
  if (!length(hit) || hit[1] == "g__") "Unclassified" else sub("^g__", "", hit[1])
}, character(1))

sample_counts <- as.matrix(mpa[, base_meta$ids, drop = FALSE])
storage.mode(sample_counts) <- "double"
genus_by_sample <- rowsum(sample_counts, group = genus, reorder = TRUE)
all_subject_by_genus <- rowsum(t(genus_by_sample), group = base_meta$subject, reorder = TRUE)
all_subject_by_genus <- all_subject_by_genus[, colSums(all_subject_by_genus) > 0, drop = FALSE]
storage.mode(all_subject_by_genus) <- "integer"
stopifnot(nrow(all_subject_by_genus) == 60L, all(rowSums(all_subject_by_genus) > 0))

write.csv(data.frame(subject = rownames(all_subject_by_genus), all_subject_by_genus, check.names = FALSE),
          file.path(table_dir, "R1_3_DMM_baseline_genus_counts.csv"), row.names = FALSE)

# DMM estimation excludes only near-zero components with cohort mean abundance
# below 0.01%; the full QC genus matrix is retained above and remains available
# for the separate all-feature ecological analyses requested in Comment 7.
all_rel <- all_subject_by_genus / rowSums(all_subject_by_genus)
dmm_keep <- colMeans(all_rel) >= 1e-4
dmm_rel <- all_rel[, dmm_keep, drop = FALSE]
# Standardize the multinomial library size to avoid numerical instability from
# the approximately 90-million-read aggregated baseline libraries while
# preserving each participant's genus composition.
subject_by_genus <- t(apply(dmm_rel, 1, function(v) as.integer(round(v / sum(v) * 100000))))
colnames(subject_by_genus) <- colnames(dmm_rel)
rownames(subject_by_genus) <- rownames(dmm_rel)
write.csv(data.frame(genus = colnames(all_subject_by_genus),
                     mean_relative_abundance = colMeans(all_rel),
                     prevalence_subjects = colSums(all_subject_by_genus > 0),
                     included_in_DMM = dmm_keep),
          file.path(table_dir, "R1_3_DMM_genus_inclusion_audit.csv"), row.names = FALSE)

# Multiple deterministic starts per k; the best Laplace fit is retained.
k_grid <- 2:6
n_starts <- 1L
fits <- vector("list", length(k_grid))
selection_rows <- list()
for (k in k_grid) {
  candidates <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    fit_seed <- seed0 + 100L * k + s
    message("Fitting DMM k=", k, ", start=", s, ", seed=", fit_seed)
    candidates[[s]] <- dmn(subject_by_genus, k = k, verbose = FALSE, seed = fit_seed)
    if (is.null(candidates[[s]])) stop("DMM returned NULL for k=", k, ", start=", s)
    selection_rows[[length(selection_rows) + 1L]] <- data.frame(
      k = k, start = s, seed = fit_seed, laplace = laplace(candidates[[s]]),
      AIC = AIC(candidates[[s]]), BIC = BIC(candidates[[s]])
    )
  }
  vals <- vapply(candidates, laplace, numeric(1))
  fits[[k]] <- candidates[[which.min(vals)]]
}
selection <- do.call(rbind, selection_rows)
best_by_k <- do.call(rbind, lapply(split(selection, selection$k), function(x) x[which.min(x$laplace), ]))
best_k <- best_by_k$k[which.min(best_by_k$laplace)]
best_fit <- fits[[best_k]]
write.csv(selection, file.path(table_dir, "R1_3_DMM_model_selection_all_starts.csv"), row.names = FALSE)
write.csv(best_by_k, file.path(table_dir, "R1_3_DMM_model_selection_best_by_k.csv"), row.names = FALSE)

post <- mixture(best_fit)
if (is.null(dim(post))) post <- matrix(post, ncol = 1)
cluster_num <- max.col(post, ties.method = "first")
cluster_conf <- apply(post, 1, max)

mix <- fitted(best_fit)
if (is.null(dim(mix))) mix <- matrix(mix, ncol = 1)
component_markers <- vapply(seq_len(ncol(mix)), function(j) {
  if (!is.null(rownames(mix))) rownames(mix)[which.max(mix[, j])]
  else colnames(subject_by_genus)[which.max(mix[, j])]
}, character(1))
cluster_label <- paste0("DMM", cluster_num, ":", component_markers[cluster_num])

metrics <- read.csv(file.path(root_dir, "reviewer_response", "03_tables", "core_subject_recovery_metrics.csv"), stringsAsFactors = FALSE)
metrics <- metrics[metrics$dataset_branch == "rarefied_primary", ]
idx <- match(rownames(subject_by_genus), metrics$subject)
stopifnot(!anyNA(idx))
d <- metrics[idx, ]
d$dmm_cluster <- factor(cluster_label)
d$dmm_component <- cluster_num
d$dmm_posterior <- cluster_conf
d$z_shannon <- z(d$baseline_shannon)
d$z_recovery <- z(-d$recovery_auc_bray_actual_day)
d$endpoint_group <- factor(d$original_endpoint_group, levels = c("LDR", "Middle", "HDR"))

assignments <- data.frame(subject = d$subject, dmm_component = cluster_num,
                          dmm_enterotype = cluster_label, posterior_probability = cluster_conf,
                          baseline_shannon = d$baseline_shannon,
                          recovery_score_D5 = -d$recovery_auc_bray_actual_day,
                          endpoint_group = d$original_endpoint_group)
write.csv(assignments, file.path(table_dir, "R1_3_DMM_subject_assignments.csv"), row.names = FALSE)

fit_shannon <- lm(z_recovery ~ z_shannon, data = d)
fit_dmm <- lm(z_recovery ~ dmm_cluster, data = d)
fit_both <- lm(z_recovery ~ z_shannon + dmm_cluster, data = d)
models <- rbind(hc3_table(fit_shannon, "Shannon only"),
                hc3_table(fit_dmm, "DMM only"),
                hc3_table(fit_both, "Shannon + DMM"))
write.csv(models, file.path(table_dir, "R1_3_DMM_conditional_models_HC3.csv"), row.names = FALSE)
blocks <- rbind(block_test(fit_dmm, "^dmm_cluster", "DMM enterotype"),
                block_test(fit_both, "^dmm_cluster", "DMM enterotype | baseline Shannon"))
write.csv(blocks, file.path(table_dir, "R1_3_DMM_block_tests_HC3.csv"), row.names = FALSE)

ext <- d[d$endpoint_group %in% c("LDR", "HDR"), , drop = FALSE]
ext$endpoint_group <- droplevels(ext$endpoint_group)
cross <- table(ext$endpoint_group, ext$dmm_cluster)
if (ncol(cross) == 1L) {
  fisher_p <- 1
  fisher_method <- "Single DMM component; HDR/LDR association is not defined and is reported as no stratification"
} else {
  ft <- fisher.test(cross, simulate.p.value = TRUE, B = 99999)
  fisher_p <- ft$p.value
  fisher_method <- ft$method
}
write.csv(as.data.frame(cross), file.path(table_dir, "R1_3_DMM_HDR_LDR_crossclassification.csv"), row.names = FALSE)
write.csv(data.frame(n_subjects = nrow(ext), DMM_components = ncol(cross),
                     test = fisher_method, p_value = fisher_p),
          file.path(table_dir, "R1_3_DMM_HDR_LDR_association_test.csv"), row.names = FALSE)

# Relative abundance profiles for the most abundant genera; all genera remain in the fitted DMM.
rel <- subject_by_genus / rowSums(subject_by_genus)
top <- names(sort(colMeans(rel), decreasing = TRUE))[seq_len(min(12L, ncol(rel)))]
profile <- aggregate(rel[, top, drop = FALSE], list(dmm_enterotype = d$dmm_cluster), mean)
write.csv(profile, file.path(table_dir, "R1_3_DMM_component_profiles_top_genera.csv"), row.names = FALSE)

# PCA is only a visualization; DMM fitting uses the full count matrix.
pseudocount <- 0.5 / median(rowSums(subject_by_genus))
clr <- log(rel + pseudocount)
clr <- clr - rowMeans(clr)
pca <- prcomp(clr, center = TRUE, scale. = FALSE)
pca_df <- data.frame(subject = d$subject, PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     dmm_enterotype = d$dmm_cluster, endpoint_group = d$endpoint_group)

cols <- c("#356A9A", "#B5524E", "#3E8C89", "#C7973A", "#7A5195", "#6B6B6B")
theme_set(theme_classic(base_size = 9))
p1 <- ggplot(best_by_k, aes(k, laplace)) + geom_line(color = "#777777") +
  geom_point(aes(color = k == best_k), size = 2.7) +
  scale_color_manual(values = c(`FALSE` = "#356A9A", `TRUE` = "#B5524E"), guide = "none") +
  scale_x_continuous(breaks = k_grid) + labs(x = "Number of DMM components", y = "Laplace approximation",
  subtitle = paste0("Minimum criterion selects k = ", best_k))

p2 <- ggplot(pca_df, aes(PC1, PC2, color = dmm_enterotype, shape = endpoint_group)) +
  geom_point(size = 2.3, alpha = 0.9) + scale_color_manual(values = cols) +
  labs(x = sprintf("CLR PC1 (%.1f%%)", 100 * summary(pca)$importance[2,1]),
       y = sprintf("CLR PC2 (%.1f%%)", 100 * summary(pca)$importance[2,2]), color = "DMM enterotype", shape = "Original group",
       subtitle = "PCA display only; DMM fitted to QC genera") +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 6.2), legend.title = element_text(size = 6.5),
        plot.subtitle = element_text(size = 7.3))

cross_df <- as.data.frame(cross)
p3 <- ggplot(cross_df, aes(Var1, Freq, fill = Var2)) + geom_col(position = "fill", width = 0.65) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), size = 3) +
  scale_fill_manual(values = cols) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Original endpoint-Shannon group", y = "Proportion", fill = "DMM enterotype",
       subtitle = paste0("n=15/group; association P = ", formatC(fisher_p, digits = 3, format = "f"))) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.7))

forest <- models[models$model_id == "Shannon + DMM" & models$term != "(Intercept)", ]
forest$label <- ifelse(forest$term == "z_shannon", "Baseline Shannon | DMM", sub("dmm_cluster", "", forest$term))
p4 <- ggplot(forest, aes(estimate, reorder(label, estimate))) + geom_vline(xintercept = 0, linetype = 2, color = "#888888") +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = 0.15, color = "#666666") +
  geom_point(size = 2.5, color = "#356A9A") + labs(x = "Standardized coefficient (HC3 95% CI)", y = NULL,
  subtitle = "Continuous recovery; participant-level inference") +
  theme(plot.subtitle = element_text(size = 7.3))

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.07, 0.465, 0.465), "npc"))))
  grid.text("DMM enterotype analysis of baseline community structure",
            vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
            gp = gpar(fontface = "bold", fontsize = 13))
  print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
  print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2), newpage = FALSE)
  print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1), newpage = FALSE)
  print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 2), newpage = FALSE)
}
pdf(file.path(figure_dir, "R1_3_DMM_enterotype_analysis.pdf"), width = 225/25.4, height = 180/25.4, useDingbats = FALSE); draw_combined(); dev.off()
tiff(file.path(figure_dir, "R1_3_DMM_enterotype_analysis.tif"), width = 225, height = 180, units = "mm", res = 600, compression = "lzw"); draw_combined(); dev.off()
png(file.path(figure_dir, "R1_3_DMM_enterotype_analysis_word.png"), width = 225, height = 180, units = "mm", res = 180); draw_combined(); dev.off()

writeLines(c(
  paste0("DirichletMultinomial version: ", packageVersion("DirichletMultinomial")),
  paste0("R: ", R.version.string),
  paste0("Subjects: ", nrow(subject_by_genus)),
  paste0("Baseline samples: ", nrow(base_meta)),
  paste0("All detected QC genus features retained in audit matrix: ", ncol(all_subject_by_genus)),
  paste0("DMM-model genera (cohort mean abundance >=0.01%): ", ncol(subject_by_genus)),
  paste0("DMM standardized multinomial library: 100,000 counts/participant (rounding range ",
         min(rowSums(subject_by_genus)), "-", max(rowSums(subject_by_genus)), ")"),
  paste0("Selected DMM components by minimum Laplace: ", best_k),
  paste0("HDR/LDR versus DMM association P: ", signif(fisher_p, 6)),
  paste0("Input: raw/combine_brk.mpa + raw/meta.csv"),
  paste0("Seed base: ", seed0)
), file.path(log_dir, "R1_3_DMM_run_summary.txt"))
capture.output(sessionInfo(), file = file.path(log_dir, "R1_3_DMM_sessionInfo.txt"))
message("DMM analysis complete. Selected k=", best_k)
