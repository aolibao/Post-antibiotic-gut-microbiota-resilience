#!/home/qqt/anaconda3/bin/Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)), mustWork = TRUE)
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
table_dir <- file.path(base_dir, "tables"); figure_dir <- file.path(base_dir, "figures"); log_dir <- file.path(base_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261013)

input <- file.path(root_dir, "ABX课题", "Analysis_Results", "Diversity_Tables", "alpha_diversity_selected_days.csv")
d <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
names(d)[1] <- "sample_id"
mapping <- c(A = "ABX-Trp", B = "ABX-Trp-During", C = "ABX-NMN-During", D = "ABX-NMN", E = "ABX-Con")
d$group_code <- substr(d$sample_id, 1, 1)
d$treatment <- unname(mapping[d$group_code])
d$animal_number <- sub("^[A-E](?:Neg)?([0-9]+).*", "\\1", d$sample_id, perl = TRUE)
d$animal_id <- paste(d$group_code, d$animal_number, sep = "-")
d <- d[d$days %in% c(-10, 7, 14, 28) & !is.na(d$treatment), ]
d$treatment <- factor(d$treatment, levels = c("ABX-Con", "ABX-Trp", "ABX-Trp-During", "ABX-NMN", "ABX-NMN-During"))
stopifnot(!anyDuplicated(d$sample_id))

desc_fun <- function(x) c(n = length(x), mean = mean(x), sd = sd(x), median = median(x), q1 = quantile(x, .25), q3 = quantile(x, .75))
desc <- do.call(rbind, lapply(split(d, list(d$days, d$treatment), drop = TRUE), function(x) {
  a <- desc_fun(x$Shannon); b <- desc_fun(x$Evenness)
  data.frame(day = unique(x$days), treatment = as.character(unique(x$treatment)),
             Shannon_n = a["n"], Shannon_mean = a["mean"], Shannon_sd = a["sd"], Shannon_median = a["median"],
             Evenness_n = b["n"], Evenness_mean = b["mean"], Evenness_sd = b["sd"], Evenness_median = b["median"])
}))
write.csv(desc, file.path(table_dir, "mouse_alpha_descriptives.csv"), row.names = FALSE)

hedges_g <- function(a, b) {
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  J * (mean(a) - mean(b)) / sp
}

welch_contrast <- function(a, b, outcome, day, contrast, family) {
  tt <- t.test(a, b, var.equal = FALSE)
  data.frame(outcome = outcome, day = day, contrast = contrast, family = family,
             n1 = length(a), n2 = length(b), mean1 = mean(a), mean2 = mean(b),
             mean_difference = mean(a) - mean(b), CI95_low = tt$conf.int[1], CI95_high = tt$conf.int[2],
             statistic = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value,
             hedges_g = hedges_g(a, b))
}

# Reviewer-requested day-28 evenness family: four treatment arms versus ABX-Con.
day28 <- d[d$days == 28, ]
control <- day28$Evenness[day28$treatment == "ABX-Con"]
day28_tests <- do.call(rbind, lapply(levels(d$treatment)[-1], function(g) {
  a <- day28$Evenness[day28$treatment == g]
  welch_contrast(a, control, "Evenness", 28, paste(g, "vs ABX-Con"), "Day-28 evenness: four treatment-vs-control contrasts")
}))
day28_tests$q_BH <- p.adjust(day28_tests$p_value, method = "BH")
write.csv(day28_tests, file.path(table_dir, "mouse_day28_evenness_treatment_vs_control_BH.csv"), row.names = FALSE)

# Direct timing family requested by the reviewer: pretreatment vs during-exposure
# tryptophan at the three recovery observations available in the source table.
timing_tests <- do.call(rbind, lapply(c(7, 14, 28), function(day) {
  x <- d[d$days == day, ]
  a <- x$Shannon[x$treatment == "ABX-Trp"]
  b <- x$Shannon[x$treatment == "ABX-Trp-During"]
  welch_contrast(a, b, "Shannon", day, "ABX-Trp vs ABX-Trp-During", "Direct timing: Shannon at days 7, 14 and 28")
}))
timing_tests$q_BH <- p.adjust(timing_tests$p_value, method = "BH")
write.csv(timing_tests, file.path(table_dir, "mouse_trp_pre_vs_during_direct_BH.csv"), row.names = FALSE)

# Supportive day-14 intervention contrasts versus control, corrected as one family.
day14 <- d[d$days == 14, ]
control14 <- day14$Shannon[day14$treatment == "ABX-Con"]
day14_tests <- do.call(rbind, lapply(levels(d$treatment)[-1], function(g) {
  a <- day14$Shannon[day14$treatment == g]
  welch_contrast(a, control14, "Shannon", 14, paste(g, "vs ABX-Con"), "Day-14 Shannon: four treatment-vs-control contrasts")
}))
day14_tests$q_BH <- p.adjust(day14_tests$p_value, method = "BH")
write.csv(day14_tests, file.path(table_dir, "mouse_day14_shannon_treatment_vs_control_BH.csv"), row.names = FALSE)

cols <- c("ABX-Con" = "#6B7280", "ABX-Trp" = "#B5524E", "ABX-Trp-During" = "#D28E5A",
          "ABX-NMN" = "#356A9A", "ABX-NMN-During" = "#63A6C7")
theme_set(theme_classic(base_size = 9))

p1 <- ggplot(d[d$days >= 7, ], aes(days, Shannon, color = treatment, group = treatment)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2.2) +
  stat_summary(fun.data = function(x) {
    m <- mean(x); se <- sd(x) / sqrt(length(x)); crit <- qt(.975, df = length(x) - 1)
    data.frame(y = m, ymin = m - crit * se, ymax = m + crit * se)
  }, geom = "errorbar", width = 0.7, linewidth = 0.45) +
  geom_point(alpha = 0.25, position = position_jitter(width = 0.35, height = 0), size = 1) +
  scale_color_manual(values = cols) + scale_x_continuous(breaks = c(7,14,28)) +
  labs(x = "Recovery day", y = "Shannon diversity", color = NULL,
       subtitle = "Animal-level values and means with 95% CI") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.2),
        plot.subtitle = element_text(size = 7.3))

p2 <- ggplot(day14, aes(treatment, Shannon, color = treatment)) +
  geom_boxplot(width = .55, outlier.shape = NA, color = "#555555", fill = "white") +
  geom_jitter(width = .11, size = 2) + scale_color_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "Day-14 Shannon", subtitle = "Raw animal-level observations") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

p3 <- ggplot(timing_tests, aes(mean_difference, factor(day))) +
  geom_vline(xintercept = 0, linetype = 2, color = "#888888") +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = .15, color = "#666666") +
  geom_point(size = 2.6, color = cols["ABX-Trp"]) +
  labs(x = "Shannon difference: pretreatment - during", y = "Recovery day",
       subtitle = "Timing contrasts (positive = pretreatment)") +
  theme(plot.subtitle = element_text(size = 7.3))

day28_tests$label <- sub(" vs ABX-Con", "", day28_tests$contrast)
p4 <- ggplot(day28_tests, aes(mean_difference, reorder(label, mean_difference))) +
  geom_vline(xintercept = 0, linetype = 2, color = "#888888") +
  geom_errorbar(aes(xmin = CI95_low, xmax = CI95_high), orientation = "y", width = .15, color = "#666666") +
  geom_point(size = 2.5, color = "#356A9A") +
  labs(x = "Day-28 evenness difference vs control", y = NULL,
       subtitle = "Day-28 contrasts; BH-corrected source table") +
  theme(plot.subtitle = element_text(size = 7.3))

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.07, 0.465, 0.465), "npc"))))
  grid.text("Animal-level reanalysis of mouse recovery outcomes",
            vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
            gp = gpar(fontface = "bold", fontsize = 13))
  print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
  print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2), newpage = FALSE)
  print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1), newpage = FALSE)
  print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 2), newpage = FALSE)
  grid.text("A", x = unit(0.012, "npc"), y = unit(0.915, "npc"),
            just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
  grid.text("B", x = unit(0.505, "npc"), y = unit(0.915, "npc"),
            just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
  grid.text("C", x = unit(0.012, "npc"), y = unit(0.455, "npc"),
            just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
  grid.text("D", x = unit(0.505, "npc"), y = unit(0.455, "npc"),
            just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 11))
}
pdf(file.path(figure_dir, "R1_9_mouse_alpha_reanalysis.pdf"), width = 225/25.4, height = 180/25.4, useDingbats = FALSE); draw_combined(); dev.off()
tiff(file.path(figure_dir, "R1_9_mouse_alpha_reanalysis.tif"), width = 225, height = 180, units = "mm", res = 600, compression = "lzw"); draw_combined(); dev.off()
png(file.path(figure_dir, "R1_9_mouse_alpha_reanalysis_word.png"), width = 225, height = 180, units = "mm", res = 180); draw_combined(); dev.off()

writeLines(c(
  paste0("Input: ", input),
  paste0("Rows: ", nrow(d), "; unique sample IDs: ", length(unique(d$sample_id))),
  "Treatment mapping recovered from supplied mouse analysis scripts: A=ABX-Trp, B=ABX-Trp-During, C=ABX-NMN-During, D=ABX-NMN, E=ABX-Con.",
  paste0("Direct timing estimates favoring pretreatment at days 7/14/28: ", paste(round(timing_tests$mean_difference, 4), collapse = ", ")),
  paste0("Direct timing BH q-values: ", paste(signif(timing_tests$q_BH, 4), collapse = ", ")),
  "No cage-allocation, randomisation/blinding, ARRIVE checklist, prespecified sample-size rationale, or mouse KO-by-animal table was found; those items remain author input.",
  paste0("R: ", R.version.string)
), file.path(log_dir, "R1_9_mouse_run_summary.txt"))
capture.output(sessionInfo(), file = file.path(log_dir, "R1_9_mouse_sessionInfo.txt"))
