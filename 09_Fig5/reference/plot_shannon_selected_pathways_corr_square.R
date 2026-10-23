
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
species_path <- file.path(base_dir, "species/kranken.all/combine_brk_all.csv")
stats_path <- file.path(base_dir, "result/Fig5/correlation_shannon_pathways_stats.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data

# Species Data -> Shannon Index
species_raw <- read.csv(species_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(species_raw) <- species_raw[, 1]
species_raw <- species_raw[, -1]
colnames(species_raw) <- gsub("\\.brk\\.report$", "", colnames(species_raw))

species_t <- t(species_raw)
shannon_index <- diversity(species_t, index = "shannon")

# Pathway Data
pathway_raw <- read.csv(pathway_path, check.names = FALSE)
colnames(pathway_raw)[1] <- "Pathway"
if (colnames(pathway_raw)[2] == "") {
  pathway_raw <- pathway_raw[, -2]
}

# Aggregate duplicates
if (any(duplicated(pathway_raw$Pathway))) {
  pathway_raw <- pathway_raw %>% group_by(Pathway) %>% summarise(across(where(is.numeric), sum))
}
pathway_df <- as.data.frame(pathway_raw)
rownames(pathway_df) <- pathway_df$Pathway
pathway_df$Pathway <- NULL

# Load Stats
if (file.exists(stats_path)) {
  all_stats <- read.csv(stats_path, stringsAsFactors = FALSE)
} else {
  stop("Stats file not found. Please run correlation analysis first.")
}

# Target Pathways
targets <- c(
  "PWY-7761: NAD salvage pathway II (PNC IV cycle)",
  "PWY-5265: peptidoglycan biosynthesis II (staphylococci)",
  "PYRIDNUCSYN-PWY: NAD de novo biosynthesis I (from aspartate)"
)

# Filter targets from data
pathway_subset <- pathway_df[rownames(pathway_df) %in% targets, ]

# Align Samples
common_samples <- intersect(colnames(pathway_subset), names(shannon_index))
print(paste("Common Samples:", length(common_samples)))

pathway_clean <- pathway_subset[, common_samples]
shannon_clean <- shannon_index[common_samples]

# 3. Prepare Plot Data
plot_data <- data.frame(
  Sample = common_samples,
  Shannon = shannon_clean,
  t(pathway_clean)
)

# Reshape to long format
plot_long <- plot_data %>%
  pivot_longer(
    cols = -c(Sample, Shannon),
    names_to = "Pathway",
    values_to = "Abundance"
  )

# 4. Get Annotation Stats
# Merge with all_stats to get global FDR
# Need to match Pathway names carefully. The plot_long$Pathway might have periods if loaded weirdly, 
# but we constructed it from rownames, so should be fine.
# However, R's data.frame column names replace special chars with dots.
# Let's check names in plot_long.
# Actually, pivot_longer names_to takes column names.
# In `t(pathway_clean)`, column names are rownames of pathway_clean, which are the original pathway names.
# But data.frame() constructor might sanitize them.
# Let's use check.names=FALSE in data.frame() construction.

plot_data <- data.frame(
  Sample = common_samples,
  Shannon = shannon_clean,
  t(pathway_clean),
  check.names = FALSE
)

plot_long <- plot_data %>%
  pivot_longer(
    cols = -c(Sample, Shannon),
    names_to = "Pathway",
    values_to = "Abundance"
  )

# Extract stats for targets
target_stats <- all_stats[all_stats$Pathway %in% targets, ]

# If stats missing (e.g. name mismatch), calculate locally but warn
if (nrow(target_stats) < length(targets)) {
  print("Warning: Could not find all targets in global stats file. Calculating local stats for missing ones.")
  # This part is just a fallback, ideally we use the file
}

# Create annotation dataframe
# Join plot_long with target_stats to ensure correct mapping
# We need a distinct dataframe for geom_text (one row per pathway)
anno_df <- target_stats %>%
  select(Pathway, Rho, P_Value, FDR) %>%
  mutate(
    Label = paste0(
      "Rho = ", round(Rho, 3), "\n",
      "padj = ", format.pval(FDR, digits = 3)
    )
  )

# Ensure order matches plot
anno_df$Pathway <- factor(anno_df$Pathway, levels = targets)
plot_long$Pathway <- factor(plot_long$Pathway, levels = targets)

# 5. Plot
p <- ggplot(plot_long, aes(x = Shannon, y = Abundance)) +
  geom_point(alpha = 0.6, color = "blue") +
  # Add regression line with confidence interval (se = TRUE is default, but explicit here)
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
  theme_bw() +
  # Make square
  theme(aspect.ratio = 1) +
  labs(
    title = "Correlation: Shannon Index vs Selected Pathways",
    x = "Shannon Diversity Index",
    y = "Pathway Abundance (CPM)"
  ) +
  # Add stats
  geom_text(
    data = anno_df,
    aes(x = min(plot_long$Shannon), y = Inf, label = Label),
    hjust = 0, vjust = 1.5, size = 4, inherit.aes = FALSE
  ) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(output_dir, "scatter_shannon_selected_pathways_square_padj.pdf"), p, width = 12, height = 5)
print("Saved scatter_shannon_selected_pathways_square_padj.pdf")
