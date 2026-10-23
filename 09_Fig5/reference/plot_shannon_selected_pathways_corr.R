
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
species_path <- file.path(base_dir, "species/kranken.all/combine_brk_all.csv")

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

# Target Pathways
targets <- c(
  "PWY-7761: NAD salvage pathway II (PNC IV cycle)",
  "PWY-5265: peptidoglycan biosynthesis II (staphylococci)",
  "PYRIDNUCSYN-PWY: NAD de novo biosynthesis I (from aspartate)"
)

# Filter targets
# Use partial matching if exact not found, but prefer exact
pathway_subset <- pathway_df[rownames(pathway_df) %in% targets, ]

# If missing, try partial matching
if (nrow(pathway_subset) < length(targets)) {
  print("Warning: Some exact targets not found. Checking partial matches...")
  found_rows <- rownames(pathway_df)[rownames(pathway_df) %in% targets]
  missing <- setdiff(targets, found_rows)
  print(paste("Missing:", paste(missing, collapse=", ")))
}

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

# 4. Calculate Correlation Stats for Annotation
cor_stats <- plot_long %>%
  group_by(Pathway) %>%
  summarise(
    Rho = cor.test(Abundance, Shannon, method = "spearman")$estimate,
    P_Value = cor.test(Abundance, Shannon, method = "spearman")$p.value
  ) %>%
  mutate(
    Label = paste0("Rho = ", round(Rho, 3), "\nP = ", format.pval(P_Value, digits = 3))
  )

# 5. Plot
p <- ggplot(plot_long, aes(x = Shannon, y = Abundance)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title = "Correlation: Shannon Index vs Selected Pathways",
    x = "Shannon Diversity Index",
    y = "Pathway Abundance (CPM)"
  ) +
  # Add correlation stats as text
  geom_text(
    data = cor_stats,
    aes(x = min(plot_long$Shannon), y = Inf, label = Label),
    hjust = 0, vjust = 1.5, size = 4, inherit.aes = FALSE
  ) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(output_dir, "scatter_shannon_selected_pathways.pdf"), p, width = 12, height = 5)
print("Saved scatter_shannon_selected_pathways.pdf")
