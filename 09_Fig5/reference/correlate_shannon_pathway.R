
library(vegan)
library(dplyr)
library(ggplot2)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
species_path <- file.path(base_dir, "species/kranken.all/combine_brk_all.csv")

# 2. Load Data
# Pathway Data
pathway_raw <- read.csv(pathway_path, check.names = FALSE)
colnames(pathway_raw)[1] <- "Pathway"
if (colnames(pathway_raw)[2] == "") pathway_raw <- pathway_raw[, -2]
if (any(duplicated(pathway_raw$Pathway))) {
  pathway_raw <- pathway_raw %>% group_by(Pathway) %>% summarise(across(where(is.numeric), sum))
}
pathway_df <- as.data.frame(pathway_raw)
rownames(pathway_df) <- pathway_df$Pathway
pathway_df$Pathway <- NULL

# Calculate Total Pathway Abundance (Sum of CPM per sample)
total_pathway_abundance <- colSums(pathway_df)

# Species Data
species_raw <- read.csv(species_path, check.names = FALSE, stringsAsFactors = FALSE)
# Assuming first column is ID or taxonomy
rownames(species_raw) <- species_raw[, 1]
species_raw <- species_raw[, -1]

# Clean column names to match pathway data
colnames(species_raw) <- gsub("\\.brk\\.report$", "", colnames(species_raw))

# Calculate Shannon Index
# Shannon index requires species as columns, samples as rows
species_t <- t(species_raw)
shannon_index <- diversity(species_t, index = "shannon")

# 3. Align Data
common_samples <- intersect(names(total_pathway_abundance), names(shannon_index))
print(paste("Aligned Samples:", length(common_samples)))

if (length(common_samples) == 0) stop("No common samples found.")

plot_data <- data.frame(
  Sample = common_samples,
  Shannon = shannon_index[common_samples],
  TotalPathway = total_pathway_abundance[common_samples]
)

# 4. Correlation Analysis
cor_res <- cor.test(plot_data$Shannon, plot_data$TotalPathway, method = "spearman")
print(cor_res)

p_val_fmt <- format.pval(cor_res$p.value, digits = 3)
rho_val <- round(cor_res$estimate, 3)
label_text <- paste0("Spearman rho = ", rho_val, "\n", "p = ", p_val_fmt)

# 5. Plot
p <- ggplot(plot_data, aes(x = Shannon, y = TotalPathway)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  theme_bw() +
  labs(
    title = "Correlation: Shannon Index vs Total Pathway Abundance",
    x = "Shannon Diversity Index",
    y = "Total Pathway Abundance (CPM)"
  ) +
  annotate("text", x = min(plot_data$Shannon), y = max(plot_data$TotalPathway), label = label_text, hjust = 0, vjust = 1, size = 5)

ggsave(file.path(output_dir, "scatter_shannon_vs_pathway_abundance.pdf"), p, width = 6, height = 5)
print("Saved scatter_shannon_vs_pathway_abundance.pdf")
