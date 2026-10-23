
library(ggplot2)
library(dplyr)
library(ggsignif)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")

# Load Data
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
meta_endpoint_path <- file.path(base_dir, "species/kraken.P1.P4/endpoint_sample_metadata.csv")
meta_dynamic_path <- file.path(base_dir, "species/kraken.P1.P4/dynamic_sample_metadata.csv")

pathway_raw <- read.csv(pathway_path, check.names = FALSE)
colnames(pathway_raw)[1] <- "Pathway"
if (colnames(pathway_raw)[2] == "") pathway_raw <- pathway_raw[, -2]
if (any(duplicated(pathway_raw$Pathway))) {
  pathway_raw <- pathway_raw %>% group_by(Pathway) %>% summarise(across(where(is.numeric), sum))
}
pathway_df <- as.data.frame(pathway_raw)
rownames(pathway_df) <- pathway_df$Pathway
pathway_df$Pathway <- NULL

# Target Pathway
target_pwy <- "PWY-7761: NAD salvage pathway II (PNC IV cycle)"

# Check if pathway exists
if (!target_pwy %in% rownames(pathway_df)) {
  stop(paste("Pathway not found:", target_pwy))
}

# Load Metadata
meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))
meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))

# Function to plot boxplot
plot_boxplot <- function(meta_df, group_col, group1, group2, pathway_data, title_suffix) {
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group1, group2), ]
  common_samples <- meta_sub$sample_id
  
  if (length(common_samples) == 0) return(NULL)
  
  # Get data for target pathway
  val <- as.numeric(path_sub[target_pwy, common_samples])
  plot_data <- data.frame(
    Sample = common_samples,
    Group = meta_sub[[group_col]],
    Abundance = val
  )
  
  # Wilcoxon Test
  p_val <- wilcox.test(Abundance ~ Group, data = plot_data)$p.value
  
  # Plot
  p <- ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste(target_pwy, "-", title_suffix),
      y = "Abundance (CPM)",
      x = "Group"
    ) +
    theme(legend.position = "none") +
    geom_signif(
      comparisons = list(c(group1, group2)),
      map_signif_level = TRUE,
      annotations = paste0("p = ", format(p_val, digits = 3))
    )
  
  return(p)
}

# Generate Plots
p1 <- plot_boxplot(meta_end, "Group", "P1", "P4", pathway_df, "P1 vs P4")
if (!is.null(p1)) {
  ggsave(file.path(output_dir, "boxplot_PWY-7761_P1_vs_P4.pdf"), p1, width = 5, height = 6)
  print("Saved boxplot_PWY-7761_P1_vs_P4.pdf")
}

p2 <- plot_boxplot(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, "Recover vs Non_recover")
if (!is.null(p2)) {
  ggsave(file.path(output_dir, "boxplot_PWY-7761_Recover_vs_Non.pdf"), p2, width = 5, height = 6)
  print("Saved boxplot_PWY-7761_Recover_vs_Non.pdf")
}
