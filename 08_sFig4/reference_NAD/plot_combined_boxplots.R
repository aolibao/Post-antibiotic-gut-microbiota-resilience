
library(ggplot2)
library(dplyr)
library(ggsignif)
library(gridExtra)

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

# Target Pathways
targets <- c(
  "PWY-7761: NAD salvage pathway II (PNC IV cycle)",
  "PYRIDNUCSYN-PWY: NAD de novo biosynthesis I (from aspartate)"
)

# Load Metadata
meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))
meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))

# Plot Function
create_plot <- function(meta_df, group_col, group1, group2, pathway_data, pathway_name, comparison_name) {
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group1, group2), ]
  common_samples <- meta_sub$sample_id
  
  if (length(common_samples) == 0) return(NULL)
  
  # Check if pathway exists
  if (!pathway_name %in% rownames(path_sub)) return(NULL)
  
  val <- as.numeric(path_sub[pathway_name, common_samples])
  plot_data <- data.frame(
    Sample = common_samples,
    Group = meta_sub[[group_col]],
    Abundance = val
  )
  
  # Wilcoxon Test
  p_val <- wilcox.test(Abundance ~ Group, data = plot_data)$p.value
  
  # Shorten pathway name for title if too long
  short_name <- ifelse(nchar(pathway_name) > 40, paste0(substr(pathway_name, 1, 37), "..."), pathway_name)
  
  p <- ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste0(short_name, "\n(", comparison_name, ")"),
      y = "Abundance (CPM)",
      x = NULL
    ) +
    theme(legend.position = "none", plot.title = element_text(size = 10)) +
    geom_signif(
      comparisons = list(c(group1, group2)),
      map_signif_level = TRUE,
      annotations = paste0("p = ", format(p_val, digits = 3))
    )
  
  return(p)
}

# Generate List of Plots
plot_list <- list()

for (pwy in targets) {
  # Comparison 1: P1 vs P4
  p1 <- create_plot(meta_end, "Group", "P1", "P4", pathway_df, pwy, "P1 vs P4")
  if (!is.null(p1)) plot_list[[paste0(pwy, "_P1_P4")]] <- p1
  
  # Comparison 2: Recover vs Non_recover
  p2 <- create_plot(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, pwy, "Rec vs Non")
  if (!is.null(p2)) plot_list[[paste0(pwy, "_Rec_Non")]] <- p2
}

# Arrange in one row
# We have 2 pathways * 2 comparisons = 4 plots max.
# User asked for "one row".
if (length(plot_list) > 0) {
  final_plot <- do.call(grid.arrange, c(plot_list, nrow = 1))
  ggsave(file.path(output_dir, "combined_boxplots_target_pathways.pdf"), final_plot, width = 16, height = 5)
  print("Saved combined_boxplots_target_pathways.pdf")
} else {
  print("No plots generated.")
}
