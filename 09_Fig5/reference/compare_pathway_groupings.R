
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
meta_endpoint_path <- file.path(base_dir, "species/kraken.P1.P4/endpoint_sample_metadata.csv")
meta_dynamic_path <- file.path(base_dir, "species/kraken.P1.P4/dynamic_sample_metadata.csv")
output_dir <- file.path(base_dir, "result/Fig5")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data
# Pathway Data
pathway_raw <- read.csv(pathway_path, check.names = FALSE)
colnames(pathway_raw)[1] <- "Pathway"
if (colnames(pathway_raw)[2] == "") {
  pathway_raw <- pathway_raw[, -2]
}
if (any(duplicated(pathway_raw$Pathway))) {
  pathway_raw <- pathway_raw %>% group_by(Pathway) %>% summarise(across(where(is.numeric), sum))
}
pathway_df <- as.data.frame(pathway_raw)
rownames(pathway_df) <- pathway_df$Pathway
pathway_df$Pathway <- NULL

# Filter out UNMAPPED and UNINTEGRATED
ignore_pathways <- c("UNMAPPED", "UNINTEGRATED")
pathway_df <- pathway_df[!rownames(pathway_df) %in% ignore_pathways, ]
pathway_cols <- colnames(pathway_df)

# Load Metadata
meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))
# endpoint usually has sample_id, subject, Group(P1/P4), timepoint

meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))
# dynamic usually has sample_id, subject, Group(Recover/Non_recover), timepoint

# 3. Define Analysis Function
run_comparison <- function(meta_df, group_col, group1, group2, pathway_data, label) {
  print(paste("Running Comparison:", label, "-", group1, "vs", group2))
  
  # Align samples
  # Check if sample_id is in pathway columns
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  print(paste("Aligned Samples:", length(common_samples)))
  
  if (length(common_samples) == 0) return(NULL)
  
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  # Filter groups
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group1, group2), ]
  # Re-align after group filter
  common_samples <- meta_sub$sample_id
  path_sub <- path_sub[, common_samples]
  
  if (nrow(meta_sub) == 0) return(NULL)
  
  # Transpose
  path_t <- as.data.frame(t(path_sub))
  
  # Filter low abundance
  keep_feat <- colSums(path_t > 0) >= 3
  path_t <- path_t[, keep_feat, drop=FALSE]
  
  results <- data.frame(
    Pathway = colnames(path_t),
    P_Val = NA,
    Log2FC = NA,
    Mean_G1 = NA,
    Mean_G2 = NA
  )
  
  # Stats
  grp_vec <- meta_sub[[group_col]]
  
  for (i in 1:nrow(results)) {
    feat <- results$Pathway[i]
    val <- path_t[[feat]]
    
    # Wilcoxon
    p_val <- tryCatch({
      wilcox.test(val ~ grp_vec)$p.value
    }, error = function(e) NA)
    
    # Log2FC
    means <- tapply(val, grp_vec, mean)
    m1 <- if (!is.na(means[group1])) means[group1] else 0
    m2 <- if (!is.na(means[group2])) means[group2] else 0
    
    # Log2FC = log2(G1 / G2)
    # Add pseudocount
    log2fc <- log2((m1 + 1e-6) / (m2 + 1e-6))
    
    results$P_Val[i] <- p_val
    results$Log2FC[i] <- log2fc
    results$Mean_G1[i] <- m1
    results$Mean_G2[i] <- m2
  }
  
  # FDR
  results$P_Adj <- p.adjust(results$P_Val, method = "BH")
  results$Comparison <- label
  
  return(results)
}

# 4. Run Analyses
# Comparison 1: P1 vs P4 (Endpoint Meta)
# Usually P1 is resilient? Let's check Recover vs Non_recover mapping.
# Assuming P1 ~ Recover, P4 ~ Non_recover based on context (often P1 is better outcome).
# Let's verify mapping later.
res_p1_p4 <- run_comparison(meta_end, "Group", "P1", "P4", pathway_df, "P1_vs_P4")

# Comparison 2: Recover vs Non_recover (Dynamic Meta)
res_rec_non <- run_comparison(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, "Recover_vs_Non")

# 5. Merge and Plot
if (!is.null(res_p1_p4) && !is.null(res_rec_non)) {
  
  # Rename columns for merge
  df1 <- res_p1_p4[, c("Pathway", "Log2FC", "P_Adj")]
  colnames(df1) <- c("Pathway", "Log2FC_P1_P4", "Padj_P1_P4")
  
  df2 <- res_rec_non[, c("Pathway", "Log2FC", "P_Adj")]
  colnames(df2) <- c("Pathway", "Log2FC_Rec_Non", "Padj_Rec_Non")
  
  merged_df <- merge(df1, df2, by = "Pathway")
  
  # Add Significance Class
  merged_df$Sig <- "NS"
  merged_df$Sig[merged_df$Padj_P1_P4 < 0.05 & merged_df$Padj_Rec_Non < 0.05] <- "Both"
  merged_df$Sig[merged_df$Padj_P1_P4 < 0.05 & merged_df$Padj_Rec_Non >= 0.05] <- "P1_P4_Only"
  merged_df$Sig[merged_df$Padj_P1_P4 >= 0.05 & merged_df$Padj_Rec_Non < 0.05] <- "Rec_Non_Only"
  
  # Plot
  p <- ggplot(merged_df, aes(x = Log2FC_P1_P4, y = Log2FC_Rec_Non, color = Sig)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    theme_bw() +
    labs(
      x = "Log2FC (P1 / P4)",
      y = "Log2FC (Recover / Non_recover)",
      title = "Comparison of Pathway Enrichment"
    ) +
    scale_color_manual(values = c("Both" = "red", "P1_P4_Only" = "blue", "Rec_Non_Only" = "green", "NS" = "grey"))
  
  ggsave(file.path(output_dir, "scatter_pathway_log2fc_comparison.pdf"), p, width = 8, height = 6)
  
  # Correlation
  cor_val <- cor(merged_df$Log2FC_P1_P4, merged_df$Log2FC_Rec_Non, use = "complete.obs")
  print(paste("Correlation between Log2FCs:", cor_val))
  
  write.csv(merged_df, file.path(output_dir, "pathway_comparison_stats.csv"), row.names = FALSE)
  
} else {
  print("Failed to run one or both comparisons.")
}
