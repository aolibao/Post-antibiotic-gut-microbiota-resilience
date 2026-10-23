
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")

# Load Results (Assuming we use the output from compare_pathway_lefse.R but we need ALL data for volcano, not just significant ones)
# The previous script `compare_pathway_lefse.R` only returned significant results.
# We need to re-run the analysis to get ALL results, or modify the script to return everything.
# Let's modify the workflow to calculate stats for all pathways first, then filter for plotting/labeling.

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

# Filter out UNMAPPED and UNINTEGRATED
ignore_pathways <- c("UNMAPPED", "UNINTEGRATED")
pathway_df <- pathway_df[!rownames(pathway_df) %in% ignore_pathways, ]

meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))
meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))

# Function to get ALL stats
get_stats <- function(meta_df, group_col, group1, group2, pathway_data) {
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  if (length(common_samples) == 0) return(NULL)
  
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group1, group2), ]
  common_samples <- meta_sub$sample_id
  path_sub <- path_sub[, common_samples]
  
  path_t <- as.data.frame(t(path_sub))
  keep_feat <- colSums(path_t > 0) >= 3
  path_t <- path_t[, keep_feat, drop=FALSE]
  
  results <- data.frame(Pathway = colnames(path_t), P_Val = NA, Log2FC = NA, LDA = NA, Enriched = NA)
  grp_vec <- meta_sub[[group_col]]
  
  for (i in 1:nrow(results)) {
    feat <- results$Pathway[i]
    val <- path_t[[feat]]
    
    # Wilcoxon
    p_val <- tryCatch({ wilcox.test(val ~ grp_vec)$p.value }, error = function(e) NA)
    
    # Log2FC
    means <- tapply(val, grp_vec, mean)
    m1 <- if (!is.na(means[group1])) means[group1] else 0
    m2 <- if (!is.na(means[group2])) means[group2] else 0
    log2fc <- log2((m1 + 1e-6) / (m2 + 1e-6))
    
    # LDA Score (Log10 diff)
    diff <- abs(m1 - m2)
    lda <- log10(1 + diff)
    
    enriched <- names(means)[which.max(means)]
    
    results$P_Val[i] <- p_val
    results$Log2FC[i] <- log2fc
    results$LDA[i] <- lda
    results$Enriched[i] <- enriched
    
    # Calculate Signed LDA
    # If Enriched == group1 -> LDA
    # If Enriched == group2 -> -LDA
    if (!is.na(enriched) && enriched == group1) {
      results$SignedLDA[i] <- lda
    } else if (!is.na(enriched) && enriched == group2) {
      results$SignedLDA[i] <- -lda
    } else {
      results$SignedLDA[i] <- 0
    }
  }
  
  results$P_Adj <- p.adjust(results$P_Val, method = "BH")
  return(results)
}

# Run Analysis
res_p1_p4 <- get_stats(meta_end, "Group", "P1", "P4", pathway_df)
res_rec_non <- get_stats(meta_dyn, "Group", "Recover", "Non_recover", pathway_df)

# Plot Function
plot_volcano <- function(res, title, group1_name, group2_name, target_label_group) {
  
  res$Significance <- "NS"
  # User request: LDA > 1.5, p < 0.05
  res$Significance[res$P_Val < 0.05 & abs(res$SignedLDA) > 1.5] <- "Significant"
  
  # Direction color
  res$Color <- "NS"
  # SignedLDA > 0 -> group1 enriched (Red)
  # SignedLDA < 0 -> group2 enriched (Blue)
  
  # Ensure we color based on enrichment, regardless of significance for now to visualize distribution,
  # but usually only significant ones get colored.
  # Let's color significant ones by their direction.
  res$Color[res$Significance == "Significant" & res$SignedLDA > 0] <- group1_name
  res$Color[res$Significance == "Significant" & res$SignedLDA < 0] <- group2_name
  
  # Top 10 Significant Enriched in Target Group
  # If target_label_group is group1, look for SignedLDA > 0
  # If target_label_group is group2, look for SignedLDA < 0
  
  if (target_label_group == group1_name) {
    target_hits <- res %>% filter(Significance == "Significant" & SignedLDA > 0)
  } else {
    target_hits <- res %>% filter(Significance == "Significant" & SignedLDA < 0)
  }
  
  top_hits <- target_hits %>%
    arrange(P_Val) %>%
    head(10)
  
  p <- ggplot(res, aes(x = SignedLDA, y = -log10(P_Val), color = Color)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = setNames(c("red", "blue", "grey"), c(group1_name, group2_name, "NS"))) +
    theme_bw() +
    labs(title = title, x = paste0("LDA Score (", group1_name, " / ", group2_name, ")"), y = "-Log10 P-value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey") +
    geom_text_repel(data = top_hits, aes(label = Pathway), size = 3, max.overlaps = 20, box.padding = 0.5)
  
  return(p)
}

# Generate Plots
if (!is.null(res_p1_p4)) {
  # Comparison: P1 vs P4.
  # Original SignedLDA: P1 (>0), P4 (<0).
  # User wants P4 / P1.
  # So we flip the sign: P4 (>0), P1 (<0).
  res_p1_p4_flipped <- res_p1_p4
  res_p1_p4_flipped$SignedLDA <- -res_p1_p4_flipped$SignedLDA
  
  # Now Positive = P4, Negative = P1.
  # Label top 10 enriched in P4.
  
  p1 <- plot_volcano(res_p1_p4_flipped, "Volcano: P4 vs P1 (LDA>1.5, P<0.05)", "P4", "P1", "P4")
  ggsave(file.path(output_dir, "volcano_pathway_P4_vs_P1.pdf"), p1, width = 10, height = 8)
  print("Saved volcano_pathway_P4_vs_P1.pdf")
}

if (!is.null(res_rec_non)) {
  # Comparison: Recover vs Non_recover.
  # Original SignedLDA: Recover (>0), Non (<0).
  # User wants to label "Recover".
  
  p2 <- plot_volcano(res_rec_non, "Volcano: Recover vs Non_recover (LDA>1.5, P<0.05)", "Recover", "Non_recover", "Recover")
  ggsave(file.path(output_dir, "volcano_pathway_Recover_vs_Non.pdf"), p2, width = 10, height = 8)
  print("Saved volcano_pathway_Recover_vs_Non.pdf")
}
