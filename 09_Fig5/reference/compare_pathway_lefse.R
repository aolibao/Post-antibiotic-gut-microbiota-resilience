
library(MASS)
library(dplyr)
library(tidyr)

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

# Load Metadata
meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))

meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))

# 3. LEfSe Analysis Function
run_lefse <- function(meta_df, group_col, group1, group2, pathway_data, label) {
  print(paste("Running LEfSe Analysis:", label, "-", group1, "vs", group2))
  
  # Align samples
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
  
  # Filter low abundance features
  keep_feat <- colSums(path_t > 0) >= 3
  path_t <- path_t[, keep_feat, drop=FALSE]
  
  results <- data.frame()
  grp_vec <- meta_sub[[group_col]]
  
  p_values <- numeric(ncol(path_t))
  names(p_values) <- colnames(path_t)
  
  # Calculate Wilcoxon P-values
  for (feat in colnames(path_t)) {
    val <- path_t[[feat]]
    p_val <- tryCatch({
      wilcox.test(val ~ grp_vec)$p.value
    }, error = function(e) NA)
    p_values[feat] <- p_val
  }
  
  # FDR Correction
  p_adjusted <- p.adjust(p_values, method = "BH")
  
  # Identify significant features (Padj < 0.05)
  sig_feats <- names(p_adjusted)[!is.na(p_adjusted) & p_adjusted < 0.05]
  
  for (feat in sig_feats) {
    val <- path_t[[feat]]
    
    # LDA Score Calculation
    means <- tapply(val, grp_vec, mean)
    # Ensure both groups exist
    if (length(means) < 2) next
    
    diff <- abs(means[group1] - means[group2])
    lda_score <- log10(1 + diff) 
    
    enriched <- names(means)[which.max(means)]
    
    # Cutoff LDA > 1.5
    if (lda_score > 1.5) {
      results <- rbind(results, data.frame(
        Pathway = feat,
        P_Value = p_values[feat],
        P_Adjusted = p_adjusted[feat],
        LDA = lda_score,
        Enriched = enriched,
        Comparison = label
      ))
    }
  }
  
  return(results)
}

# 4. Run Comparisons
# 1. P1 vs P4 (Endpoint Meta)
res_p1_p4 <- run_lefse(meta_end, "Group", "P1", "P4", pathway_df, "P1_vs_P4")

# 2. Recover vs Non_recover (Dynamic Meta)
res_rec_non <- run_lefse(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, "Recover_vs_Non")

# 5. Save Results
if (!is.null(res_p1_p4)) {
  write.csv(res_p1_p4, file.path(output_dir, "lefse_pathway_P1_vs_P4_overall_padj0.05_lda1.5.csv"), row.names = FALSE)
  print(paste("P1 vs P4 Significant Pathways:", nrow(res_p1_p4)))
}

if (!is.null(res_rec_non)) {
  write.csv(res_rec_non, file.path(output_dir, "lefse_pathway_Recover_vs_Non_overall_padj0.05_lda1.5.csv"), row.names = FALSE)
  print(paste("Recover vs Non_recover Significant Pathways:", nrow(res_rec_non)))
}
