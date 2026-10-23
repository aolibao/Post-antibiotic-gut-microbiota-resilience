
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

# Load Metadata
meta_end <- read.csv(meta_endpoint_path, stringsAsFactors = FALSE)
colnames(meta_end) <- trimws(colnames(meta_end))
meta_dyn <- read.csv(meta_dynamic_path, stringsAsFactors = FALSE)
colnames(meta_dyn) <- trimws(colnames(meta_dyn))

# 3. Wilcoxon Analysis Function
run_wilcox <- function(meta_df, group_col, group1, group2, pathway_data, label) {
  print(paste("Running Wilcoxon Analysis:", label, "-", group1, "vs", group2))
  
  # Align samples
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  print(paste("Aligned Samples:", length(common_samples)))
  
  if (length(common_samples) == 0) return(NULL)
  
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  # Filter groups
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group1, group2), ]
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
    Mean_G1 = NA,
    Mean_G2 = NA,
    Log2FC = NA,
    P_Value = NA,
    stringsAsFactors = FALSE
  )
  
  grp_vec <- meta_sub[[group_col]]
  
  for (i in 1:nrow(results)) {
    feat <- results$Pathway[i]
    val <- path_t[[feat]]
    
    # Means
    means <- tapply(val, grp_vec, mean)
    m1 <- if (!is.na(means[group1])) means[group1] else 0
    m2 <- if (!is.na(means[group2])) means[group2] else 0
    
    # Log2FC (G1 / G2)
    log2fc <- log2((m1 + 1e-6) / (m2 + 1e-6))
    
    # Wilcoxon Test
    p_val <- tryCatch({
      wilcox.test(val ~ grp_vec)$p.value
    }, error = function(e) NA)
    
    results$Mean_G1[i] <- m1
    results$Mean_G2[i] <- m2
    results$Log2FC[i] <- log2fc
    results$P_Value[i] <- p_val
  }
  
  # FDR Correction
  results$P_Adjusted <- p.adjust(results$P_Value, method = "BH")
  results$Comparison <- label
  
  # Rename Mean columns to actual group names
  colnames(results)[colnames(results) == "Mean_G1"] <- paste0("Mean_", group1)
  colnames(results)[colnames(results) == "Mean_G2"] <- paste0("Mean_", group2)
  
  return(results)
}

# 4. Run Comparisons
# Comparison 1: P1 vs P4 (Log2FC = P1 / P4)
res_p1_p4 <- run_wilcox(meta_end, "Group", "P1", "P4", pathway_df, "P1_vs_P4")

# Comparison 2: Recover vs Non_recover (Log2FC = Recover / Non_recover)
res_rec_non <- run_wilcox(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, "Recover_vs_Non")

# 5. Save Results
if (!is.null(res_p1_p4)) {
  write.csv(res_p1_p4, file.path(output_dir, "wilcox_pathway_P1_vs_P4_log2fc.csv"), row.names = FALSE)
  print(paste("P1 vs P4 Results Saved. Rows:", nrow(res_p1_p4)))
}

if (!is.null(res_rec_non)) {
  write.csv(res_rec_non, file.path(output_dir, "wilcox_pathway_Recover_vs_Non_log2fc.csv"), row.names = FALSE)
  print(paste("Recover vs Non_recover Results Saved. Rows:", nrow(res_rec_non)))
}
