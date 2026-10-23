
library(limma)
library(dplyr)

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

# 3. Limma Analysis Function
run_limma <- function(meta_df, group_col, group_test, group_ref, pathway_data, label) {
  print(paste("Running Limma Analysis:", label, "-", group_test, "vs", group_ref))
  
  # Align samples
  common_samples <- intersect(meta_df$sample_id, colnames(pathway_data))
  print(paste("Aligned Samples:", length(common_samples)))
  
  if (length(common_samples) == 0) return(NULL)
  
  meta_sub <- meta_df[meta_df$sample_id %in% common_samples, ]
  path_sub <- pathway_data[, common_samples]
  
  # Filter groups
  meta_sub <- meta_sub[meta_sub[[group_col]] %in% c(group_test, group_ref), ]
  common_samples <- meta_sub$sample_id
  path_sub <- path_sub[, common_samples]
  
  if (nrow(meta_sub) == 0) return(NULL)
  
  # Filter low abundance
  keep_feat <- rowSums(path_sub > 0) >= 3
  path_sub <- path_sub[keep_feat, ]
  
  # Log Transformation (CPM is similar to Count, often Log2(CPM+1) is used for Limma)
  # Limma voom is for counts, but here we have CPM. Standard Limma on log-CPM is appropriate.
  log_cpm <- log2(path_sub + 1)
  
  # Design Matrix
  group_factor <- factor(meta_sub[[group_col]], levels = c(group_ref, group_test))
  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", paste0(group_test, "_vs_", group_ref))
  
  # Fit Model
  fit <- lmFit(log_cpm, design)
  fit <- eBayes(fit)
  
  # Get Results
  # Coefficient 2 corresponds to the comparison
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$Pathway <- rownames(res)
  res$Comparison <- label
  
  return(res)
}

# 4. Run Comparisons
# Comparison 1: P4 vs P1 (P4 is Test, P1 is Ref)
res_p4_p1 <- run_limma(meta_end, "Group", "P4", "P1", pathway_df, "P4_vs_P1")

# Comparison 2: Recover vs Non_recover (Recover is Test, Non_recover is Ref)
# Usually we want to see what is enriched in Recover.
# Or if we want to match P4 (Non-resilient?) vs P1 (Resilient?),
# Non_recover vs Recover might be more analogous?
# User input: "Recover / Non". Usually means Recover vs Non_recover.
res_rec_non <- run_limma(meta_dyn, "Group", "Recover", "Non_recover", pathway_df, "Recover_vs_Non")

# 5. Save Results
if (!is.null(res_p4_p1)) {
  write.csv(res_p4_p1, file.path(output_dir, "limma_pathway_P4_vs_P1.csv"), row.names = FALSE)
  print(paste("P4 vs P1 Results Saved. Top P-value:", min(res_p4_p1$P.Value)))
}

if (!is.null(res_rec_non)) {
  write.csv(res_rec_non, file.path(output_dir, "limma_pathway_Recover_vs_Non.csv"), row.names = FALSE)
  print(paste("Recover vs Non_recover Results Saved. Top P-value:", min(res_rec_non$P.Value)))
}
