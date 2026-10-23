
library(MASS)
library(dplyr)
library(tidyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
meta_group_path <- file.path(base_dir, "species/kraken.P1.P4/endpoint_sample_metadata.csv")
output_dir <- file.path(base_dir, "result/Fig5")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data
# Metadata
meta <- read.csv(meta_group_path, stringsAsFactors = FALSE)
colnames(meta) <- trimws(colnames(meta))
# Ensure columns exist
if (!all(c("sample_id", "Group", "timepoint") %in% colnames(meta))) {
  stop("Metadata missing required columns: sample_id, Group, timepoint")
}

# Filter P1/P4
meta <- subset(meta, Group %in% c("P1", "P4"))
rownames(meta) <- meta$sample_id

print("Metadata loaded:")
print(table(meta$timepoint, meta$Group))

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

# Filter out UNMAPPED/UNINTEGRATED
ignore_pathways <- c("UNMAPPED", "UNINTEGRATED")
pathway_df <- pathway_df[!rownames(pathway_df) %in% ignore_pathways, ]

# Align Samples
common_samples <- intersect(colnames(pathway_df), rownames(meta))
print(paste("Common samples:", length(common_samples)))

if (length(common_samples) == 0) {
  stop("No matching samples found between pathway data and metadata.")
}

pathway_clean <- pathway_df[, common_samples]
meta_clean <- meta[common_samples, ]
pathway_t <- as.data.frame(t(pathway_clean))

# 3. Analysis per Timepoint
results <- data.frame()
timepoints <- sort(unique(meta_clean$timepoint))

for (tp in timepoints) {
  print(paste("Analyzing Timepoint:", tp))
  
  tp_samples <- rownames(meta_clean)[meta_clean$timepoint == tp]
  tp_meta <- meta_clean[tp_samples, ]
  
  if (length(unique(tp_meta$Group)) < 2) {
    print(paste("Skipping", tp, "- Less than 2 groups."))
    next
  }
  
  tp_data <- pathway_t[tp_samples, ]
  
  # Filter low prevalence
  keep_feat <- colSums(tp_data > 0) >= 3
  tp_data <- tp_data[, keep_feat, drop=FALSE]
  
  p_values <- numeric(ncol(tp_data))
  names(p_values) <- colnames(tp_data)
  
  for (feat in colnames(tp_data)) {
    val <- tp_data[[feat]]
    grp <- tp_meta$Group
    
    p_val <- tryCatch({
      wilcox.test(val ~ grp)$p.value
    }, error = function(e) NA)
    
    p_values[feat] <- p_val
  }
  
  # Filter Significant (P < 0.05, Unadjusted)
  sig_feats <- names(p_values)[!is.na(p_values) & p_values < 0.05]
  
  for (feat in sig_feats) {
    val <- tp_data[[feat]]
    grp <- tp_meta$Group
    
    means <- tapply(val, grp, mean)
    diff <- abs(means["P1"] - means["P4"])
    
    # LDA Score approximation: log10(1 + abs(mean_diff))
    # Note: Real LEfSe uses bootstrapping and LDA. This is a simplified approximation often used in R scripts.
    lda_score <- log10(1 + diff) 
    
    enriched <- names(means)[which.max(means)]
    
    if (lda_score > 1.0) {
      results <- rbind(results, data.frame(
        Timepoint = tp,
        Pathway = feat,
        P_Value = p_values[feat],
        LDA = lda_score,
        Enriched = enriched
      ))
    }
  }
}

# Save Results
if (nrow(results) > 0) {
  write.csv(results, file.path(output_dir, "lefse_pathway_P1_vs_P4_timepoints_lda1_pval0.05.csv"), row.names = FALSE)
  print(paste("Saved results with", nrow(results), "rows."))
} else {
  print("No significant results found.")
}
