
library(MASS)
library(dplyr)
library(tidyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
meta_full_path <- file.path(base_dir, "species/kranken.all/meta.csv")
meta_group_path <- file.path(base_dir, "species/kraken.P1.P4/endpoint_sample_metadata.csv")
output_dir <- file.path(base_dir, "result/Fig5")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data
# Metadata
meta_full <- read.csv(meta_full_path, stringsAsFactors = FALSE)
colnames(meta_full) <- trimws(colnames(meta_full))
meta_full$sample_id <- trimws(meta_full$ids)

meta_group <- read.csv(meta_group_path, stringsAsFactors = FALSE)
colnames(meta_group) <- trimws(colnames(meta_group))
subject_group_map <- unique(meta_group[, c("subject", "Group")])
subject_group_map <- subject_group_map[!duplicated(subject_group_map$subject), ]

meta_merged <- merge(meta_full, subject_group_map, by = "subject", all.x = TRUE)
# Filter for P1/P4 groups
meta_clean <- subset(meta_merged, Group %in% c("P1", "P4"))
rownames(meta_clean) <- meta_clean$sample_id

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
print(paste("Pathways after filtering UNMAPPED/UNINTEGRATED:", nrow(pathway_df)))

# Align Samples
pathway_cols <- colnames(pathway_df)
common_ids <- intersect(pathway_cols, meta_clean$sample_id)
print(paste("Aligned Samples:", length(common_ids)))

if (length(common_ids) == 0) {
  stop("No common samples found between pathway data and metadata.")
}

pathway_clean <- pathway_df[, common_ids]
meta_clean <- meta_clean[common_ids, ]

# Transpose for analysis (Samples as rows)
pathway_t <- as.data.frame(t(pathway_clean))

# 3. Overall P1 vs P4 Comparison (All Timepoints Aggregated)
print("Analyzing P1 vs P4 (All Timepoints Aggregated)")

results <- data.frame()

# Filter low abundance features
keep_feat <- colSums(pathway_t > 0) >= 3
pathway_data <- pathway_t[, keep_feat, drop=FALSE]

p_values <- numeric(ncol(pathway_data))
names(p_values) <- colnames(pathway_data)

# Calculate Wilcoxon P-values
for (feat in colnames(pathway_data)) {
  val <- pathway_data[[feat]]
  grp <- meta_clean$Group
  
  p_val <- tryCatch({
    wilcox.test(val ~ grp)$p.value
  }, error = function(e) NA)
  
  p_values[feat] <- p_val
}

# FDR Correction (Not strictly required by user, but good to have)
p_adjusted <- p.adjust(p_values, method = "BH")

# Calculate LDA for ALL features
print("Calculating LDA for all pathways...")
all_feats <- colnames(pathway_data)

for (feat in all_feats) {
  val <- pathway_data[[feat]]
  grp <- meta_clean$Group
  
  # LDA Score Calculation
  means <- tapply(val, grp, mean)
  # Handle case where means might be NA if group is missing (unlikely here)
  if(any(is.na(means))) next
  
  diff <- abs(means["P1"] - means["P4"])
  lda_score <- log10(1 + diff) 
  
  enriched <- names(means)[which.max(means)]
  
  results <- rbind(results, data.frame(
    Pathway = feat,
    P_Value = p_values[feat],
    P_Adjusted = p_adjusted[feat],
    LDA = lda_score,
    Enriched = enriched
  ))
}

# Save Results (All)
if (nrow(results) > 0) {
  # Sort by P-value
  results <- results[order(results$P_Value), ]
  write.csv(results, file.path(output_dir, "lefse_pathway_P1_vs_P4_overall_all.csv"), row.names = FALSE)
  print(paste("Analysis complete. Processed", nrow(results), "pathways."))
} else {
  print("No results generated.")
}
