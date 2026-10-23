
library(dplyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
output_dir <- file.path(base_dir, "result/Fig5")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Pathway Data
pathway_raw <- read.csv(pathway_path, check.names = FALSE)
colnames(pathway_raw)[1] <- "Pathway"
if (colnames(pathway_raw)[2] == "") {
  pathway_raw <- pathway_raw[, -2]
}

# Aggregate duplicates if any
if (any(duplicated(pathway_raw$Pathway))) {
  pathway_raw <- pathway_raw %>% group_by(Pathway) %>% summarise(across(where(is.numeric), sum))
}
pathway_df <- as.data.frame(pathway_raw)
rownames(pathway_df) <- pathway_df$Pathway
pathway_df$Pathway <- NULL

# Filter out UNMAPPED and UNINTEGRATED
ignore_pathways <- c("UNMAPPED", "UNINTEGRATED")
pathway_df <- pathway_df[!rownames(pathway_df) %in% ignore_pathways, ]

# Calculate Mean Abundance
mean_abundance <- rowMeans(pathway_df)
pathway_list <- data.frame(
  Pathway = rownames(pathway_df),
  Mean_CPM = mean_abundance
)

# Sort by abundance
pathway_list <- pathway_list[order(pathway_list$Mean_CPM, decreasing = TRUE), ]

# Save to CSV
write.csv(pathway_list, file.path(output_dir, "all_pathways_abundance.csv"), row.names = FALSE)

print(paste("Total pathways found:", nrow(pathway_list)))
print("Top 10 most abundant pathways:")
print(head(pathway_list, 10))
