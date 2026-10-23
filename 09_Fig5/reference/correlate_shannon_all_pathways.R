
library(vegan)
library(dplyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
species_path <- file.path(base_dir, "species/kranken.all/combine_brk_all.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data

# Species Data -> Shannon Index
species_raw <- read.csv(species_path, check.names = FALSE, stringsAsFactors = FALSE)
# Assuming first column is ID or taxonomy
rownames(species_raw) <- species_raw[, 1]
species_raw <- species_raw[, -1]
# Clean column names
colnames(species_raw) <- gsub("\\.brk\\.report$", "", colnames(species_raw))

species_t <- t(species_raw)
shannon_index <- diversity(species_t, index = "shannon")

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

# 3. Align Samples
common_samples <- intersect(colnames(pathway_df), names(shannon_index))
print(paste("Common Samples:", length(common_samples)))

if (length(common_samples) < 3) {
  stop("Not enough common samples for correlation analysis.")
}

pathway_clean <- pathway_df[, common_samples]
shannon_clean <- shannon_index[common_samples]

# 4. Correlation Analysis Loop
results <- data.frame(
  Pathway = rownames(pathway_clean),
  Rho = numeric(nrow(pathway_clean)),
  P_Value = numeric(nrow(pathway_clean)),
  stringsAsFactors = FALSE
)

# Transpose pathway for easier access (rows = samples) if using apply, but here we loop over rows
print("Calculating correlations...")
counter <- 0

for (i in 1:nrow(pathway_clean)) {
  path_vals <- as.numeric(pathway_clean[i, ])
  
  # Check if variance is zero
  if (sd(path_vals) == 0) {
    results$Rho[i] <- NA
    results$P_Value[i] <- NA
  } else {
    test <- cor.test(path_vals, shannon_clean, method = "spearman", exact = FALSE)
    results$Rho[i] <- test$estimate
    results$P_Value[i] <- test$p.value
  }
  
  counter <- counter + 1
  if (counter %% 100 == 0) cat(".")
}
cat("\n")

# Remove NAs
results <- na.omit(results)

# FDR Correction
results$FDR <- p.adjust(results$P_Value, method = "BH")

# Sort by P-value
results <- results[order(results$P_Value), ]

# 5. Save Results
output_file <- file.path(output_dir, "correlation_shannon_pathways_stats.csv")
write.csv(results, output_file, row.names = FALSE)

print(paste("Correlation analysis complete. Saved to:", output_file))
print("Top 10 correlated pathways:")
print(head(results, 10))
