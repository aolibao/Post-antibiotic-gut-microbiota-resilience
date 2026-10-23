
library(vegan)
library(dplyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
pathway_path <- file.path(base_dir, "pathway/pathway_cpm_unstratified.csv")
species_path <- file.path(base_dir, "species/kranken.all/combine_brk_all.csv")
meta_path <- file.path(base_dir, "species/kraken.P1.P4/endpoint_sample_metadata.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Data

# Metadata
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
colnames(meta) <- trimws(colnames(meta))
# Ensure we have sample_id and timepoint
if (!all(c("sample_id", "timepoint") %in% colnames(meta))) {
  stop("Metadata missing required columns: sample_id, timepoint")
}
rownames(meta) <- meta$sample_id

# Species Data -> Shannon Index
species_raw <- read.csv(species_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(species_raw) <- species_raw[, 1]
species_raw <- species_raw[, -1]
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

# 3. Analysis per Timepoint
results_all <- data.frame()
timepoints <- sort(unique(meta$timepoint))

print(paste("Timepoints found:", paste(timepoints, collapse=", ")))

for (tp in timepoints) {
  print(paste("Analyzing Timepoint:", tp))
  
  # Get samples for this timepoint
  tp_samples_meta <- rownames(meta)[meta$timepoint == tp]
  
  # Intersect with pathway and species data
  common_samples <- intersect(tp_samples_meta, intersect(colnames(pathway_df), names(shannon_index)))
  
  print(paste("  Samples:", length(common_samples)))
  
  if (length(common_samples) < 3) {
    print("  Skipping: Not enough samples.")
    next
  }
  
  # Subset data
  pathway_sub <- pathway_df[, common_samples]
  shannon_sub <- shannon_index[common_samples]
  
  # Correlation Loop
  tp_results <- data.frame(
    Timepoint = tp,
    Pathway = rownames(pathway_sub),
    Rho = numeric(nrow(pathway_sub)),
    P_Value = numeric(nrow(pathway_sub)),
    stringsAsFactors = FALSE
  )
  
  counter <- 0
  for (i in 1:nrow(pathway_sub)) {
    path_vals <- as.numeric(pathway_sub[i, ])
    
    if (sd(path_vals) == 0) {
      tp_results$Rho[i] <- NA
      tp_results$P_Value[i] <- NA
    } else {
      test <- cor.test(path_vals, shannon_sub, method = "spearman", exact = FALSE)
      tp_results$Rho[i] <- test$estimate
      tp_results$P_Value[i] <- test$p.value
    }
  }
  
  # Remove NAs
  tp_results <- na.omit(tp_results)
  
  # FDR Correction (within timepoint)
  tp_results$FDR <- p.adjust(tp_results$P_Value, method = "BH")
  
  results_all <- rbind(results_all, tp_results)
}

# 4. Save Combined Results
if (nrow(results_all) > 0) {
  output_file <- file.path(output_dir, "correlation_shannon_pathways_by_timepoint.csv")
  write.csv(results_all, output_file, row.names = FALSE)
  print(paste("Analysis complete. Saved to:", output_file))
  print("Sample of results:")
  print(head(results_all))
} else {
  print("No valid results generated.")
}
