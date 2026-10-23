
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Setup
base_dir <- "/home/lihui/qiaohan/16s测序--paa处理/ABX课题"
output_dir <- file.path(base_dir, "result/Fig5")
data_path <- file.path(output_dir, "intersect_pathways_stats.csv")

# 2. Load Data
df <- read.csv(data_path, stringsAsFactors = FALSE)

# 3. Define Highlight Pathways
highlight_pathways <- c(
  "PWY-7761: NAD salvage pathway II (PNC IV cycle)",
  "PYRIDNUCSYN-PWY: NAD de novo biosynthesis I (from aspartate)"
)

# Create Label Column
df$Label <- NA
df$Label[df$Pathway %in% highlight_pathways] <- df$Pathway[df$Pathway %in% highlight_pathways]

# Create Color Column
df$Color <- "Other"
df$Color[df$Pathway %in% highlight_pathways] <- "Highlight"

# 4. Plot
p <- ggplot(df, aes(x = LogFC_Meta3, y = LogFC_Meta2)) +
  geom_point(aes(color = Color, size = Color), alpha = 0.7) +
  scale_color_manual(values = c("Highlight" = "red", "Other" = "grey")) +
  scale_size_manual(values = c("Highlight" = 3, "Other" = 1)) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20, box.padding = 0.5, color = "black") +
  theme_bw() +
  labs(
    title = "Comparison of LogFC (Meta3 vs Meta2)",
    x = "LogFC (Recover / Non_recover)", # Assuming Meta3 is Rec/Non
    y = "LogFC (P1 / P4)" # Assuming Meta2 is P1/P4
  ) +
  theme(legend.position = "none")

# Save
ggsave(file.path(output_dir, "scatter_pathway_log2fc_highlight.pdf"), p, width = 8, height = 6)
print("Saved scatter_pathway_log2fc_highlight.pdf")
