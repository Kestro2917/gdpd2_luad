############################################################
# AIM 4.4 - Immune Landscape Visualization (FINAL CLEAN)
# TCGA-LUAD GDPD2 Project
############################################################
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)

# =========================
# 1. Load input files
# =========================
immune <- read.csv("Aim4_1_quanTIseq_immune_matrix.csv")
stat_results <- read.csv("Aim4_2_StatResults.csv")
cor_results  <- read.csv("Aim4_3_Correlation_Results.csv")

# =========================
# 2. Identify immune columns
# =========================
immune_cols <- setdiff(colnames(immune), "SampleID")

# =========================
# 3. CLEAN immune matrix (CRITICAL FIX)
# =========================
immune_clean <- immune %>%
  dplyr::select(all_of(immune_cols)) %>%
  mutate(across(everything(), ~as.numeric(as.character(.))))

# Check how many NAs were introduced BEFORE zeroing them out
na_counts <- colSums(is.na(immune_clean))
if (any(na_counts > 0)) {
  cat("Warning: non-numeric values found and coerced to NA in these columns:\n")
  print(na_counts[na_counts > 0])
}

# replace NA safely
immune_clean[is.na(immune_clean)] <- 0

# reattach SampleID
immune_clean$SampleID <- immune$SampleID

cat("Clean immune matrix dimensions:\n")
print(dim(immune_clean))

# =========================
# 4. AIM 4.4A - Immune significance plot
# =========================
stat_results <- stat_results %>%
  mutate(neglogP = -log10(adj_p + 1e-10)) %>%
  arrange(adj_p)

p1 <- ggplot(stat_results, aes(
  x = reorder(CellType, neglogP),
  y = neglogP
)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "AIM 4.4A: Immune Differences (GDPD2 High vs Low)",
    x = "Immune Cell Type",
    y = "-log10(FDR)"
  )
ggsave("Aim4_4A_Immune_Significance.png", p1, width = 8, height = 5)

# =========================
# 5. AIM 4.4B - Correlation plot
# =========================
cor_results <- cor_results %>%
  mutate(Sign = ifelse(Correlation > 0, "Positive", "Negative")) %>%
  arrange(desc(abs(Correlation)))

p2 <- ggplot(cor_results, aes(
  x = reorder(CellType, Correlation),
  y = Correlation,
  fill = Sign
)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue")) +
  theme_classic() +
  labs(
    title = "AIM 4.4B: GDPD2–Immune Correlation",
    x = "Immune Cell Type",
    y = "Spearman Correlation"
  )
ggsave("Aim4_4B_Correlation.png", p2, width = 8, height = 5)

# =========================
# 6. AIM 4.4C - HEATMAP (SAFE NUMERIC MATRIX)
# =========================
mat_for_heatmap <- immune_clean[, immune_cols]

# Identify and drop zero-variance columns (these break scale() -> hclust())
col_sd <- apply(mat_for_heatmap, 2, sd, na.rm = TRUE)
zero_var_cols <- names(col_sd[col_sd == 0 | is.na(col_sd)])

if (length(zero_var_cols) > 0) {
  cat("Dropping zero-variance columns from heatmap (cannot be scaled):\n")
  print(zero_var_cols)
  mat_for_heatmap <- mat_for_heatmap[, !(colnames(mat_for_heatmap) %in% zero_var_cols)]
}

heatmap_matrix <- scale(mat_for_heatmap)
rownames(heatmap_matrix) <- immune_clean$SampleID

# Final safety net: replace any remaining NaN/Inf (shouldn't occur now, but just in case)
heatmap_matrix[!is.finite(heatmap_matrix)] <- 0

pheatmap(
  heatmap_matrix,
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "AIM 4.4C: Immune Landscape (quanTIseq)",
  filename = "Aim4_4C_Immune_Heatmap.png",
  width = 12,
  height = 8
)

# =========================
# 7. AIM 4.4D - Top immune boxplot (FIXED pivot_longer)
# =========================
top_cells <- cor_results %>%
  arrange(adj_p) %>%
  head(5) %>%
  pull(CellType)

long_df <- immune_clean %>%
  dplyr::select(SampleID, all_of(top_cells)) %>%
  pivot_longer(
    cols = all_of(top_cells),
    names_to = "CellType",
    values_to = "Score"
  ) %>%
  mutate(Score = as.numeric(Score))

p3 <- ggplot(long_df, aes(x = CellType, y = Score)) +
  geom_boxplot(fill = "lightgreen") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "AIM 4.4D: Top Immune Cell Distribution",
    x = "Immune Cell Type",
    y = "Score"
  )
ggsave("Aim4_4D_TopImmune_Boxplot.png", p3, width = 7, height = 5)

# =========================
# 8. FINAL COMBINED TABLE
# =========================
final_summary <- cor_results %>%
  inner_join(stat_results, by = "CellType")

write.csv(final_summary,
          "Aim4_4_Final_Integrated_Results.csv",
          row.names = FALSE)

cat("\n=================================\n")
cat(" AIM 4.4 COMPLETED SUCCESSFULLY\n")
cat("=================================\n")