library(tidyverse)
library(ggplot2)

cat("============================================================\n")
cat("AIM 4.3 - GDPD2 vs Immune Correlation Analysis (CLEAN)\n")
cat("============================================================\n\n")

# ==============================================================================
# 1. LOAD EXPRESSION AND IMMUNE DATA
# ==============================================================================

expr <- read.csv(
  "TCGA_LUAD_WhiteBlack_VST.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_matrix <- as.matrix(expr)

if (!"GDPD2" %in% rownames(expr_matrix)) {
  stop("CRITICAL ERROR: GDPD2 not found in expression matrix.")
}

# Extract GDPD2 vector and establish clean tracking dataframe
gdpd2_df <- data.frame(
  SampleID = colnames(expr_matrix),
  GDPD2 = as.numeric(expr_matrix["GDPD2", ]),
  stringsAsFactors = FALSE
)

# Load your previously generated quanTIseq matrix
immune <- read.csv("Aim4_1_quanTIseq_immune_matrix.csv", check.names = FALSE)

# ==============================================================================
# 2. MERGE DATASETS AND PURGE NA METRICS
# ==============================================================================

merged <- merge(immune, gdpd2_df, by = "SampleID")
cat("Merged patient dimensions:", dim(merged)[1], "patients x", dim(merged)[2], "features\n")

# Dynamically identify cell column labels
immune_cols <- setdiff(colnames(merged), c("SampleID", "GDPD2", "Other"))

# FIXED: Removed the disruptive character-coercion block.
# Ensure column structures are directly registered as purely numeric.
merged[immune_cols] <- lapply(merged[immune_cols], as.numeric)

# Drop any entirely empty columns if they exist
immune_cols <- immune_cols[sapply(merged[immune_cols], function(x) !all(is.na(x)))]
cat("Valid immune cell features evaluated:", length(immune_cols), "\n")

# ==============================================================================
# 3. SPEARMAN CORRELATION ANALYSIS LOOP
# ==============================================================================

cor_results <- data.frame()

for (cell in immune_cols) {
  x_val <- merged[[cell]]
  y_val <- merged$GDPD2

  cor_test <- cor.test(
    x_val,
    y_val,
    method = "spearman",
    use = "complete.obs",
    exact = FALSE
  )

  cor_results <- rbind(
    cor_results,
    data.frame(
      CellType    = cell,
      Correlation = as.numeric(cor_test$estimate),
      p_value     = as.numeric(cor_test$p.value),
      stringsAsFactors = FALSE
    )
  )
}

# Adjust p-values using False Discovery Rate (Benjamini-Hochberg)
cor_results$adj_p <- p.adjust(cor_results$p_value, method = "BH")
write.csv(cor_results, "Aim4_3_Correlation_Results.csv", row.names = FALSE)

# ==============================================================================
# 4. REPLACED HEATMAP: GGPLOT2 SINGLE-COLUMN CORRELATION TILE PLOT
# ==============================================================================

# Replaced the failing base heatmap() function with an explicit ggplot matrix tile.
# This avoids the 'x must have at least 2 rows and columns' dimension crash.
heatmap_plot <- ggplot(cor_results, aes(x = "GDPD2", y = reorder(CellType, Correlation))) +
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, adj_p)), size = 3) +
  scale_fill_gradient2(low = "#2E9FDF", mid = "white", high = "#E7B800", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "GDPD2 vs Immune Infiltration Correlation",
    subtitle = "Spearman Rho & BH-Adjusted P-values",
    x = "Target Feature",
    y = "Immune Cell Type Fractions",
    fill = "Spearman Rho"
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    axis.text.y  = element_text(face = "bold", size = 10),
    axis.text.x  = element_text(face = "bold", size = 12)
  )

ggsave("Aim4_3_Correlation_Heatmap.png", heatmap_plot, width = 5, height = 7, dpi = 300)
cat("-> Heatmap generated successfully via ggplot2 engine.\n")

# ==============================================================================
# 5. SCATTER PLOTS FOR TOP TARGETED ASSOCIATIONS
# ==============================================================================

# Arrange and grab cells showing the strongest statistical evidence
top_cells <- cor_results %>%
  arrange(p_value) %>%
  head(6) %>%
  pull(CellType)

cat("\nStrongest immune associations isolated:\n")
print(top_cells)

for (cell in top_cells) {
  p_scatter <- ggplot(merged, aes(x = GDPD2, y = .data[[cell]])) +
    geom_point(alpha = 0.4, color = "#2E9FDF") +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", alpha = 0.2) +
    theme_bw() +
    labs(
      title = paste("GDPD2 Correlation with", cell),
      x = "GDPD2 Expression (VST)",
      y = "quanTIseq Immune Fraction Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.title = element_text(size = 9)
    )

  ggsave(
    filename = paste0("Aim4_3_", gsub("\\.", "_", cell), "_scatter.png"),
    plot     = p_scatter,
    width    = 5,
    height   = 4,
    dpi      = 300
  )
}

# ==============================================================================
# 6. REGISTER MERGED DATASETS
# ==============================================================================

write.csv(merged, "Aim4_3_Merged_GDPD2_Immune.csv", row.names = FALSE)
cat("\nAIM 4.3 COMPLETED SUCCESSFULLY\n")
