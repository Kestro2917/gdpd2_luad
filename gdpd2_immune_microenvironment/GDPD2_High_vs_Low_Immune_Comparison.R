############################################################
# AIM 4.2 - GDPD2 High vs Low Immune Comparison (ROBUST)
# TCGA-LUAD GDPD2 Project
############################################################

# =========================
# 1. Load libraries
# =========================

library(tidyverse)
library(ggplot2)

# =========================
# 2. Load expression matrix
# =========================

expr <- read.csv(
  "TCGA_LUAD_WhiteBlack_VST.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_matrix <- as.matrix(expr)

# =========================
# 3. Extract GDPD2 expression
# =========================

if (!"GDPD2" %in% rownames(expr_matrix)) {
  stop("GDPD2 not found in expression matrix")
}

gdpd2 <- as.numeric(expr_matrix["GDPD2", ])

gdpd2_df <- data.frame(
  SampleID = colnames(expr_matrix),
  GDPD2 = gdpd2
)

# =========================
# 4. Define High / Low groups
# =========================

median_val <- median(gdpd2_df$GDPD2, na.rm = TRUE)

gdpd2_df$Group <- ifelse(
  gdpd2_df$GDPD2 >= median_val,
  "High",
  "Low"
)

cat("Group distribution:\n")
print(table(gdpd2_df$Group))

# =========================
# 5. Load immune matrix (from Aim 4.1)
# =========================

immune <- read.csv("Aim4_1_quanTIseq_immune_matrix.csv")

# =========================
# 6. Merge datasets
# =========================

merged <- merge(immune, gdpd2_df, by = "SampleID")

cat("Merged dimensions:\n")
print(dim(merged))

# =========================
# 7. Identify immune columns
# =========================

immune_cols <- setdiff(
  colnames(merged),
  c("SampleID", "GDPD2", "Group")
)

# =========================
# 8. FORCE numeric conversion (CRITICAL FIX)
# =========================

merged[immune_cols] <- lapply(merged[immune_cols], function(x) {
  as.numeric(as.character(x))
})

# Remove invalid columns (all NA)
immune_cols <- immune_cols[
  sapply(merged[immune_cols], function(x) !all(is.na(x)))
]

cat("Valid immune features:", length(immune_cols), "\n")

# =========================
# 9. Wilcoxon test (SAFE)
# =========================

stat_results <- data.frame()

for (cell in immune_cols) {

  x <- merged[[cell]]

  if (!is.numeric(x)) next

  test <- wilcox.test(x ~ merged$Group)

  stat_results <- rbind(
    stat_results,
    data.frame(
      CellType = cell,
      p_value = test$p.value
    )
  )
}

# Adjust p-values
stat_results$adj_p <- p.adjust(stat_results$p_value, method = "BH")

write.csv(stat_results, "Aim4_2_StatResults.csv", row.names = FALSE)

# =========================
# 10. Identify significant features
# =========================

sig_cells <- stat_results %>%
  filter(adj_p < 0.05) %>%
  pull(CellType)

cat("Significant immune features:", length(sig_cells), "\n")

# =========================
# 11. Boxplots for significant features
# =========================

if (length(sig_cells) > 0) {

  for (cell in sig_cells) {

    p <- ggplot(merged, aes(x = Group, y = .data[[cell]], fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
      theme_classic() +
      labs(
        title = paste("GDPD2 High vs Low -", cell),
        y = "Immune Score"
      ) +
      theme(legend.position = "none")

    ggsave(
      paste0("Aim4_2_", cell, "_boxplot.png"),
      p,
      width = 5,
      height = 4,
      dpi = 300
    )
  }
}

# =========================
# 12. Full immune landscape plot
# =========================

long_df <- merged %>%
  pivot_longer(
    cols = all_of(immune_cols),
    names_to = "CellType",
    values_to = "Score"
  )

p_all <- ggplot(long_df, aes(x = CellType, y = Score, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GDPD2 High vs Low Immune Landscape")

ggsave("Aim4_2_AllImmune_Boxplot.png", p_all, width = 10, height = 6)

# =========================
# 13. Save merged dataset
# =========================

write.csv(merged, "Aim4_2_Merged_GDPD2_Immune.csv", row.names = FALSE)

# =========================
# DONE
# =========================

cat("\nAIM 4.2 COMPLETED SUCCESSFULLY\n")