############################################################
# AIM2 STEP 2
# OVERLAP ANALYSIS
############################################################

de_genes <- read.csv(
  "Aim2_DESeq2_FDR05.csv",
  check.names = FALSE
)

lasso_genes <- read.csv(
  "WhiteBlack_LASSO_SelectedGenes.csv",
  check.names = FALSE
)

############################################################
# OVERLAP
############################################################

overlap <- merge(
  de_genes,
  lasso_genes,
  by = "Gene"
)

############################################################
# SORT
############################################################

overlap <- overlap[
  order(overlap$padj),
]

############################################################
# SAVE
############################################################

write.csv(
  overlap,
  "Aim2_Race_Prognostic_Genes.csv",
  row.names = FALSE
)

############################################################
# REPORT
############################################################

report <- data.frame(
  Metric = c(
    "Race_Associated_Genes",
    "LASSO_Genes",
    "Overlap_Genes"
  ),
  Value = c(
    nrow(de_genes),
    nrow(lasso_genes),
    nrow(overlap)
  )
)

write.csv(
  report,
  "Aim2_Race_Prognostic_Report.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat("\n=====================================\n")
cat("AIM2 OVERLAP COMPLETE\n")
cat("=====================================\n")

cat("\nRace-associated genes:\n")
print(nrow(de_genes))

cat("\nLASSO genes:\n")
print(nrow(lasso_genes))

cat("\nOverlap genes:\n")
print(nrow(overlap))

cat("\nTop overlap genes:\n")
print(
  head(
    overlap[, c(
      "Gene",
      "log2FoldChange",
      "padj",
      "Coefficient"
    )],
    20
  )
)

cat("\n=====================================\n")