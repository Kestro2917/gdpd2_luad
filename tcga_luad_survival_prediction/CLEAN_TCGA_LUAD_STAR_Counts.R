############################################################
# CLEAN TCGA_LUAD_STAR_Counts.csv
############################################################

library(dplyr)

############################################################
# LOAD FILE
############################################################

expr_raw <- read.csv(
  "TCGA_LUAD_STAR_Counts.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Original dimensions:\n")
print(dim(expr_raw))

############################################################
# IDENTIFY GENE COLUMN
############################################################

gene_col <- expr_raw[,1]

############################################################
# REMOVE EMPTY GENE NAMES
############################################################

expr_raw <- expr_raw[
  !is.na(gene_col) &
  gene_col != "",
]

############################################################
# REMOVE NON-GENE ROWS
############################################################

expr_raw <- expr_raw[
  !grepl(
    "^__",
    expr_raw[,1]
  ),
]

############################################################
# REMOVE POSSIBLE METADATA ROWS
############################################################

expr_raw <- expr_raw[
  !tolower(expr_raw[,1]) %in%
    c(
      "gene",
      "gene_name",
      "symbol",
      "description"
    ),
]

############################################################
# RENAME FIRST COLUMN
############################################################

colnames(expr_raw)[1] <- "Gene"

############################################################
# CONVERT EXPRESSION TO NUMERIC
############################################################

expr_raw[,-1] <- lapply(
  expr_raw[,-1],
  function(x) as.numeric(as.character(x))
)

############################################################
# CALCULATE MEAN EXPRESSION
############################################################

expr_raw$MeanExpression <- rowMeans(
  expr_raw[,-1],
  na.rm = TRUE
)

############################################################
# COLLAPSE DUPLICATE GENES
############################################################

expr_clean <- expr_raw %>%
  arrange(desc(MeanExpression)) %>%
  distinct(Gene, .keep_all = TRUE)

############################################################
# REMOVE HELPER COLUMN
############################################################

expr_clean$MeanExpression <- NULL

############################################################
# REPORT
############################################################

cat("\nGenes before cleaning:\n")
print(nrow(expr_raw))

cat("\nGenes after duplicate removal:\n")
print(nrow(expr_clean))

cat("\nDuplicates removed:\n")
print(
  nrow(expr_raw) -
  nrow(expr_clean)
)

############################################################
# CHECK DUPLICATES
############################################################

cat("\nRemaining duplicate genes:\n")
print(
  sum(
    duplicated(expr_clean$Gene)
  )
)

############################################################
# SET ROW NAMES
############################################################

rownames(expr_clean) <- expr_clean$Gene

expr_clean$Gene <- NULL

############################################################
# SAVE CLEAN FILE
############################################################

write.csv(
  expr_clean,
  "TCGA_LUAD_STAR_Counts_Clean.csv",
  quote = FALSE
)

cat("\nCleaned files saved:\n")
cat("TCGA_LUAD_STAR_Counts_Clean.csv\n")
