############################################################
# AIM2_SCRIPT_01_DESEQ2_RACE_ANALYSIS.R
#
# INPUT:
# TCGA_LUAD_WhiteBlack_Counts_Unique.csv
# TCGA_LUAD_WhiteBlack_Clinical.csv
#
# OUTPUT:
# Aim2_DESeq2_AllGenes.csv
# Aim2_DESeq2_FDR05.csv
# Aim2_DESeq2_Top100.csv
# Aim2_DESeq2_Report.csv
############################################################

library(DESeq2)

cat("=====================================\n")
cat("AIM 2 - DESEQ2 RACE ANALYSIS\n")
cat("=====================================\n\n")

############################################################
# LOAD DATA
############################################################

counts <- read.csv(
  "TCGA_LUAD_WhiteBlack_Counts_Unique.csv",
  check.names = FALSE
)

clinical <- read.csv(
  "TCGA_LUAD_WhiteBlack_Clinical.csv",
  check.names = FALSE
)

############################################################
# BUILD COUNT MATRIX
############################################################

gene_names <- counts[,1]

count_matrix <- as.matrix(
  counts[, -1]
)

rownames(count_matrix) <- gene_names

storage.mode(count_matrix) <- "integer"

############################################################
# MATCH CLINICAL DATA
############################################################

patient_ids <- substr(
  colnames(count_matrix),
  1,
  12
)

clinical <- clinical[
  match(
    patient_ids,
    clinical$PatientID
  ),
]

############################################################
# KEEP WHITE AND BLACK
############################################################

keep <- clinical$Race %in% c(
  "White",
  "Black"
)

clinical <- clinical[keep, ]

count_matrix <- count_matrix[
  ,
  keep
]

############################################################
# REMOVE LOW COUNT GENES
############################################################

keep_genes <- rowSums(
  count_matrix >= 10
) >= 10

count_matrix <- count_matrix[
  keep_genes,
]

############################################################
# DESEQ2 OBJECT
############################################################

clinical$Race <- factor(
  clinical$Race,
  levels = c(
    "White",
    "Black"
  )
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = clinical,
  design = ~ Race
)

############################################################
# RUN DESEQ2
############################################################

cat("Running DESeq2...\n\n")

dds <- DESeq(dds)

############################################################
# RESULTS
############################################################

res <- results(
  dds,
  contrast = c(
    "Race",
    "Black",
    "White"
  )
)

res_df <- as.data.frame(res)

res_df$Gene <- rownames(res_df)

res_df <- res_df[
  order(res_df$padj),
]

############################################################
# CLEAN RESULTS
############################################################

res_df <- res_df[
  !is.na(res_df$padj),
]

############################################################
# SIGNIFICANT GENES
############################################################

sig_genes <- res_df[
  res_df$padj < 0.05,
]

############################################################
# SAVE FILES
############################################################

write.csv(
  res_df,
  "Aim2_DESeq2_AllGenes.csv",
  row.names = FALSE
)

write.csv(
  sig_genes,
  "Aim2_DESeq2_FDR05.csv",
  row.names = FALSE
)

write.csv(
  head(sig_genes,100),
  "Aim2_DESeq2_Top100.csv",
  row.names = FALSE
)

############################################################
# REPORT
############################################################

report <- data.frame(
  Metric = c(
    "Patients",
    "White",
    "Black",
    "Genes_Tested",
    "FDR05_Genes"
  ),
  Value = c(
    nrow(clinical),
    sum(clinical$Race=="White"),
    sum(clinical$Race=="Black"),
    nrow(res_df),
    nrow(sig_genes)
  )
)

write.csv(
  report,
  "Aim2_DESeq2_Report.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat("=====================================\n")
cat("DESEQ2 COMPLETE\n")
cat("=====================================\n\n")

cat("Patients:\n")
print(nrow(clinical))

cat("\nRace Distribution:\n")
print(table(clinical$Race))

cat("\nGenes Tested:\n")
print(nrow(res_df))

cat("\nFDR < 0.05 Genes:\n")
print(nrow(sig_genes))

cat("\nFiles Created:\n")
cat("Aim2_DESeq2_AllGenes.csv\n")
cat("Aim2_DESeq2_FDR05.csv\n")
cat("Aim2_DESeq2_Top100.csv\n")
cat("Aim2_DESeq2_Report.csv\n")

cat("\n=====================================\n")