############################################################

# WHITE-BLACK VST DIAGNOSTIC SCRIPT

############################################################

library(DESeq2)

############################################################

# LOAD FILES

############################################################

counts <- read.csv(
"TCGA_LUAD_WhiteBlack_Counts_Unique.csv",
check.names = FALSE,
stringsAsFactors = FALSE
)

clinical <- read.csv(
"TCGA_LUAD_WhiteBlack_Clinical.csv",
check.names = FALSE,
stringsAsFactors = FALSE
)

############################################################

# RAW FILE CHECK

############################################################

cat("\n=====================================\n")
cat("RAW COUNT FILE CHECK\n")
cat("=====================================\n")

cat("\nDimensions:\n")
print(dim(counts))

cat("\nFirst 5 column names:\n")
print(colnames(counts)[1:5])

cat("\nFirst 10 values in first column:\n")
print(head(counts[,1],10))

############################################################

# FIRST COLUMN CHECK

############################################################

gene_col <- colnames(counts)[1]

cat("\nGene column name:\n")
print(gene_col)

############################################################

# SET ROW NAMES

############################################################

rownames(counts) <- counts[[gene_col]]

counts <- counts[, -1]

cat("\n=====================================\n")
cat("AFTER ROWNAME ASSIGNMENT\n")
cat("=====================================\n")

cat("\nDimensions:\n")
print(dim(counts))

cat("\nNumber of rownames:\n")
print(length(rownames(counts)))

cat("\nFirst 10 rownames:\n")
print(head(rownames(counts),10))

############################################################

# CONVERT TO MATRIX

############################################################

count_matrix <- as.matrix(counts)

storage.mode(count_matrix) <- "integer"

cat("\n=====================================\n")
cat("COUNT MATRIX CHECK\n")
cat("=====================================\n")

cat("\nDimensions:\n")
print(dim(count_matrix))

cat("\nNumber of rownames:\n")
print(length(rownames(count_matrix)))

cat("\nFirst 10 rownames:\n")
print(head(rownames(count_matrix),10))

############################################################

# SAMPLE CHECK

############################################################

sample_ids <- colnames(count_matrix)

patient_ids <- substr(
sample_ids,
1,
12
)

cat("\n=====================================\n")
cat("SAMPLE CHECK\n")
cat("=====================================\n")

cat("\nNumber of samples:\n")
print(length(sample_ids))

cat("\nFirst 10 sample IDs:\n")
print(head(sample_ids,10))

############################################################

# MATCH CLINICAL

############################################################

clinical <- clinical[
match(
patient_ids,
clinical$PatientID
),
]

cat("\n=====================================\n")
cat("CLINICAL MATCHING\n")
cat("=====================================\n")

cat("\nClinical dimensions:\n")
print(dim(clinical))

cat("\nMatched patients:\n")
print(sum(!is.na(clinical$PatientID)))

############################################################

# DESEQ OBJECT

############################################################

dds <- DESeqDataSetFromMatrix(
countData = count_matrix,
colData = clinical,
design = ~ 1
)

cat("\n=====================================\n")
cat("DESEQ OBJECT CHECK\n")
cat("=====================================\n")

cat("\nNumber of genes in dds:\n")
print(nrow(dds))

cat("\nLength of rownames(dds):\n")
print(length(rownames(dds)))

cat("\nFirst 10 genes:\n")
print(head(rownames(dds),10))

############################################################

# VST NORMALIZATION

############################################################

cat("\nRunning VST...\n")

vst_obj <- vst(
dds,
blind = TRUE
)

vst_matrix <- assay(vst_obj)

############################################################

# VST CHECK

############################################################

cat("\n=====================================\n")
cat("VST MATRIX CHECK\n")
cat("=====================================\n")

cat("\nDimensions:\n")
print(dim(vst_matrix))

cat("\nLength of rownames(vst_matrix):\n")
print(length(rownames(vst_matrix)))

cat("\nFirst 10 rownames:\n")
print(head(rownames(vst_matrix),10))

############################################################

# FINAL DIAGNOSIS

############################################################

cat("\n=====================================\n")
cat("FINAL DIAGNOSIS\n")
cat("=====================================\n")

cat("\nCount matrix genes:\n")
print(nrow(count_matrix))

cat("\nVST matrix genes:\n")
print(nrow(vst_matrix))

cat("\nLength of VST rownames:\n")
print(length(rownames(vst_matrix)))

cat("\nFirst 20 VST rownames:\n")
print(head(rownames(vst_matrix),20))

############################################################

# SAVE REPORT

############################################################

sink("WhiteBlack_VST_Diagnostic_Report.txt")

cat("White-Black VST Diagnostic Completed\n")

sink()

cat("\nDiagnostic report saved:\n")
cat("WhiteBlack_VST_Diagnostic_Report.txt\n")
############################################################

# END

############################################################
