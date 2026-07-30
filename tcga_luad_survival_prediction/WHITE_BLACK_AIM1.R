############################################################

# WHITE-BLACK AIM1

# STEP 1: VST NORMALIZATION

############################################################

library(DESeq2)

############################################################

# LOAD DATA

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

# FIX GENE NAMES

############################################################

gene_names <- counts[,1]

rownames(counts) <- gene_names

counts <- counts[, -1]

############################################################

# CONVERT TO INTEGER MATRIX

############################################################

count_matrix <- as.matrix(counts)

rownames(count_matrix) <- gene_names

storage.mode(count_matrix) <- "integer"

############################################################

# MATCH CLINICAL TO EXPRESSION

############################################################

sample_ids <- colnames(count_matrix)

patient_ids <- substr(
sample_ids,
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

# QC CHECKS

############################################################

cat("\n=====================================\n")
cat("WHITE-BLACK VST QC\n")
cat("=====================================\n")

cat("\nExpression dimensions:\n")
print(dim(count_matrix))

cat("\nClinical dimensions:\n")
print(dim(clinical))

cat("\nMatched patients:\n")
print(sum(!is.na(clinical$PatientID)))

cat("\nGene names retained:\n")
print(length(rownames(count_matrix)))

cat("\nFirst 10 genes:\n")
print(head(rownames(count_matrix),10))

############################################################

# DESEQ OBJECT

############################################################

dds <- DESeqDataSetFromMatrix(
countData = count_matrix,
colData = clinical,
design = ~ 1
)

############################################################

# VST NORMALIZATION

############################################################

cat("\nRunning VST normalization...\n")

vst_obj <- vst(
dds,
blind = TRUE
)

vst_matrix <- assay(vst_obj)

############################################################

# VERIFY GENE NAMES SURVIVED

############################################################

cat("\nGenes after VST:\n")
print(length(rownames(vst_matrix)))

cat("\nFirst 10 VST genes:\n")
print(head(rownames(vst_matrix),10))

############################################################

# SAVE VST FILE

############################################################

vst_df <- data.frame(
Gene = rownames(vst_matrix),
vst_matrix,
check.names = FALSE
)

write.csv(
vst_df,
"TCGA_LUAD_WhiteBlack_VST.csv",
row.names = FALSE
)

############################################################

# SAVE QC REPORT

############################################################

sink("WhiteBlack_VST_Report.txt")

cat("WHITE-BLACK VST REPORT\n\n")

cat("Genes:\n")
print(nrow(vst_df))

cat("\nSamples:\n")
print(ncol(vst_df)-1)

cat("\nClinical patients:\n")
print(nrow(clinical))

cat("\nRace distribution:\n")
print(table(clinical$Race))

sink()

############################################################

# FINAL OUTPUT

############################################################

cat("\n=====================================\n")
cat("VST COMPLETE\n")
cat("=====================================\n")

cat("\nGenes:\n")
print(nrow(vst_df))

cat("\nSamples:\n")
print(ncol(vst_df)-1)

cat("\nFiles Created:\n")
cat("TCGA_LUAD_WhiteBlack_VST.csv\n")
cat("WhiteBlack_VST_Report.txt\n")
############################################################

# END

############################################################
