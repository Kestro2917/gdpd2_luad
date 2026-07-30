############################################################

# WHITE-BLACK COX DIAGNOSTIC

############################################################

library(survival)

cat("\n=====================================\n")
cat("STEP 1: LOAD FILES\n")
cat("=====================================\n")

expr <- read.csv(
"TCGA_LUAD_WhiteBlack_VST.csv",
stringsAsFactors = FALSE,
check.names = FALSE
)

clinical <- read.csv(
"TCGA_LUAD_WhiteBlack_Clinical.csv",
stringsAsFactors = FALSE,
check.names = FALSE
)

cat("\nExpression Dimensions:\n")
print(dim(expr))

cat("\nClinical Dimensions:\n")
print(dim(clinical))

############################################################

# CHECK COLUMN NAMES

############################################################

cat("\n=====================================\n")
cat("STEP 2: COLUMN CHECK\n")
cat("=====================================\n")

print(colnames(expr)[1:5])

############################################################

# FIX FIRST COLUMN

############################################################

colnames(expr)[1] <- "Gene"

cat("\nFirst Column Renamed Successfully\n")

############################################################

# EXTRACT GENE NAMES

############################################################

cat("\n=====================================\n")
cat("STEP 3: GENE EXTRACTION\n")
cat("=====================================\n")

gene_names <- expr$Gene

cat("\nNumber of Genes:\n")
print(length(gene_names))

cat("\nFirst 10 Genes:\n")
print(head(gene_names,10))

############################################################

# REMOVE GENE COLUMN

############################################################

expr <- expr[, -1]

cat("\nExpression Dimensions After Removal:\n")
print(dim(expr))

############################################################

# MATRIX CREATION

############################################################

cat("\n=====================================\n")
cat("STEP 4: MATRIX CREATION\n")
cat("=====================================\n")

expr_matrix <- as.matrix(expr)

rownames(expr_matrix) <- gene_names

mode(expr_matrix) <- "numeric"

cat("\nMatrix Dimensions:\n")
print(dim(expr_matrix))

cat("\nFirst 10 Row Names:\n")
print(head(rownames(expr_matrix),10))

############################################################

# SAMPLE IDS

############################################################

cat("\n=====================================\n")
cat("STEP 5: SAMPLE IDS\n")
cat("=====================================\n")

sample_ids <- colnames(expr_matrix)

patient_ids <- substr(
sample_ids,
1,
12
)

cat("\nSamples:\n")
print(length(sample_ids))

cat("\nPatients:\n")
print(length(patient_ids))

cat("\nFirst 10 Patient IDs:\n")
print(head(patient_ids,10))

############################################################

# CLINICAL MATCHING

############################################################

cat("\n=====================================\n")
cat("STEP 6: CLINICAL MATCHING\n")
cat("=====================================\n")

clinical <- clinical[
match(
patient_ids,
clinical$PatientID
),
]

cat("\nMatched Clinical Dimensions:\n")
print(dim(clinical))

cat("\nMatched Patients:\n")
print(sum(!is.na(clinical$PatientID)))

############################################################

# SURVIVAL CHECK

############################################################

cat("\n=====================================\n")
cat("STEP 7: SURVIVAL VARIABLES\n")
cat("=====================================\n")

print(table(clinical$OS_event))

print(summary(clinical$OS_time))

############################################################

# REMOVE MISSING SURVIVAL

############################################################

keep <- complete.cases(
clinical[, c("OS_time","OS_event")]
)

clinical <- clinical[keep, ]

expr_matrix <- expr_matrix[, keep]

cat("\nPatients Remaining:\n")
print(nrow(clinical))

cat("\nSamples Remaining:\n")
print(ncol(expr_matrix))

############################################################

# TEST ONE GENE

############################################################

cat("\n=====================================\n")
cat("STEP 8: SINGLE GENE COX TEST\n")
cat("=====================================\n")

gene_exp <- as.numeric(
expr_matrix[1, ]
)

cat("\nGene Name:\n")
print(rownames(expr_matrix)[1])

cat("\nGene Length:\n")
print(length(gene_exp))

cat("\nClinical Rows:\n")
print(nrow(clinical))

############################################################

# COX MODEL

############################################################

fit <- coxph(
Surv(OS_time, OS_event) ~ gene_exp,
data = clinical
)

cat("\nCOX MODEL SUCCESSFUL\n")

print(summary(fit))

############################################################

# TEST FIRST 10 GENES

############################################################

cat("\n=====================================\n")
cat("STEP 9: LOOP TEST\n")
cat("=====================================\n")

for(i in 1:10)
{
gene_exp <- as.numeric(
expr_matrix[i, ]
)

fit <- coxph(
Surv(OS_time, OS_event) ~ gene_exp,
data = clinical
)

cat(
"Gene",
i,
"OK\n"
)
}

############################################################

# FINAL

############################################################

cat("\n=====================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("=====================================\n")

cat("\nIf you reached this point,\n")
cat("the data are valid and the error is inside the original Cox script.\n")
############################################################

# END

############################################################
