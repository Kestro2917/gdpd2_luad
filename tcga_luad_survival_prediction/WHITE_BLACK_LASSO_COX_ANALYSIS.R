############################################################

# WHITE-BLACK LASSO COX ANALYSIS

# INPUT:

# TCGA_LUAD_WhiteBlack_Clinical.csv

# TCGA_LUAD_WhiteBlack_VST.csv

# WhiteBlack_Cox_FDR10_Genes.csv

#

# OUTPUT:

# WhiteBlack_LASSO_SelectedGenes.csv

# WhiteBlack_LASSO_RiskScores.csv

# WhiteBlack_LASSO_Report.csv

# WhiteBlack_LASSO_Top20Genes.csv

############################################################

library(glmnet)
library(survival)

set.seed(123)

cat("=====================================\n")
cat("WHITE-BLACK LASSO ANALYSIS\n")
cat("=====================================\n\n")

############################################################

# LOAD FILES

############################################################

clinical <- read.csv(
"TCGA_LUAD_WhiteBlack_Clinical.csv",
check.names = FALSE
)

expr <- read.csv(
"TCGA_LUAD_WhiteBlack_VST.csv",
check.names = FALSE
)

cox_genes <- read.csv(
"WhiteBlack_Cox_FDR10_Genes.csv",
check.names = FALSE
)

cat("Files loaded successfully\n\n")

############################################################

# BUILD EXPRESSION MATRIX

############################################################

gene_names <- expr[, 1]

expr_matrix <- as.matrix(
expr[, -1]
)

rownames(expr_matrix) <- gene_names

storage.mode(expr_matrix) <- "numeric"

cat("Expression matrix dimensions:\n")
print(dim(expr_matrix))

############################################################

# MATCH PATIENTS

############################################################

patient_ids <- substr(
colnames(expr_matrix),
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

# REMOVE MISSING SURVIVAL DATA

############################################################

keep <- complete.cases(
clinical[, c("OS_time", "OS_event")]
)

clinical <- clinical[keep, ]

expr_matrix <- expr_matrix[, keep]

cat("\nPatients after QC:\n")
print(nrow(clinical))

############################################################

# EXTRACT FDR SIGNIFICANT GENES

############################################################

candidate_genes <- intersect(
cox_genes$Gene,
rownames(expr_matrix)
)

expr_subset <- expr_matrix[
candidate_genes,
]

cat("\nCandidate genes:\n")
print(length(candidate_genes))

############################################################

# BUILD LASSO INPUT

############################################################

x <- t(expr_subset)

y <- Surv(
clinical$OS_time,
clinical$OS_event
)

cat("\nLASSO matrix dimensions:\n")
print(dim(x))

############################################################

# CROSS-VALIDATED LASSO

############################################################

cat("\nRunning cross-validation...\n")

cvfit <- cv.glmnet(
x = x,
y = y,
family = "cox",
alpha = 1,
nfolds = 10,
standardize = TRUE,
maxit = 100000
)

############################################################

# FINAL MODEL

############################################################

cat("Building final model...\n")

lasso_model <- glmnet(
x = x,
y = y,
family = "cox",
alpha = 1,
lambda = cvfit$lambda.min,
standardize = TRUE,
maxit = 100000
)

############################################################

# EXTRACT NON-ZERO COEFFICIENTS

############################################################

coef_matrix <- as.matrix(
coef(lasso_model)
)

selected <- coef_matrix[
coef_matrix[, 1] != 0,
,
drop = FALSE
]

selected_genes <- rownames(selected)

cat("\nSelected genes:\n")
print(length(selected_genes))

############################################################

# SAVE SELECTED GENES

############################################################

lasso_results <- data.frame(
Gene = rownames(selected),
Coefficient = selected[, 1],
stringsAsFactors = FALSE
)

lasso_results <- lasso_results[
order(
abs(lasso_results$Coefficient),
decreasing = TRUE
),
]

write.csv(
lasso_results,
"WhiteBlack_LASSO_SelectedGenes.csv",
row.names = FALSE
)

############################################################

# SAVE TOP 20 GENES

############################################################

write.csv(
head(lasso_results, 20),
"WhiteBlack_LASSO_Top20Genes.csv",
row.names = FALSE
)

############################################################

# CALCULATE RISK SCORES

############################################################

risk_score <- predict(
lasso_model,
newx = x,
type = "link"
)

risk_df <- data.frame(
PatientID = clinical$PatientID,
RiskScore = as.numeric(risk_score),
OS_time = clinical$OS_time,
OS_event = clinical$OS_event,
stringsAsFactors = FALSE
)

write.csv(
risk_df,
"WhiteBlack_LASSO_RiskScores.csv",
row.names = FALSE
)

############################################################

# SUMMARY REPORT

############################################################

report_df <- data.frame(
Metric = c(
"Patients",
"Candidate_Genes",
"Selected_Genes",
"Lambda_Min"
),
Value = c(
nrow(clinical),
length(candidate_genes),
length(selected_genes),
cvfit$lambda.min
),
stringsAsFactors = FALSE
)

write.csv(
report_df,
"WhiteBlack_LASSO_Report.csv",
row.names = FALSE
)

############################################################

# FINAL SUMMARY

############################################################

cat("\n=====================================\n")
cat("LASSO ANALYSIS COMPLETE\n")
cat("=====================================\n")

cat("\nPatients:\n")
print(nrow(clinical))

cat("\nCandidate genes:\n")
print(length(candidate_genes))

cat("\nSelected genes:\n")
print(length(selected_genes))

cat("\nLambda Min:\n")
print(cvfit$lambda.min)

cat("\nFiles Created:\n")
cat("WhiteBlack_LASSO_SelectedGenes.csv\n")
cat("WhiteBlack_LASSO_Top20Genes.csv\n")
cat("WhiteBlack_LASSO_RiskScores.csv\n")
cat("WhiteBlack_LASSO_Report.csv\n")

cat("\n=====================================\n")
