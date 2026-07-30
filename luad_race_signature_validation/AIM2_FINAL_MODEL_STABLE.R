library(survival)
library(dplyr)
library(survminer)

cat("=====================================\n")
cat("AIM 2.5 - FINAL MODEL (STABLE)\n")
cat("=====================================\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

clinical <- read.csv("TCGA_LUAD_WhiteBlack_Clinical.csv", check.names = FALSE, stringsAsFactors = FALSE)
expr     <- read.csv("TCGA_LUAD_WhiteBlack_VST.csv", check.names = FALSE, stringsAsFactors = FALSE)
genes    <- read.csv("WhiteBlack_LASSO_SelectedGenes.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Build a clean numeric expression matrix
expr_matrix <- as.matrix(expr[, -1])
rownames(expr_matrix) <- expr[, 1]
class(expr_matrix) <- "numeric"

# Match clinical observations perfectly to the expression matrix columns
patient_ids <- substr(colnames(expr_matrix), 1, 12)
clinical_matched <- clinical[match(patient_ids, clinical$PatientID), ]

# ==============================================================================
# 2. QUALITY CONTROL & MISSING DATA REMOVAL
# ==============================================================================

# Ensure all critical survival and covariate data is non-missing
qc_vars <- c("OS_time", "OS_event", "Race", "Age", "Sex")
keep <- complete.cases(clinical_matched[, qc_vars])

clinical_clean <- clinical_matched[keep, , drop = FALSE]
expr_matrix_clean <- expr_matrix[, keep, drop = FALSE]

cat("Patients passing complete-case QC:", nrow(clinical_clean), "\n")

# ==============================================================================
# 3. GENE SELECTION & EXPRESSION TRANSFORMATION
# ==============================================================================

selected_genes <- intersect(genes$Gene, rownames(expr_matrix_clean))
expr_sel <- expr_matrix_clean[selected_genes, , drop = FALSE]

# Transpose so rows are patients and columns are gene features
x_features <- as.data.frame(t(expr_sel))
y_survival <- Surv(clinical_clean$OS_time, clinical_clean$OS_event)

cat("Genes successfully matched:", length(selected_genes), "\n")

# ==============================================================================
# 4. MULTIVARIATE RISK GENERATION (BUGBUSTER VERSION)
# ==============================================================================

# Fit multivariate model strictly on gene expressions
cox_multi <- coxph(
  y_survival ~ .,
  data = x_features,
  control = coxph.control(iter.max = 100)
)

# FIXED: Use predict() to safely extract PI (Prognostic Index / Risk Score).
# This completely bypasses matrix matching errors and natively handles dropped genes!
risk_score <- predict(cox_multi, type = "lp")
clinical_clean$risk_score <- risk_score

# Stratify patients into risk tiers by the median score
median_score <- median(risk_score, na.rm = TRUE)
clinical_clean$risk_group <- ifelse(risk_score >= median_score, "High", "Low")

# Ensure risk group is stored as a factor variable for survival routines
clinical_clean$risk_group <- factor(clinical_clean$risk_group, levels = c("Low", "High"))

# ==============================================================================
# 5. KAPLAN-MEIER ANALYSIS & GRAPH EXPORT
# ==============================================================================

km_fit <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = clinical_clean)

pdf("AIM2_5_KM.pdf", width = 7, height = 6, onefile = FALSE)
print(
  ggsurvplot(
    km_fit,
    data       = clinical_clean,
    pval       = TRUE,
    risk.table = TRUE,
    conf.int   = TRUE,
    palette    = c("#2E9FDF", "#E7B800"),
    title      = "AIM 2.5 Risk Score Survival Stratification"
  )
)
dev.off()
cat("-> Kaplan-Meier Curve PDF generated successfully.\n")

# ==============================================================================
# 6. ADJUSTED CLINICAL COVARIATE TESTING
# ==============================================================================

# Evaluate the independent prognostic power of your signature against clinical metrics
final_model <- coxph(
  Surv(OS_time, OS_event) ~ risk_score + Race + Age + Sex,
  data = clinical_clean
)

# ==============================================================================
# 7. EXPORT DATASETS & COMPREHENSIVE SUMMARIES
# ==============================================================================

# Structure the comprehensive summary framework file
risk_df <- clinical_clean %>%
  select(PatientID, risk_score, risk_group, OS_time, OS_event, Race, Age, Sex)

write.csv(risk_df, "AIM2_5_RiskScores.csv", row.names = FALSE)
write.csv(summary(final_model)$coefficients, "AIM2_5_MultivariateCox.csv")

cat("\n=====================================\n")
cat("ANALYSIS RUN COMPLETE\n")
cat("=====================================\n")
cat("  Total Sample Space :", nrow(clinical_clean), "\n")
cat("  Generated Matrices : AIM2_5_RiskScores.csv, AIM2_5_MultivariateCox.csv\n")
cat("=====================================\n")
