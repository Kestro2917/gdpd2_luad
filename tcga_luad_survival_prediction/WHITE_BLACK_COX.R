library(survival)
library(dplyr)

cat("==> SCRIPT STARTED\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

clinical <- read.csv("TCGA_LUAD_WhiteBlack_Clinical.csv", check.names = FALSE)
expr_raw <- read.csv("TCGA_LUAD_WhiteBlack_VST.csv", check.names = FALSE)
cat("-> Files loaded successfully.\n")

# Build expression matrix
expr_matrix <- as.matrix(expr_raw[, -1])
rownames(expr_matrix) <- expr_raw[, 1]
class(expr_matrix) <- "numeric"

# Match clinical data to expression sample IDs
patient_ids <- substr(colnames(expr_matrix), 1, 12)
clinical_matched <- clinical[match(patient_ids, clinical$PatientID), ]
cat("-> Data matrix created and samples aligned.\n")

# Remove missing survival data (QC)
keep <- complete.cases(clinical_matched[, c("OS_time", "OS_event")])
clinical_matched <- clinical_matched[keep, ]
expr_matrix <- expr_matrix[, keep]
cat("-> Survival QC complete.\n")

# Print dataset dimensions to console
cat(
  "\n[Dataset Dimensions]:",
  "\n  Patients :", nrow(clinical_matched),
  "\n  Genes    :", nrow(expr_matrix),
  "\n  Events   :", sum(clinical_matched$OS_event), "\n\n"
)


# ==============================================================================
# 2. COX PROPORTIONAL HAZARDS SCREENING
# ==============================================================================

n_genes <- nrow(expr_matrix)
results_list <- vector("list", n_genes)

cat("==> STARTING COX SCREENING\n")

for (i in seq_len(n_genes)) {
  gene_exp <- expr_matrix[i, ]

  fit <- tryCatch(
    coxph(Surv(OS_time, OS_event) ~ gene_exp, data = clinical_matched),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    s <- summary(fit)
    results_list[[i]] <- data.frame(
      Gene   = rownames(expr_matrix)[i],
      Beta   = s$coef[1, "coef"],
      HR     = s$coef[1, "exp(coef)"],
      Pvalue = s$coef[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }

  if (i %% 5000 == 0) {
    cat("  Processed:", i, "of", n_genes, "genes\n")
  }
}


# ==============================================================================
# 3. PROCESS RESULTS AND ADJUST FDR
# ==============================================================================

cat("\n==> PROCESSING AND FILTERING RESULTS\n")

# Combine results, filter out failed fits, calculate FDR, and sort by P-value
cox_results <- bind_rows(results_list) %>%
  filter(complete.cases(.)) %>%
  mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
  arrange(Pvalue)

# Isolate significant feature sets
sig_p05   <- filter(cox_results, Pvalue < 0.05)
sig_fdr10 <- filter(cox_results, FDR < 0.10)


# ==============================================================================
# 4. SAVE RESULTS AND EXPORT TXT SUMMARY REPORT
# ==============================================================================

# Export CSV tables
write.csv(cox_results, "WhiteBlack_Cox_AllGenes.csv", row.names = FALSE)
write.csv(sig_p05,     "WhiteBlack_Cox_P05_Genes.csv", row.names = FALSE)
write.csv(sig_fdr10,   "WhiteBlack_Cox_FDR10_Genes.csv", row.names = FALSE)

# Generate and export the summary text file
sink("WhiteBlack_Cox_Summary.txt")
cat("WHITE-BLACK COX ANALYSIS SUMMARY REPORT\n")
cat("=========================================\n\n")
cat("Patients Evaluated : ", nrow(clinical_matched), "\n")
cat("Total Genes Tested : ", nrow(cox_results), "\n")
cat("Genes with P < 0.05: ", nrow(sig_p05), "\n")
cat("Genes with FDR < 0.10: ", nrow(sig_fdr10), "\n\n")
cat("TOP 20 SIGNIFICANT GENES:\n")
cat("-----------------------------------------\n")
print(head(cox_results, 20))
sink()


# ==============================================================================
# 5. FINAL TERMINAL DISPLAY
# ==============================================================================

cat(
  "\n======================================",
  "\n  ANALYSIS RUN COMPLETE",
  "\n======================================",
  "\n  Total Genes Tested :", nrow(cox_results),
  "\n  Genes (P < 0.05)   :", nrow(sig_p05),
  "\n  Genes (FDR < 0.10) :", nrow(sig_fdr10),
  "\n--------------------------------------",
  "\n  Generated Outputs:",
  "\n   - WhiteBlack_Cox_AllGenes.csv",
  "\n   - WhiteBlack_Cox_P05_Genes.csv",
  "\n   - WhiteBlack_Cox_FDR10_Genes.csv",
  "\n   - WhiteBlack_Cox_Summary.txt",
  "\n======================================\n"
)
