############################################################
# AIM 2 STEP 3
# SURVIVAL VALIDATION OF RACE + PROGNOSTIC GENES
# (KM PLOTS + COX + LOG-RANK TEST)
############################################################

library(survival)
library(survminer)
library(dplyr)

cat("=====================================\n")
cat("AIM 2 STEP 3 - SURVIVAL VALIDATION\n")
cat("=====================================\n\n")

############################################################
# LOAD DATA
############################################################

clinical <- read.csv(
  "TCGA_LUAD_WhiteBlack_Clinical.csv",
  check.names = FALSE
)

expr <- read.csv(
  "TCGA_LUAD_WhiteBlack_VST.csv",
  check.names = FALSE
)

overlap_genes <- read.csv(
  "Aim2_Race_Prognostic_Genes.csv",
  check.names = FALSE
)

genes <- unique(overlap_genes$Gene)

cat("Genes to validate:\n")
print(genes)

############################################################
# BUILD EXPRESSION MATRIX
############################################################

gene_names <- expr[,1]
expr_matrix <- as.matrix(expr[,-1])
rownames(expr_matrix) <- gene_names
storage.mode(expr_matrix) <- "numeric"

############################################################
# ALIGN SAMPLES
############################################################

patient_ids <- substr(colnames(expr_matrix), 1, 12)

clinical <- clinical[
  match(patient_ids, clinical$PatientID),
]

keep <- complete.cases(clinical[, c("OS_time", "OS_event")])

clinical <- clinical[keep, ]
expr_matrix <- expr_matrix[, keep]

############################################################
# STORAGE OBJECTS
############################################################

cox_results <- list()
km_results  <- list()

############################################################
# LOOP OVER GENES
############################################################

for (g in genes) {

  if (!(g %in% rownames(expr_matrix))) next

  cat("\nProcessing gene:", g, "\n")

  exp <- expr_matrix[g, ]

  ##########################################################
  # GROUPING (HIGH vs LOW)
  ##########################################################

  median_val <- median(exp, na.rm = TRUE)

  clinical$Group <- ifelse(exp >= median_val, "High", "Low")
  clinical$Group <- factor(clinical$Group, levels = c("Low", "High"))

  ##########################################################
  # COX MODEL
  ##########################################################

  cox_fit <- coxph(
    Surv(OS_time, OS_event) ~ Group,
    data = clinical
  )

  cox_sum <- summary(cox_fit)

  cox_results[[g]] <- data.frame(
    Gene = g,
    HR = cox_sum$coef[1,2],
    Pvalue = cox_sum$coef[1,5]
  )

  ##########################################################
  # KM CURVE + LOG-RANK TEST
  ##########################################################

  km_fit <- survfit(
    Surv(OS_time, OS_event) ~ Group,
    data = clinical
  )

  surv_test <- survdiff(
    Surv(OS_time, OS_event) ~ Group,
    data = clinical
  )

  km_p <- 1 - pchisq(surv_test$chisq, df = 1)

  km_results[[g]] <- data.frame(
    Gene = g,
    KM_Pvalue = km_p
  )

  ##########################################################
  # HIGH-QUALITY KM PLOT
  ##########################################################

  pdf(
    paste0("KM_", g, ".pdf"),
    width = 7,
    height = 6
  )

  print(
    ggsurvplot(
      km_fit,
      data = clinical,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = FALSE,
      palette = c("#2E86C1", "#E74C3C"),
      legend.title = "Expression",
      legend.labs = c("Low", "High"),
      title = paste("Kaplan-Meier Survival:", g),
      xlab = "Time (days)",
      ylab = "Overall Survival Probability",
      risk.table.height = 0.25,
      risk.table.y.text = FALSE,
      ggtheme = theme_minimal(base_size = 14)
    )
  )

  dev.off()
}

############################################################
# MERGE RESULTS
############################################################

cox_df <- bind_rows(cox_results)
km_df  <- bind_rows(km_results)

final_df <- merge(cox_df, km_df, by = "Gene")

############################################################
# SAVE RESULTS
############################################################

write.csv(
  final_df,
  "Aim2_Step3_All_Survival_Results.csv",
  row.names = FALSE
)

sig_df <- final_df %>%
  filter(KM_Pvalue < 0.05)

write.csv(
  sig_df,
  "Aim2_Step3_Significant_Genes.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat("\n=====================================\n")
cat("STEP 3 COMPLETE\n")
cat("=====================================\n")

cat("\nGenes tested:\n")
print(nrow(final_df))

cat("\nSignificant genes (KM P < 0.05):\n")
print(nrow(sig_df))

cat("\n=====================================\n")
cat("KM PLOTS GENERATED AS PDF FILES\n")
cat("=====================================\n")