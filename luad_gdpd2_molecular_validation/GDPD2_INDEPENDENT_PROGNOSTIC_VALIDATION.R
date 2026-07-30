############################################################
# AIM 3.1
# GDPD2 INDEPENDENT PROGNOSTIC VALIDATION
############################################################

library(survival)
library(ggplot2)

cat("=====================================\n")
cat("AIM 3.1 - GDPD2 MULTIVARIABLE COX\n")
cat("=====================================\n\n")

############################################################
# LOAD DATA
############################################################

clinical <- read.csv(
  "TCGA_LUAD_WhiteBlack_Clinical.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- read.csv(
  "TCGA_LUAD_WhiteBlack_VST.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################################################
# BUILD EXPRESSION MATRIX
############################################################

expr_matrix <- as.matrix(expr[, -1])

rownames(expr_matrix) <- expr[, 1]

storage.mode(expr_matrix) <- "numeric"

############################################################
# MATCH PATIENTS
############################################################

patient_ids <- substr(
  colnames(expr_matrix),
  1,
  12
)

clinical <- clinical[
  match(patient_ids,
        clinical$PatientID),
]

############################################################
# EXTRACT GDPD2
############################################################

if(!"GDPD2" %in% rownames(expr_matrix))
{
  stop("GDPD2 not found in expression matrix")
}

clinical$GDPD2 <- as.numeric(
  expr_matrix["GDPD2", ]
)

############################################################
# QUALITY CONTROL
############################################################

required_vars <- c(
  "OS_time",
  "OS_event",
  "GDPD2",
  "Race",
  "Age",
  "Sex"
)

keep <- complete.cases(
  clinical[, required_vars]
)

clinical <- clinical[keep, ]

cat("Patients after QC:", nrow(clinical), "\n\n")

############################################################
# CONVERT VARIABLES
############################################################

clinical$Race <- factor(clinical$Race)
clinical$Sex  <- factor(clinical$Sex)

############################################################
# UNIVARIATE COX
############################################################

uni_cox <- coxph(
  Surv(OS_time, OS_event) ~ GDPD2,
  data = clinical
)

uni_sum <- summary(uni_cox)

uni_df <- data.frame(
  Variable = "GDPD2",
  HR = uni_sum$coefficients[,"exp(coef)"],
  Lower95 = uni_sum$conf.int[,"lower .95"],
  Upper95 = uni_sum$conf.int[,"upper .95"],
  P_value = uni_sum$coefficients[,"Pr(>|z|)"]
)

write.csv(
  uni_df,
  "AIM3_1_GDPD2_UnivariateCox.csv",
  row.names = FALSE
)

############################################################
# MULTIVARIABLE COX
############################################################

multi_cox <- coxph(
  Surv(OS_time, OS_event) ~
    GDPD2 +
    Age +
    Sex +
    Race,
  data = clinical
)

multi_sum <- summary(multi_cox)

multi_df <- data.frame(
  Variable = rownames(multi_sum$coefficients),
  HR = multi_sum$coefficients[,"exp(coef)"],
  Lower95 = multi_sum$conf.int[,"lower .95"],
  Upper95 = multi_sum$conf.int[,"upper .95"],
  P_value = multi_sum$coefficients[,"Pr(>|z|)"]
)

write.csv(
  multi_df,
  "AIM3_1_GDPD2_MultivariateCox.csv",
  row.names = FALSE
)

############################################################
# FOREST PLOT
############################################################

plot_df <- multi_df

pdf(
  "AIM3_1_GDPD2_ForestPlot.pdf",
  width = 8,
  height = 5
)

ggplot(
  plot_df,
  aes(
    x = Variable,
    y = HR,
    ymin = Lower95,
    ymax = Upper95
  )
) +
  geom_pointrange() +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  ylab("Hazard Ratio (95% CI)") +
  xlab("") +
  ggtitle("Multivariable Cox Model")

dev.off()

############################################################
# SUMMARY
############################################################

cat("=====================================\n")
cat("AIM 3.1 COMPLETE\n")
cat("=====================================\n\n")

cat("Univariate GDPD2 P-value:\n")
print(uni_df$P_value)

cat("\nMultivariable GDPD2 P-value:\n")
print(
  multi_df[
    multi_df$Variable == "GDPD2",
    "P_value"
  ]
)

cat("\nFiles Generated:\n")
cat("AIM3_1_GDPD2_UnivariateCox.csv\n")
cat("AIM3_1_GDPD2_MultivariateCox.csv\n")
cat("AIM3_1_GDPD2_ForestPlot.pdf\n")