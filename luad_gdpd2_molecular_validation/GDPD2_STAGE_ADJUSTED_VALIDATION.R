############################################################
# AIM 3.2
# GDPD2 STAGE-ADJUSTED VALIDATION
############################################################

library(survival)
library(ggplot2)

cat("=====================================\n")
cat("AIM 3.2 - STAGE ADJUSTED GDPD2 MODEL\n")
cat("=====================================\n\n")

############################################################
# LOAD FILES
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
  match(
    patient_ids,
    clinical$PatientID
  ),
]

############################################################
# EXTRACT GDPD2
############################################################

clinical$GDPD2 <- as.numeric(
  expr_matrix["GDPD2", ]
)

############################################################
# DETECT STAGE COLUMN
############################################################

stage_col <- grep(
  "stage",
  colnames(clinical),
  ignore.case = TRUE,
  value = TRUE
)[1]

cat("Stage column detected:\n")
print(stage_col)

############################################################
# QC
############################################################

required <- c(
  "OS_time",
  "OS_event",
  "GDPD2",
  "Age",
  "Sex",
  "Race",
  stage_col
)

keep <- complete.cases(
  clinical[, required]
)

clinical <- clinical[
  keep,
]

cat("\nPatients after QC:\n")
print(nrow(clinical))

############################################################
# FACTORS
############################################################

clinical$Race <- factor(clinical$Race)

clinical$Sex <- factor(clinical$Sex)

clinical[[stage_col]] <- factor(
  clinical[[stage_col]]
)

############################################################
# COX MODEL
############################################################

formula_text <- paste(
  "Surv(OS_time, OS_event) ~ GDPD2 + Age + Sex + Race +",
  stage_col
)

cox_model <- coxph(
  as.formula(formula_text),
  data = clinical
)

cox_sum <- summary(cox_model)

############################################################
# RESULTS TABLE
############################################################

results <- data.frame(
  Variable =
    rownames(cox_sum$coefficients),

  HR =
    cox_sum$coefficients[, "exp(coef)"],

  Lower95 =
    cox_sum$conf.int[, "lower .95"],

  Upper95 =
    cox_sum$conf.int[, "upper .95"],

  P_value =
    cox_sum$coefficients[, "Pr(>|z|)"],

  stringsAsFactors = FALSE
)

write.csv(
  results,
  "AIM3_2_GDPD2_StageAdjustedCox.csv",
  row.names = FALSE
)

############################################################
# PUBLICATION FOREST PLOT
############################################################

plot_df <- results

plot_df$Label <- plot_df$Variable

plot_df$Label[
  plot_df$Variable == "GDPD2"
] <- "GDPD2 Expression"

plot_df$Label[
  grepl("Sex", plot_df$Variable)
] <- "Male vs Female"

plot_df$Label[
  grepl("Race", plot_df$Variable)
] <- "White vs Black"

pdf(
  "AIM3_2_GDPD2_StageAdjusted_ForestPlot.pdf",
  width = 8,
  height = 5
)

ggplot(
  plot_df,
  aes(
    x = Label,
    y = HR,
    ymin = Lower95,
    ymax = Upper95
  )
) +
  geom_pointrange(
    size = 0.6
  ) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Stage-Adjusted Cox Model",
    y = "Hazard Ratio (95% CI)",
    x = NULL
  )

dev.off()

############################################################
# GDPD2 RESULT
############################################################

gdpd2_row <- results[
  results$Variable == "GDPD2",
]

cat("\n=====================================\n")
cat("GDPD2 RESULT\n")
cat("=====================================\n")

print(gdpd2_row)

cat("\nFiles Generated:\n")
cat("AIM3_2_GDPD2_StageAdjustedCox.csv\n")
cat("AIM3_2_GDPD2_StageAdjusted_ForestPlot.pdf\n")

cat("\n=====================================\n")
cat("AIM 3.2 COMPLETE\n")
cat("=====================================\n")