############################################################
# AIM 2.6 - RACE ROBUSTNESS VALIDATION
############################################################

library(survival)
library(dplyr)
library(survminer)
library(boot)

cat("=====================================\n")
cat("AIM 2.6 - PUBLICATION GRADE ANALYSIS\n")
cat("=====================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

risk_data <- read.csv("AIM2_5_RiskScores.csv", stringsAsFactors = FALSE)

risk_data$OS_time  <- as.numeric(risk_data$OS_time)
risk_data$OS_event <- as.numeric(risk_data$OS_event)

cat("Total patients:", nrow(risk_data), "\n")

# ==============================================================================
# 2. SPLIT BY RACE
# ==============================================================================

white <- risk_data %>% filter(Race == "White")
black <- risk_data %>% filter(Race == "Black")

cat("White:", nrow(white), "Black:", nrow(black), "\n")

# ==============================================================================
# 3. FUNCTION: CINDEX WITH BOOTSTRAP CI
# ==============================================================================

boot_cindex <- function(data, indices) {

  d <- data[indices, ]

  model <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = d)

  return(summary(model)$concordance[1])
}

compute_bootstrap_cindex <- function(df, label) {

  set.seed(123)

  boot_res <- boot(
    data = df,
    statistic = boot_cindex,
    R = 500
  )

  ci <- boot.ci(boot_res, type = "perc")

  c_index <- mean(boot_res$t)

  cat("\n", label, "C-index:", c_index, "\n")

  return(data.frame(
    Race = label,
    C_index = c_index,
    CI_Lower = ci$percent[4],
    CI_Upper = ci$percent[5]
  ))
}

# ==============================================================================
# 4. RUN BOOTSTRAP C-INDEX
# ==============================================================================

white_cindex <- compute_bootstrap_cindex(white, "White")
black_cindex <- compute_bootstrap_cindex(black, "Black")

cindex_results <- rbind(white_cindex, black_cindex)

write.csv(cindex_results, "AIM2_6_Bootstrap_CIndex.csv", row.names = FALSE)

# ==============================================================================
# 5. KM PLOTS (RACE STRATIFIED)
# ==============================================================================

make_km <- function(df, race_name) {

  df$risk_group <- ifelse(
    df$risk_score >= median(df$risk_score, na.rm = TRUE),
    "High",
    "Low"
  )

  fit <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = df)

  pdf(paste0("AIM2_6_KM_", race_name, ".pdf"), width = 7, height = 6)

  print(
    ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = TRUE,
      palette = c("#2E9FDF", "#E7B800"),
      title = paste("AIM 2.6 KM -", race_name)
    )
  )

  dev.off()
}

make_km(white, "White")
make_km(black, "Black")

# ==============================================================================
# 6. RISK SCORE DISTRIBUTION TEST (STATISTICAL COMPARISON)
# ==============================================================================

risk_test <- wilcox.test(
  white$risk_score,
  black$risk_score
)

risk_dist <- data.frame(
  Test = "Wilcoxon Rank Sum",
  Statistic = risk_test$statistic,
  P_value = risk_test$p.value
)

write.csv(risk_dist, "AIM2_6_RiskScore_Distribution_Test.csv", row.names = FALSE)

# ==============================================================================
# 7. RISK SCORE SUMMARY TABLE
# ==============================================================================

summary_table <- risk_data %>%
  group_by(Race) %>%
  summarise(
    Patients = n(),
    Mean_Risk = mean(risk_score),
    Median_Risk = median(risk_score),
    SD_Risk = sd(risk_score)
  )

write.csv(summary_table, "AIM2_6_RiskScore_Summary.csv", row.names = FALSE)

# ==============================================================================
# 8. FINAL SUMMARY
# ==============================================================================

cat("\n=====================================\n")
cat("AIM 2.6 COMPLETE (PUBLICATION GRADE)\n")
cat("=====================================\n")

print(cindex_results)
print(summary_table)
print(risk_dist)

cat("\nFILES GENERATED:\n")
cat("- AIM2_6_Bootstrap_CIndex.csv\n")
cat("- AIM2_6_KM_White.pdf\n")
cat("- AIM2_6_KM_Black.pdf\n")
cat("- AIM2_6_RiskScore_Distribution_Test.csv\n")
cat("- AIM2_6_RiskScore_Summary.csv\n")

cat("=====================================\n")