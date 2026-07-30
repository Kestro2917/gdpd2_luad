library(tidyverse)
library(immunedeconv)
library(quantiseqr)

# ==============================================================================
# 1. LOAD EXPRESSION DATA
# ==============================================================================

expr <- read.csv(
  "TCGA_LUAD_WhiteBlack_VST.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw VST matrix dimensions:", dim(expr), "\n")

# ==============================================================================
# 2. STERILE NUMERIC MATRIX TRANSFORMATION
# ==============================================================================

expr_matrix <- as.matrix(expr)
storage.mode(expr_matrix) <- "numeric"
attr(expr_matrix, "class") <- NULL

# Handle missing data metrics safely
expr_matrix[is.na(expr_matrix) | is.nan(expr_matrix) | is.infinite(expr_matrix)] <- 0

# Convert VST back to linear scale (2^x) as required by quanTIseq
cat("Converting VST values back to linear scale for quanTIseq...\n")
expr_linear_matrix <- 2^expr_matrix

# ==============================================================================
# 3. DIRECT RUN VIA THE QUANTISEQR ENGINE
# ==============================================================================

cat("\nRunning quanTIseq Deconvolution Engine directly...\n")

quantiseq_raw_output <- quantiseqr::run_quantiseq(
  expression_data = expr_linear_matrix,
  is_arraydata    = FALSE,
  is_tumordata    = TRUE,
  scale_mRNA      = TRUE
)

# Convert raw data matrix to a structured baseline data frame for saving
quantiseq_results <- as.data.frame(quantiseq_raw_output) %>%
  rownames_to_column(var = "SampleID")

cat("Output dimensions calculated:", dim(quantiseq_results), "\n")

# Save raw output table
write.csv(quantiseq_results, "Aim4_1_quanTIseq_raw_results.csv", row.names = FALSE)

# ==============================================================================
# 4. TRANSFORMATION CODE (FIXED: REMOVED NON-EXISTENT COLUMN DROP)
# ==============================================================================

# Because quantiseq output rows are ALREADY samples and columns are cell types,
# we simply assign the row names to SampleID without trying to drop cell_type.
immune_df <- as.data.frame(quantiseq_raw_output) %>%
  rownames_to_column(var = "SampleID")

write.csv(immune_df, "Aim4_1_quanTIseq_immune_matrix.csv", row.names = FALSE)

# ==============================================================================
# 5. SUMMARY STATS GENERATION
# ==============================================================================

summary_stats <- immune_df %>%
  dplyr::select(-SampleID) %>%
  summarise(across(
    everything(),
    list(
      mean   = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      sd     = ~sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))

write.csv(summary_stats, "Aim4_1_quanTIseq_summary.csv", row.names = FALSE)

cat("\n===============================================\n")
cat(" AIM 4.1 COMPLETED SUCCESSFULLY\n")
cat("===============================================\n")
