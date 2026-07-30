############################################################
# CREATE WHITE-BLACK COHORT
############################################################

clinical <- read.csv(
  "TCGA_LUAD_Aim1_Clinical_Final.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- read.csv(
  "TCGA_LUAD_STAR_Counts_Clean.csv",
  check.names = FALSE
)

############################################################
# KEEP WHITE + BLACK ONLY
############################################################

clinical <- clinical[
  clinical$Race %in% c("White","Black"),
]

cat("White + Black patients:\n")
print(table(clinical$Race))

############################################################
# EXTRACT PATIENT IDS FROM COUNTS
############################################################

sample_ids <- colnames(counts)[-1]

patient_ids <- substr(
  sample_ids,
  1,
  12
)

############################################################
# MATCH PATIENTS
############################################################

keep_patients <- intersect(
  clinical$PatientID,
  patient_ids
)

clinical <- clinical[
  clinical$PatientID %in% keep_patients,
]

############################################################
# KEEP MATCHED SAMPLES
############################################################

keep_cols <- c(
  TRUE,
  patient_ids %in% keep_patients
)

counts <- counts[, keep_cols]

############################################################
# SAVE
############################################################

write.csv(
  clinical,
  "TCGA_LUAD_WhiteBlack_Clinical.csv",
  row.names = FALSE
)

write.csv(
  counts,
  "TCGA_LUAD_WhiteBlack_Counts.csv",
  row.names = FALSE
)

############################################################
# REPORT
############################################################

sink("WhiteBlack_Cohort_Report.txt")

cat("White-Black Cohort Report\n\n")

cat("Clinical patients:\n")
print(nrow(clinical))

cat("\nRace distribution:\n")
print(table(clinical$Race))

cat("\nExpression samples:\n")
print(ncol(counts)-1)

sink()

cat("DONE\n")