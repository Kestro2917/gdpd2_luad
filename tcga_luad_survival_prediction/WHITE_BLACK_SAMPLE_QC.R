############################################################
# WHITE-BLACK SAMPLE QC
# REMOVE DUPLICATE PATIENTS
############################################################

counts <- read.csv(
  "TCGA_LUAD_WhiteBlack_Counts.csv",
  check.names = FALSE
)

############################################################
# EXTRACT SAMPLE IDS
############################################################

sample_ids <- colnames(counts)[-1]

patient_ids <- substr(
  sample_ids,
  1,
  12
)

sample_type <- substr(
  sample_ids,
  14,
  15
)

############################################################
# KEEP PRIMARY TUMOR ONLY
############################################################

tumor_idx <- which(sample_type == "01")

counts_tumor <- counts[
  ,
  c(
    1,
    tumor_idx + 1
  )
]

############################################################
# REMOVE DUPLICATE PATIENTS
############################################################

sample_ids <- colnames(counts_tumor)[-1]

patient_ids <- substr(
  sample_ids,
  1,
  12
)

keep <- !duplicated(patient_ids)

counts_unique <- counts_tumor[
  ,
  c(
    TRUE,
    keep
  )
]

############################################################
# REPORT
############################################################

cat("\nPatients after deduplication:\n")
print(sum(keep))

cat("\nSamples retained:\n")
print(ncol(counts_unique)-1)

############################################################
# SAVE
############################################################

write.csv(
  counts_unique,
  "TCGA_LUAD_WhiteBlack_Counts_Unique.csv",
  row.names = FALSE
)

cat(
  "\nSaved:\n",
  "TCGA_LUAD_WhiteBlack_Counts_Unique.csv\n"
)