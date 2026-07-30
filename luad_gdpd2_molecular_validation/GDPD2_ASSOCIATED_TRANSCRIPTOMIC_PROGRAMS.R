############################################################

# AIM 3.3

# GDPD2-ASSOCIATED TRANSCRIPTOMIC PROGRAMS

#

# TOP QUARTILE vs BOTTOM QUARTILE GDPD2

# DESEQ2 ANALYSIS

############################################################

library(DESeq2)

cat("=====================================\n")
cat("AIM 3.3 - GDPD2 TRANSCRIPTOMIC ANALYSIS\n")
cat("=====================================\n\n")

############################################################

# LOAD DATA

############################################################

clinical <- read.csv(
"TCGA_LUAD_WhiteBlack_Clinical.csv",
check.names = FALSE,
stringsAsFactors = FALSE
)

counts <- read.csv(
"TCGA_LUAD_WhiteBlack_Counts_Unique.csv",
check.names = FALSE,
stringsAsFactors = FALSE
)

############################################################

# BUILD COUNT MATRIX

############################################################

gene_names <- counts[,1]

count_matrix <- as.matrix(
counts[, -1]
)

rownames(count_matrix) <- gene_names

storage.mode(count_matrix) <- "integer"

############################################################

# MATCH PATIENTS

############################################################

patient_ids <- substr(
colnames(count_matrix),
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

# EXTRACT GDPD2 EXPRESSION

############################################################

if(!"GDPD2" %in% rownames(count_matrix))
{
stop("GDPD2 not found in count matrix")
}

gdpd2_expr <- as.numeric(
count_matrix["GDPD2", ]
)

############################################################

# BUILD QUARTILES

############################################################

q25 <- quantile(
gdpd2_expr,
0.25,
na.rm = TRUE
)

q75 <- quantile(
gdpd2_expr,
0.75,
na.rm = TRUE
)

group <- rep(
NA,
length(gdpd2_expr)
)

group[gdpd2_expr <= q25] <- "Low"

group[gdpd2_expr >= q75] <- "High"

############################################################

# KEEP ONLY TOP/BOTTOM QUARTILES

############################################################

keep <- !is.na(group)

count_matrix <- count_matrix[, keep]

clinical <- clinical[keep, ]

clinical$GDPD2_Group <- factor(
group[keep],
levels = c("Low","High")
)

############################################################

# FILTER LOW COUNTS

############################################################

dds <- DESeqDataSetFromMatrix(
countData = count_matrix,
colData = clinical,
design = ~ GDPD2_Group
)

dds <- dds[
rowSums(counts(dds)) >= 10,
]

############################################################

# RUN DESEQ2

############################################################

cat("Running DESeq2...\n\n")

dds <- DESeq(dds)

############################################################

# EXTRACT RESULTS

############################################################

res <- results(
dds,
contrast = c(
"GDPD2_Group",
"High",
"Low"
)
)

res_df <- as.data.frame(res)

res_df$Gene <- rownames(res_df)

res_df <- res_df[
order(res_df$padj),
]

############################################################

# SIGNIFICANT GENES

############################################################

sig_genes <- subset(
res_df,
!is.na(padj) &
padj < 0.05
)

############################################################

# SAVE OUTPUTS

############################################################

write.csv(
res_df,
"AIM3_3_GDPD2_DESeq2_AllGenes.csv",
row.names = FALSE
)

write.csv(
sig_genes,
"AIM3_3_GDPD2_DESeq2_FDR05.csv",
row.names = FALSE
)

write.csv(
head(sig_genes,100),
"AIM3_3_GDPD2_Top100.csv",
row.names = FALSE
)

report <- data.frame(
Metric = c(
"Patients_Used",
"GDPD2_High",
"GDPD2_Low",
"Genes_Tested",
"FDR05_Genes"
),
Value = c(
ncol(count_matrix),
sum(clinical$GDPD2_Group=="High"),
sum(clinical$GDPD2_Group=="Low"),
nrow(res_df),
nrow(sig_genes)
)
)

write.csv(
report,
"AIM3_3_GDPD2_Report.csv",
row.names = FALSE
)

############################################################

# SUMMARY

############################################################

cat("=====================================\n")
cat("AIM 3.3 COMPLETE\n")
cat("=====================================\n\n")

cat("Patients Used:\n")
print(ncol(count_matrix))

cat("\nGDPD2 High:\n")
print(sum(clinical$GDPD2_Group=="High"))

cat("\nGDPD2 Low:\n")
print(sum(clinical$GDPD2_Group=="Low"))

cat("\nGenes Tested:\n")
print(nrow(res_df))

cat("\nFDR < 0.05 Genes:\n")
print(nrow(sig_genes))

cat("\nFiles Generated:\n")
cat("AIM3_3_GDPD2_DESeq2_AllGenes.csv\n")
cat("AIM3_3_GDPD2_DESeq2_FDR05.csv\n")
cat("AIM3_3_GDPD2_Top100.csv\n")
cat("AIM3_3_GDPD2_Report.csv\n")

cat("\n=====================================\n")
