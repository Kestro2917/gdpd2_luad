############################################################
# AIM 3.4
# FUNCTIONAL ENRICHMENT OF GDPD2 SIGNATURE
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

cat("=====================================\n")
cat("AIM 3.4 - FUNCTIONAL ENRICHMENT ANALYSIS\n")
cat("=====================================\n\n")

############################################################
# LOAD DESEQ2 RESULTS
############################################################

res <- read.csv(
  "AIM3_3_GDPD2_DESeq2_AllGenes.csv",
  check.names = FALSE
)

############################################################
# CLEAN DATA
############################################################

res <- res %>%
  filter(
    !is.na(padj)
  )

############################################################
# FILTER SIGNIFICANT GENES
############################################################

sig <- res %>%
  filter(
    padj < 0.05 &
    abs(log2FoldChange) > 1
  )

cat("Significant genes for enrichment:\n")
print(nrow(sig))

############################################################
# GENE SYMBOLS
############################################################

gene_list <- sig$Gene

############################################################
# CONVERT TO ENTREZ IDS
############################################################

gene_map <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

############################################################
# GO ENRICHMENT
############################################################

go_enrich <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_df <- as.data.frame(go_enrich)

############################################################
# KEGG ENRICHMENT
############################################################

kegg_enrich <- enrichKEGG(
  gene = gene_map$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg_df <- as.data.frame(kegg_enrich)

############################################################
# SAVE TABLES
############################################################

write.csv(
  go_df,
  "AIM3_4_GO_Enrichment.csv",
  row.names = FALSE
)

write.csv(
  kegg_df,
  "AIM3_4_KEGG_Enrichment.csv",
  row.names = FALSE
)

############################################################
# PLOT TOP PATHWAYS
############################################################

pdf(
  "AIM3_4_Enrichment_DotPlot.pdf",
  width = 10,
  height = 6
)

if(nrow(go_df) > 0)
{
  print(
    dotplot(go_enrich, showCategory = 15) +
      ggtitle("GO Biological Process - GDPD2 Signature")
  )
}

if(nrow(kegg_df) > 0)
{
  print(
    dotplot(kegg_enrich, showCategory = 15) +
      ggtitle("KEGG Pathways - GDPD2 Signature")
  )
}

dev.off()

############################################################
# SUMMARY REPORT
############################################################

report <- data.frame(
  Metric = c(
    "Genes_Input",
    "GO_Terms_Enriched",
    "KEGG_Pathways_Enriched"
  ),
  Value = c(
    nrow(sig),
    nrow(go_df),
    nrow(kegg_df)
  )
)

write.csv(
  report,
  "AIM3_4_Enrichment_Summary_Report.csv",
  row.names = FALSE
)

############################################################
# FINAL OUTPUT
############################################################

cat("\n=====================================\n")
cat("AIM 3.4 COMPLETE\n")
cat("=====================================\n")

cat("\nGenes used:\n")
print(nrow(sig))

cat("\nGO terms:\n")
print(nrow(go_df))

cat("\nKEGG pathways:\n")
print(nrow(kegg_df))

cat("\n=====================================\n")