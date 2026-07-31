================================================================================

GDPD2 in Lung Adenocarcinoma - Analysis Code

================================================================================



This repository contains the analysis code for a study investigating whether

self-reported race adds prognostic value beyond tumor transcriptomics in

lung adenocarcinoma (LUAD), and identifying GDPD2 as a prognostic biomarker

using a race-aware discovery pipeline in TCGA-LUAD.



The pipeline includes:

&#x20; - Genome-wide survival screening and LASSO-based feature selection

&#x20; - Race-associated differential expression analysis (DESeq2)

&#x20; - Random survival forest modeling to assess the prognostic value of race

&#x20; - GDPD2-focused survival and pathway enrichment analysis

&#x20; - Immune cell deconvolution (quanTIseq)

&#x20; - External validation



Data source: TCGA-LUAD (https://portal.gdc.cancer.gov/)

