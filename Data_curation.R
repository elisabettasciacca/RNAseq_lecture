rm(list = ls())
library(TCGAbiolinks)
library(SummarizedExperiment)

library(DESeq2)
getProjectSummary(TCGA-SKCM)
query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
data <- GDCprepare(query)
counts <- assay(data)

query <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
metadata <- clinical.BCRtab.all$clinical_drug_skcm[-(1:2),]


# load data 
# let count ids match metadata ids 
colnames(counts) <- paste(
  unlist(lapply(strsplit(colnames(counts), "-"), `[[`, 1)),
  unlist(lapply(strsplit(colnames(counts), "-"), `[[`, 2)),
  unlist(lapply(strsplit(colnames(counts), "-"), `[[`, 3)),
  sep = "-")
# translate ENSEMBL genes in gene symbols 
counts$genesymbol <- AnnotationDbi::mapIds(x = org.Hs.eg.db::org.Hs.eg.db,
                                               keys = rownames(counts),
                                               keytype = 'ENSEMBL',
                                               column = 'SYMBOL')
counts <- subset(counts, !is.na(counts$genesymbol))
counts <- counts[!(duplicated(counts$genesymbol)), ]
rownames(counts) <- counts$genesymbol
counts$genesymbol <- NULL


metadata <- readRDS("../DEG-Gs/TCGA-SKCM_metadata.RDS")
table(metadata$pharmaceutical_therapy_type)
metadata <- subset(metadata, metadata$pharmaceutical_therapy_type == "Chemotherapy")
table(metadata$treatment_best_response)
metadata <- subset(metadata, metadata$treatment_best_response %in% c("Clinical Progressive Disease",
                                                                     "Complete Response"))

table(metadata$treatment_best_response)


metadata <- metadata[!duplicated(metadata$bcr_patient_barcode), ]
rownames(metadata) <- metadata$bcr_patient_barcode
metadata <- metadata[rownames(metadata) %in% colnames(counts), ]
metadata <- metadata[, c("bcr_patient_barcode", "form_completion_date",
                         "pharmaceutical_therapy_type", 
                         "pharmaceutical_therapy_drug_name",  
                         "treatment_best_response", 
                         "pharmaceutical_tx_started_days_to",
                         "pharmaceutical_tx_ended_days_to")]

metadata$treatment_best_response <- gsub(" ", "_", metadata$treatment_best_response)

counts <- counts[, rownames(metadata)] 
identical(rownames(metadata), colnames(counts))

saveRDS(counts, "TCGA-SKCM_rawCounts.RDS")
saveRDS(metadata, "TCGA-SKCM_metadata.RDS")

library(DESeq2)
counts_rlog <- rlog(as.matrix(counts))
saveRDS(counts_rlog, "TCGA-SKCM_counts_rlog.RDS")
