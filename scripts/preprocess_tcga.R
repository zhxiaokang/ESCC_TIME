# Explore the ESCC/ESCA datasets
## 90% of ESCA are ESCC

library(dplyr)
library(tibble)
library(stringr)

setwd("/Users/xiaokangzhang/github/ESCC_TIME/")

dir.esca.tcga <- "/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/"

setwd(dir.esca.tcga)

# === clinical info ===
data.bcr.clinical.data.patient <- read.table("data_bcr_clinical_data_patient.txt", sep = "\t", header = TRUE)
dim(data.bcr.clinical.data.patient)
# 185 82

data.bcr.clinical.escc <- filter(data.bcr.clinical.data.patient, HISTOLOGICAL_DIAGNOSIS=="Esophagus Squamous Cell Carcinoma")
dim(data.bcr.clinical.escc)
# 96 82

# pick out the interesting cols
df.clinical.escc <- select(data.bcr.clinical.escc, c(PATIENT_ID, TUMOR_STATUS, OS_STATUS, OS_MONTHS, DFS_STATUS, DFS_MONTHS))
df.clinical.escc <- column_to_rownames(df.clinical.escc, "PATIENT_ID")
patient.ids.clinical <- make.names(rownames(df.clinical.escc))
row.names(df.clinical.escc) <- patient.ids.clinical
# 96 5
df.clinical.escc <- filter(df.clinical.escc, OS_MONTHS > 1) # only keep patients with OS_MONTH > 1
patient.ids.clinical <- make.names(rownames(df.clinical.escc))
dim(df.clinical.escc)
# 91 5


# === gene expression ===
data.rna.seq.expression.median <- read.table("data_RNA_Seq_v2_expression_median.txt", sep = "\t", header = TRUE)
dim(data.rna.seq.expression.median)
# 20531   187

data.exp <- select(data.rna.seq.expression.median, -c(Hugo_Symbol, Entrez_Gene_Id))
df.exp <- data.frame(t(data.exp))
# give the gene names to the cols
colnames(df.exp) <- data.rna.seq.expression.median$Hugo_Symbol
# deal with patients
## remove the duplicate patient TCGA.V5.A7RC
index.TCGA.V5.A7RC <- str_which(rownames(df.exp), "TCGA.V5.A7RC")
df.exp <- df.exp[-index.TCGA.V5.A7RC, ]
## remove the number at the end of the patient ID to agree with clinical info
patient.ids.exp <- rownames(df.exp)
## remove the ".01" and ".06" from the ID
patient.ids.exp <- gsub("\\.[0-9]+$", "", patient.ids.exp)
## then re-assign the rownames
rownames(df.exp) <- patient.ids.exp
df.exp <- df.exp[, !duplicated(colnames(df.exp))] # remove duplicated genes
df.exp <- df.exp[, !is.na(colnames(df.exp))]

# intersect two matrices
patient.id.inter <- intersect(patient.ids.clinical, patient.ids.exp)

df.clinical.inter <- df.clinical.escc[patient.id.inter, ]
df.exp.inter <- df.exp[patient.id.inter, ]

df.merge <- merge(df.clinical.inter, df.exp.inter, by = "row.names")
df.merge <- column_to_rownames(df.merge, colnames(df.merge)[1])

save(df.clinical.inter, df.exp.inter, df.merge, file = paste0(dir.esca.tcga, "/escc_tcga.RData"))


