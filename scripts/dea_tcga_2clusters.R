# Do DEA on TCGA data by the clusters "good" and "bad" using DESeq2

rm(list=ls())

# library(RTCGAToolbox)
library(dplyr)
library(DESeq2)

setwd('/Users/xiaokangzhang/github/ESCC_TIME/scripts/')

dir.output <- "/Users/xiaokangzhang/github/ESCC_TIME/output/tcga_escc_DEA/"

# load the original data
load("/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/escc_tcga.RData")

# Load the clustering info
load('/Users/xiaokangzhang/github/ESCC_TIME/output/tcga_escc_best_clustering_based_on_xcell/clustering_tcga.RData')

# Pick out the common patients
id.inter <- intersect(rownames(df.exp.inter), rownames(res.deconv.sample))
res.deconv.sample <- res.deconv.sample[id.inter, ] # just to make sure that these two df have the same patients order
df.exp.inter.round <- round(df.exp.inter[id.inter, ]) # round up float numbers to integers
cluster <- res.deconv.sample[id.inter, "cluster"]

# Do DEA using DESeq2
samples <- rownames(df.exp.inter.round)
group <- factor(cluster, levels = c("good", "bad"))
design <- model.matrix(~group)

## create the DESeqDataSet
colData <- data.frame(samples, group)
dds <- DESeqDataSetFromMatrix(t(df.exp.inter.round), colData, design = design)

# generate normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Filtering
keep <- rowSums(counts(dds)>10) >= length(samples)/2
dds <- dds[keep,]

## perform DEA
dds <- DESeq(dds)

## export the results
res.dea <- results(dds)
res.dea <- res.dea[complete.cases(res.dea), ]  # remove any rows with NA

dea <- as.data.frame(res.dea)
dea <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

deg <- filter(dea, log2FoldChange > log2(1.5) | log2FoldChange < log2(1/1.5) & padj < 0.05)

gene.dea <- rownames(dea)
gene.deg <- rownames(deg)

# pick out the E2F family
index.E2F <- grep(pattern = "^E2F", x = gene.dea)
dea.E2F <- dea[index.E2F, ]

save(dea.E2F, dea, deg, gene.deg, df.exp.inter.round, cluster, file = paste0(dir.output, "/dea_tcga_deseq2.RData"))
write.csv(dea.E2F, file = paste0(dir.output, "/dea_E2F_family.csv"))
write.csv(dea, file = paste0(dir.output, "/dea.csv"))
write.csv(deg, file = paste0(dir.output, "/deg.csv"))
