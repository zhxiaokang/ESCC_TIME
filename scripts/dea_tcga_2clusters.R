# Do DEA on TCGA data by the clusters "good" and "bad" using DESeq2

rm(list=ls())

# library(RTCGAToolbox)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(stringr)

setwd('/Users/xiaokangzhang/github/ESCC_TIME/scripts/')

dir.output <- "/Users/xiaokangzhang/github/ESCC_TIME/output/tcga_esca_DEA/"

# load the original data
load("/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/esca_tcga.RData")

# Load the clustering info
load('/Users/xiaokangzhang/github/ESCC_TIME/output/tcga_esca_best_clustering_based_on_xcell/clustering_tcga.RData')

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

save(dea.E2F, dea, deg, gene.deg, df.exp.inter.round, cluster, file = file.path(dir.output, "dea_tcga_deseq2.RData"))
write.csv(dea.E2F, file = file.path(dir.output, "dea_E2F_family.csv"))
write.csv(dea, file = file.path(dir.output, "dea.csv"))
write.csv(deg, file = file.path(dir.output, "deg.csv"))

# visualization with Vocano plot and Heatmap
# volcano plot
dea.table.volcano <- dea  # for better volcano plot, 0 FDRs/padj will be changed to a very low value

# change the 0 padj to a low value (100 times smaller than the minumum non-zero value)
padj <- dea.table.volcano$padj
padj.min.non.zero <- min(padj[padj>0])
dea.table.volcano$padj[padj==0] <- padj.min.non.zero/100
## define the range of x-axis and y-axis
log2FC_lim <- c(log2(1/1.5) - 1, log2(1.5) + 1)
padj_lim <- -log10(min(dea.table.volcano$padj)) # NAs already removed from dea.table

fig.volcano <- EnhancedVolcano(dea.table.volcano, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, col = c("grey30", "orange2", "royalblue", "red2"),
                               FCcutoff = log2(1.5), xlim = log2FC_lim, ylim = c(0, padj_lim), title = NULL, subtitle = NULL,
                               # highlight E2F7 gene
                               selectLab = "E2F4")

pdf(file = file.path(dir.output, 'volcano_plot.pdf'), width = 9, height = 8)
print(fig.volcano)
dev.off()

# Heatmap
## draw heatmap
df.heatmap <- df.exp.inter[, gene.deg[1:500]]

palette <- c("#999999", "#377EB8")
palette.group <- str_replace(group, "good", "#999999")
palette.group <- str_replace(palette.group, "bad", "#377EB8")

pdf(file = file.path(dir.output, 'heatmap.pdf'), width = 15, height = 12, title = 'Heatmap using the top features')
heatmap(t(df.heatmap), ColSideColors = palette.group, margins = c(9,5.5), labRow = gene.deg[1:500], cexRow = 1.9, cexCol = 1.9)
legend("topleft", title = 'Group', legend=group, text.font = 15,
       col = palette, fill = palette, cex=1.8)
dev.off()



