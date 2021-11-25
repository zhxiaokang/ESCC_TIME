# Pick the best clustering set

rm(list=ls())

library(immunedeconv)
library(ConsensusClusterPlus)
library(tibble)
library(survival)
library(survminer)

setwd('/Users/xiaokangzhang/github/ESCC_TIME/scripts')

load("/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/escc_tcga.RData")
# Data included: df.clinical.inter, df.exp.inter, df.merge

dir.output <- "/Users/xiaokangzhang/github/ESCC_TIME/output/tcga_escc_best_clustering_based_on_xcell/"

# Pick the best set
res.deconv <- deconvolute(t(df.exp.inter), "xcell")
res.deconv <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

# Scale: z-score transformation --- 0 mean 1 deviation
res.deconv.scale <- t(scale(t(res.deconv)))

# Clustering
results <- ConsensusClusterPlus(as.matrix(res.deconv.scale), maxK=6, reps=1000, pItem=0.8, pFeature=1, title="Consensus_clustering", 
                                clusterAlg="pam", distance="pearson", seed=123456, plot=NULL)

# ======== Survival analysis with 3 clusters (but merge 2 clusters) =============
res.deconv.sample <- as.data.frame(t(res.deconv.scale))
res.deconv.sample$cluster <- unname(results[[3]][["consensusClass"]])

# Merge cluster 1 and 2
cluster <- unname(results[[3]][["consensusClass"]])
cluster[which(cluster == "1" | cluster == "3")] <- "bad"
cluster[which(cluster == "2")] <- "good"

res.deconv.sample$cluster <- cluster

df.clinical.inter.cluster <- merge(df.clinical.inter, res.deconv.sample, by = "row.names")
df.clinical.inter.cluster <- column_to_rownames(df.clinical.inter.cluster, colnames(df.clinical.inter.cluster)[1])

df.clinical.inter.cluster$OS_MONTHS <- as.numeric(df.clinical.inter.cluster$OS_MONTHS)
df.clinical.inter.cluster$OS_STATUS <- 1*(df.clinical.inter.cluster$OS_STATUS=="1:DECEASED")

# options(repr.plot.width=6, repr.plot.height=5.5)

km.fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data=df.clinical.inter.cluster)
p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)

pdf(paste0(dir.output, "/survival_curve_between_2_clusters.pdf"))
print(p.surv)
dev.off()

save(results, res.deconv, res.deconv.scale, res.deconv.sample, p.surv, file = paste0(dir.output, "/clustering_tcga.RData"))

