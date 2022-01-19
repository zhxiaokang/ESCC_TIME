# Deconvolute the TIME composition, looping through all possibilities

rm(list = ls())

library(tibble)
library(immunedeconv)
library(ConsensusClusterPlus)
library(survival)
library(survminer)

setwd("/Users/xiaokangzhang/github/ESCC_TIME/scripts/")

# load the data
dir.data <- "/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/"
dir.out <- "/Users/xiaokangzhang/github/ESCC_TIME/output/esca"
load(paste0(dir.data, "/esca_tcga.RData"))

res.deconv <- deconvolute(t(df.exp.inter), "epic", tumor = TRUE)
res.deconv.epic <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

res.deconv <- deconvolute(t(df.exp.inter), "timer", indications=rep("ESCA", nrow(df.exp.inter)))
res.deconv.timer <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

res.deconv <- deconvolute(t(df.exp.inter), "quantiseq", tumor = TRUE)
res.deconv.quantiseq <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

res.deconv <- deconvolute(t(df.exp.inter), "mcp_counter")
res.deconv.mcp_counter <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

res.deconv <- deconvolute(t(df.exp.inter), "xcell")
res.deconv.xcell <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

set_cibersort_binary("./CIBERSORT.R")
set_cibersort_mat("./LM22.txt")

res.deconv <- deconvolute(t(df.exp.inter), "cibersort")
res.deconv.cibersort <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

# ============ Loop through all the possibilities ============

setwd("/Users/xiaokangzhang/github/ESCC_TIME/output/esca/Consensus_clustering_loop_through/")

deconv.methods <- list(res.deconv.epic, res.deconv.quantiseq, res.deconv.timer, res.deconv.mcp_counter, res.deconv.xcell, res.deconv.cibersort)
deconv.methods.names <- c("res.deconv.epic", "res.deconv.quantiseq", "res.deconv.timer", "res.deconv.mcp_counter", "res.deconv.xcell", "res.deconv.cibersort")

clusterAlgs <- c("hc", "km", "pam")
distances <- c("pearson", "spearman", "euclidean")
innerLinkages <- c("complete", "average")

i <- 1
for (deconv.method in deconv.methods) {
  deconv.method.name <- deconv.methods.names[i]
  i <- i + 1
  res.deconv <- deconv.method
  
  res.deconv.scale <- t(scale(t(res.deconv)))
  
  for (clusterAlg in clusterAlgs) {
    for (distance in distances) {
      if (clusterAlg == "km") {
        if (distance != "euclidean") {
          next
        }
      }
      for (innerLinkage in innerLinkages) {
        cc.folder <- paste("Consensus_clustering", deconv.method.name, clusterAlg, distance, innerLinkage, sep = "_")
        results <- ConsensusClusterPlus(as.matrix(res.deconv.scale), maxK=6, reps=1000, pItem=0.8, pFeature=1, 
                                        title=cc.folder, 
                                        clusterAlg=clusterAlg, distance=distance, innerLinkage=innerLinkage, seed=123456, plot="pdf")
        
        res.deconv.sample <- as.data.frame(t(res.deconv.scale))
        
        for (num.cluster in c(2,3,4,5)) {
          res.deconv.sample$cluster <- unname(results[[num.cluster]][["consensusClass"]])
          
          df.clinical.inter.cluster <- merge(df.clinical.inter, res.deconv.sample, by = "row.names")
          df.clinical.inter.cluster <- column_to_rownames(df.clinical.inter.cluster, colnames(df.clinical.inter.cluster)[1])
          
          df.clinical.inter.cluster$OS_MONTHS <- as.numeric(df.clinical.inter.cluster$OS_MONTHS)
          df.clinical.inter.cluster$OS_STATUS <- 1*(df.clinical.inter.cluster$OS_STATUS=="1:DECEASED")
          
          # options(repr.plot.width=6, repr.plot.height=5.5)
          
          pdf(file = paste0(cc.folder, "/survival_clusters_", as.character(num.cluster), ".pdf"), onefile = FALSE)
          km.fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data=df.clinical.inter.cluster)
          p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)
          
          print(p.surv)
          
          dev.off()
        }  
      }
    }
  }
}
