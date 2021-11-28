# Draw survival curve for particular genes

rm(list=ls())

library(survival)
library(survminer)

setwd("/Users/xiaokangzhang/github/ESCC_TIME/scripts/")
dir.output <- "/Users/xiaokangzhang/github/ESCC_TIME/output/E2F_surv/"

load("/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/escc_tcga.RData")
# include data: df.merge, df.clinical.inter, df.exp.inter

# use df.merge
df.surv <- df.merge %>% select(c("OS_STATUS", "OS_MONTHS"), starts_with("E2F"))

# define the function to draw surv curves
time.name = "OS_MONTHS"
event.name = "OS_STATUS"

plot_surv <- function(df.surv, gene.list, time.name = "OS_MONTHS", event.name = "OS_STATUS") {
  ## best cutpoint for each gene
  res.cut <- surv_cutpoint(df.surv, time = time.name, event = event.name, gene.list)
  
  res.cat <- surv_categorize(res.cut)
  
  # draw surv plots for each gene
  for (i in seq(1, length(gene.list))) {
    gene <- gene.list[i]
    forml <- paste0("Surv(", time.name, ", ", event.name, ") ~ ", gene)
    
    km.fit <- survfit(as.formula(forml), data=res.cat)
    p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)
    
    print(p.surv)
  }
}

# prepare the df.surv
df.surv$OS_MONTHS <- as.numeric(df.surv$OS_MONTHS)
df.surv$OS_STATUS <- 1*(df.surv$OS_STATUS=="1:DECEASED")

# draw surv curve for TCGA
gene.list <- names(df.surv)[grep("^E2F", names(df.surv))]
pdf(paste0(dir.output, "/E2F_genes_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = "OS_MONTHS", event.name = "OS_STATUS")
dev.off()
