# Draw survival curve for particular genes

rm(list=ls())

library(survival)
library(survminer)

dir.output <- "./output/esca"
dir.create(dir.output)

load("./data/esca_tcga/esca_tcga.RData")
# include data: df.merge, df.clinical.inter, df.exp.inter

# use df.merge
df.surv <- df.merge %>% select(c("OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS"), starts_with("E2F"), c("CCL7", "CCL8", "CCR5"))

# define the function to draw surv curves
plot_surv <- function(df.surv, gene.list, time.name = time.name, event.name = event.name) {
  ## best cutpoint for each gene
  res.cut <- surv_cutpoint(df.surv, time = time.name, event = event.name, gene.list)
  
  res.cat <- surv_categorize(res.cut)
  
  # draw surv plots for each gene
  for (i in seq(1, length(gene.list))) {
    gene <- gene.list[i]
    forml <- paste0("Surv(", time.name, ", ", event.name, ") ~ ", gene)
    
    km.fit <- survminer::surv_fit(as.formula(forml), data=res.cat)
    p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)
    
    print(p.surv)
  }
}

# ======== overall survival ==========
time.name = "OS_MONTHS"
event.name = "OS_STATUS"
# prepare the df.surv
df.surv$OS_MONTHS <- as.numeric(df.surv$OS_MONTHS)
df.surv$OS_STATUS <- 1*(df.surv$OS_STATUS=="1:DECEASED")

# draw OS surv curve for E2F genes
gene.list <- names(df.surv)[grep("^E2F", names(df.surv))]
pdf(paste0(dir.output, "/E2F_genes_OS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()

# draw OS surv curve for genes CCL7, CCL8
gene.list <- names(df.surv)[grep("^CCL", names(df.surv))]
pdf(paste0(dir.output, "/CCL_genes_OS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()

# draw OS surv curve for genes CCR5
gene.list <- names(df.surv)[grep("^CCR5", names(df.surv))]
pdf(paste0(dir.output, "/CCR5_gene_OS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()


# ========= DFS/PFS ==========
time.name = "DFS_MONTHS"
event.name = "DFS_STATUS"
# prepare the df.surv
## remove the patients without DFS info
df.surv <- filter(df.surv, DFS_MONTHS > 0.1)
df.surv$DFS_MONTHS <- as.numeric(df.surv$DFS_MONTHS)
df.surv$DFS_STATUS <- 1*(df.surv$DFS_STATUS=="1:Recurred/Progressed")

# draw DFS surv curve for E2F genes
gene.list <- names(df.surv)[grep("^E2F", names(df.surv))]
pdf(paste0(dir.output, "/E2F_genes_DFS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()

# draw DFS surv curve for CCL7 and CCL8
gene.list <- names(df.surv)[grep("^CCL", names(df.surv))]
pdf(paste0(dir.output, "/CCL_genes_DFS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()

# draw DFS surv curve for CCR5
gene.list <- names(df.surv)[grep("^CCR5", names(df.surv))]
pdf(paste0(dir.output, "/CCR5_gene_DFS_surv.pdf"))
plot_surv(df.surv, gene.list, time.name = time.name, event.name = event.name)
dev.off()

