# Exploration of gene expression, like correlation between interesting genes

# load libs
library(dplyr)
library(ggplot2)
library(corrr)

# Load data
load("./data/esca_tcga/esca_tcga.RData")
dir.out <- "./output/esca"

dir.create(dir.out)

# # ====== gene expression of E2F family genes =======
# df.exp.e2f <- select(df.exp.inter, starts_with("E2F"))
# 
# ave.exp.e2f <- apply(df.exp.e2f, 2, function(x) {log2(sum(x)/length(x))})
# df.ave.exp.e2f <- data.frame(gene = names(ave.exp.e2f), exp = c(ave.exp.e2f))
# 
# ggplot(data = df.ave.exp.e2f, aes(x = gene, y = exp)) +
#   geom_bar(stat = "identity") +
#   scale_x_discrete("Genes") +
#   scale_y_discrete("Expression (log2)")

# ======= correlation between E2F7 and chemokine genes ============
df.chemokine <- select(df.exp.inter, starts_with("CXCL") | starts_with("CCL"))
exp.e2f7 <- select(df.exp.inter, c(E2F7))

corr.e2f7.with.chemokine <- cor(exp.e2f7, df.chemokine, method = "pearson")

df.corr.e2f7.with.chemokine <- data.frame(chem = colnames(corr.e2f7.with.chemokine), pearson = corr.e2f7.with.chemokine[1, ])

pdf(file = file.path(dir.out, "corr_chemokine_e2f7.pdf"))

ggplot(data = df.corr.e2f7.with.chemokine, aes(x = chem, y = pearson)) +
  geom_bar(stat = "identity") +
  labs(x = "Chemokine genes", y = "Pearson correlaton", title = "Correlation between E2F7 and chemokine genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dev.off()

# ======= correlation between E2F7 and interleukins genes ============
df.interleukins <- select(df.exp.inter, starts_with("IL"))
exp.e2f7 <- select(df.exp.inter, c(E2F7))

corr.e2f7.with.interleukins <- cor(exp.e2f7, df.interleukins, method = "pearson")

df.corr.e2f7.with.interleukins <- data.frame(chem = colnames(corr.e2f7.with.interleukins), pearson = corr.e2f7.with.interleukins[1, ])

pdf(file = file.path(dir.out, "corr_interleukins_e2f7.pdf"), height = 12)

ggplot(data = df.corr.e2f7.with.interleukins, aes(x = chem, y = pearson)) +
  geom_bar(stat = "identity") +
  labs(x = "interleukins genes", y = "Pearson correlaton", title = "Correlation between E2F7 and interleukins genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dev.off()


