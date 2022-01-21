# correlation between E2F7 and immune cells

rm(list=ls())

library(immunedeconv)
library(tibble)

load("./data/esca_tcga/esca_tcga.RData")
# Data included: df.clinical.inter, df.exp.inter, df.merge

dir.out <- "./output/corr_E2F7_immune_cells"
dir.create(dir.out)

# deconvolution of immune cells using quantiseq
res.deconv <- deconvolute(t(df.exp.inter), "quantiseq", tumor = TRUE)
res.deconv <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

df.immune <- as.data.frame(t(res.deconv))
exp.e2f7 <- select(df.exp.inter, c(E2F7))

corr.e2f7.with.immune <- cor(exp.e2f7, df.immune, method = "pearson")

df.corr.e2f7.with.immune <- data.frame(chem = colnames(corr.e2f7.with.immune), pearson = corr.e2f7.with.immune[1, ])

pdf(file = file.path(dir.out, "corr_e2f7_immune_quantiseq.pdf"), height = 5)

ggplot(data = df.corr.e2f7.with.immune, aes(x = chem, y = pearson)) +
  geom_bar(stat = "identity") +
  labs(x = "immune genes", y = "Pearson correlaton", title = "Correlation between E2F7 and immune cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dev.off()

# ====== deconvolution of immune cells using XCELL ======
res.deconv <- deconvolute(t(df.exp.inter), "xcell", tumor = TRUE)
res.deconv <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

df.immune <- as.data.frame(t(res.deconv))
exp.e2f7 <- select(df.exp.inter, c(E2F7))

corr.e2f7.with.immune <- cor(exp.e2f7, df.immune, method = "pearson")

df.corr.e2f7.with.immune <- data.frame(chem = colnames(corr.e2f7.with.immune), pearson = corr.e2f7.with.immune[1, ])

pdf(file = file.path(dir.out, "corr_e2f7_immune_xcell.pdf"), height = 7)

ggplot(data = df.corr.e2f7.with.immune, aes(x = chem, y = pearson)) +
  geom_bar(stat = "identity") +
  labs(x = "immune genes", y = "Pearson correlaton", title = "Correlation between E2F7 and immune cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dev.off()

# ====== deconvolution of immune cells using EPIC ======
res.deconv <- deconvolute(t(df.exp.inter), "epic", tumor = TRUE)
res.deconv <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

df.immune <- as.data.frame(t(res.deconv))
exp.e2f7 <- select(df.exp.inter, c(E2F7))

corr.e2f7.with.immune <- cor(exp.e2f7, df.immune, method = "pearson")

df.corr.e2f7.with.immune <- data.frame(chem = colnames(corr.e2f7.with.immune), pearson = corr.e2f7.with.immune[1, ])

pdf(file = file.path(dir.out, "corr_e2f7_immune_epic.pdf"), height = 5)

ggplot(data = df.corr.e2f7.with.immune, aes(x = chem, y = pearson)) +
  geom_bar(stat = "identity") +
  labs(x = "immune genes", y = "Pearson correlaton", title = "Correlation between E2F7 and immune cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dev.off()
