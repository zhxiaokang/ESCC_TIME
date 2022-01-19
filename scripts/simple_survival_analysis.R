# simple Survival analysis

load("/Users/xiaokangzhang/github/ESCC_TIME/data/esca_tcga/esca_tcga.RData")

df.e2f7 <- select(df.merge, c(OS_MONTHS, OS_STATUS, DFS_MONTHS, DFS_STATUS, E2F7))

df.e2f7.num <- df.e2f7
df.e2f7.num$DFS_MONTHS <- as.numeric(df.e2f7.num$DFS_MONTHS)
df.e2f7.num$DFS_STATUS <- 1*(df.e2f7.num$DFS_STATUS=="1:Recurred/Progressed")
res.cut <- surv_cutpoint(df.e2f7.num, time="DFS_MONTHS", event="DFS_STATUS", "E2F7")

res.cat <- surv_categorize(res.cut)

km.fit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ E2F7, data = res.cat)
ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)

df.e2f7.num <- df.e2f7
df.e2f7.num$OS_MONTHS <- as.numeric(df.e2f7.num$OS_MONTHS)
df.e2f7.num$OS_STATUS <- 1*(df.e2f7.num$OS_STATUS=="1:DECEASED")
res.cut <- surv_cutpoint(df.e2f7.num, time="OS_MONTHS", event="OS_STATUS", "E2F7")

res.cat <- surv_categorize(res.cut)

km.fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ E2F7, data = res.cat)
ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)
