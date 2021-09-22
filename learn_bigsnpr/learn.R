# 从bigsnpr包里面获取示例数据
bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")

library("bigsnpr")

tmpfile <- tempfile()
# 
# rds文件存着bigSNP object
snp_readBed(bedfile, backingfile = tmpfile)

# 在R session里面attach bigSNP object
obj.bigSNP <- snp_attach(paste0(tmpfile,".rds"))

str(obj.bigSNP, max.level = 1)

# 为slot命名别名
G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

# 对前十个SNP check counts
big_counts(G, ind.col = 1:10)

# PCA

# 半数CPU
NCORES <- nb_cores()

ind.excl <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)
ind.keep <- snp_clumping(G, infos.chr = CHR, exclude = ind.excl)

# 获取前10个PC
obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep)

plot(obj.svd)

plot(obj.svd, type = "scores")

# 使用ggplot更改plot的结果
library(ggplot2)

plot(obj.svd, type = "scores") + aes(color = pop <- rep(c("POP1", "POP2", "POP3"), c(143, 167, 207))) + labs(color = "Population")


# GWAS study
# fit 一个logistic model 在phenotype和each SNP separately
# while adding PCs as covariates to each model

y01 <- obj.bigSNP$fam$affection - 1

obj.gwas <- big_univLogReg(G, y01.train = y01, covar.train = obj.svd$u)

snp_qq(obj.gwas)

obj.gwas.gc <- snp_gc(obj.gwas)
snp_qq(obj.gwas.gc)

# 曼哈顿图
snp_manhattan(obj.gwas.gc, infos.chr = CHR, infos.pos = POS)

# Joint PRS

# 分训练集测试集
ind.train <- sample(nrow(G), 400)
ind.test <- setdiff(rows_along(G), ind.train)

# 训练模型
cmsa.logit <- big_spLogReg(X=G, y01.train = y01[ind.train], ind.train = ind.train, 
                           covar.train = obj.svd$u[ind.train, ], alphas = c(1, 0.5, 0.05, 0.001))

# 得到test set的prediction
preds <- predict(cmsa.logit, X=G, ind.row = ind.test, covar.row=obj.svd$u[ind.test, ])

AUC(preds, y01[ind.test])














