library("bigsnpr")
snp_readBed("tmp-data/public-data.bed")

obj.bigSNP <- snp_attach("tmp-data/public-data.rds")

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
sex <- obj.bigSNP$fam$sex
pop <- obj.bigSNP$fam$family.ID

big_counts(G, ind.col = 1:10)

# 分割训练集测试集的index
set.seed(1)
ind.train <- sample(nrow(G), 400)
ind.test <- setdiff(rows_along(G), ind.train)

# Population Structure: PCA

svd <- big_randomSVD(G, big_scale())
# plot(svd)

library(ggplot2)

# plot(svd, type = "scores") + aes(color = pop)
# plot(svd, type = "scores", scores = 3:4) + aes(color = pop)
# plot(svd, type = "loadings", loadings = 1:10, coeff = 0.4)

# Association: GWAS

gwas <- big_univLogReg(G, y, covar.train = svd$u, ncores=2)
# p value 
plot(gwas)

# QQplot
plot(gwas, type = "Q-Q") + xlim(1, NA)

# Manhattan plot
snp_manhattan(gwas, CHR, POS, npoints = 20e3) + geom_hline(yintercept = -log10(5e-8), color = "red")

# Polygenic Risk Score






































