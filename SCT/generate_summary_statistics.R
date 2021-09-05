library(bigsnpr)

plink <- download_plink("tmp-data")

# do QC
bed <- snp_plinkQC(plink, prefix.in = "tmp-data/hweclean",
                   geno = 0, maf = 0.05, hwe = 1e-10,
                   extra.options = " --thin 0.1")

# Read PLINK files into a "bigSNP"
# The path to the RDS file that stores the bigSNP object. Note that this function creates one other file which stores the values of the Filebacked Big Matrix.
# You shouldn't read from PLINK files more than once. Instead, use snp_attach to load the "bigSNP" object in any R session from backing files
rds <- snp_readBed(bed)
snp <- snp_attach("tmp-data/hweclean_QC.rds")

# Imputation
G <- snp$genotypes
CHR <- snp$map$chromosome
infos <- snp_fastImpute(G, CHR)
snp$genotypes$code256 <- CODE_IMPUTE_PRED
snp <- snp_save(snp)
big_counts(G, ind.col = 1:10)


counts <- big_counts(G)
# print(summary(counts))

sum(counts[4, ])
counts[, 1:10]

fam <- snp$fam
ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))

snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L

ind_norel <- which(fam$Relationship == "unrel")
maf <- snp_MAF(G, ind_norel)

bed <- snp_writeBed(snp, "tmp-data/temp.bed",
                    ind.row = ind_norel, ind.col = which(maf > 0.05))

rds <- snp_readBed(bed)
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
set.seed(1)

# devtools::install_github("privefl/paper2-PRS/pkg.paper.PRS", force=TRUE)

pheno <- pkg.paper.PRS::get_pheno(G, 0.8, 100)
snp$fam$affection <- pheno$pheno + 1
s <- scale(G[, pheno$set]) %*% pheno$effects + rnorm(nrow(G), sd = 2)
# Fast truncated SVD with initial pruning and that iteratively removes long-range LD regions.
obj.svd <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.1)
plot(obj.svd, type = "scores")

# GWAS
gwas <- big_univLinReg(G, s, covar.train = obj.svd$u, ncores = 1)

plot(gwas)
snp_manhattan(gwas, CHR, POS)
y_hat <- mean(pheno$pheno)
plot(pheno$effects, gwas$estim[pheno$set] * y_hat * (1 - y_hat), pch = 20); abline(0, 1, col = "red")


sumstats <- cbind.data.frame(
  snp$map[-3],
  beta = gwas$estim,
  beta_se = gwas$std.err,
  n_case = sum(pheno$pheno == 1),
  n_control = sum(pheno$pheno == 0),
  p = predict(gwas, log10 = FALSE))

# 写文件
snp_writeBed(snp, "tmp-data/new-data.bed")
saveRDS(pheno, file = "tmp-data/data-pheno.rds")
bigreadr::fwrite2(sumstats, "tmp-data/sumstats.txt")
zip("data-raw/public-data.zip",
    paste0("tmp-data/public-data", c(".bed", ".bim", ".fam", "-pheno.rds", "-sumstats.txt")))




