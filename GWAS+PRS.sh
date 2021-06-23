## GWAS
# Extract autosomal SNP
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' hweclean.bim > autosomal_snp.txt
plink --bfile hweclean --extract autosomal_snp.txt --make-bed --out autosomal

# Extract SNP from Nature paper: Common polygenic variation contributes to risk of schizophrenia and bipolar disorder
plink --bfile autosomal --extract SNP_ID.txt --make-bed --out nature

# generate cov file from mds
plink --bfile nature --genome --out genome
plink --bfile nature --read-genome genome.genome --cluster --mds-plot 3 --silent --out mds
awk '{print$1, $2, $4, $5, $6}' mds.mds > covar_mds.txt

# GWAS with covariates
plink --bfile nature --covar covar_mds.txt --logistic --hide-covar --out gwas 

## PRS
# LD based clumping
plink --bfile nature --clump gwas.assoc --clump-r2 0.01 --clump-kb 10 --clump-p2 1 --clump-p1 0.1

# add ORs as PRS
