# get snp list
plink --bfile ADNI1 --maf 0.0000000001 --write-snplist

# retain variants only begin with "rs"
sed -n "/rs/p" plink.snplist -i

plink --bfile ADNI1 --extract plink.snplist --make-bed --out ADNI1onlySNP

# gender check
plink --bfile ADNI1 --check-sex 

Rscript --no-save gender_check.R
# remove individuals missing genotypes greater than 0.05
plink --bfile ADNI1onlySNP --mind 0.05 --make-bed --out cleanINDIVIDUAL

# check missingness
# Investigate missingness per individual and per SNP and make histograms.
plink --bfile ADNI1onlySNP --missing    
# output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.

Rscript hist_miss.R 

# IBS
plink --bfile cleanINDIVIDUAL --cluster

# SNP filting MAF + Missing SNP

plink --bfile cleanINDIVIDUAL --maf 0.01 --geno 0.05 --make-bed --out clean

# Generate a plot of the MAF distribution.
plink --bfile cleanINDIVIDUAL --freq --out MAF_check
Rscript --no-save MAF_check.R

# Hardy Weinberg E
plink --bfile clean --hardy
awk '{ if ($9 <0.000001) print $0 }' plink.hwe>plinkzoomhwe.hwe

Rscript --no-save hwe.R

plink --bfile clean --hwe 1e-6 --hwe-all --make-bed --out hweclean
