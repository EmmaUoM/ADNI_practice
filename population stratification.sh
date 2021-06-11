# make sure to merge files which are in same build
# Extract the variants present in OWN dataset from the hapmap dataset.
awk '{print$2}' hweclean.bim > OWN_SNPs.txt
plink --bfile CEU --extract OWN_SNPs.txt --make-bed --out CEU_COMMON

# Extract the variants present in hapmap dataset from the OWN dataset.
awk '{print$2}' CEU.bim > CEU_SNPs.txt
plink --bfile hweclean --extract CEU_SNPs.txt --recode --make-bed --out OWN_COMMON

## The datasets must have the same build. Change the build hapmap data build.
awk '{print$2,$4}' OWN_COMMON.map > build.txt
# buildhapmap.txt contains one SNP-id and physical position per line.
plink --bfile CEU_COMMON --update-map build.txt --make-bed --out CEU_new_build

## Merge the HapMap and OWN data sets

# Prior to merging OWN data with the HapMap data we want to make sure that the files are mergeable, for this we conduct 3 steps:
# 1) Make sure the reference genome is similar in the HapMap and OWN datasets.
# 2) Resolve strand issues.
# 3) Remove the SNPs which after the previous two steps still differ between datasets. (didn't do since no differ SNPs)

# The following steps are maybe quite technical in terms of commands, but we just compare the two data sets and make sure they correspond.

# 1) set reference genome 
awk '{print$2,$5}' CEU_new_build.bim > ref-list.txt
plink --bfile OWN_COMMON --reference-allele ref-list.txt --make-bed --out OWN-adj

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' CEU_new_build.bim > CEU_new_build_tmp
awk '{print$2,$5,$6}' OWN-adj.bim > OWN-adj_tmp
sort CEU_new_build_tmp OWN-adj_tmp |uniq -u > all_differences.txt

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt

plink --bfile OWN-adj --flip flip_list.txt --reference-allele ref-list.txt --make-bed --out corrected_OWN

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' corrected_OWN.bim > corrected_OWN_tmp
sort CEU_new_build_tmp corrected_OWN_tmp |uniq -u  > uncorresponding_SNPs.txt 
# since nothing in uncorresponding_SNPs.txt, we skip this step

# Merge
plink --bfile CEU_new_build --bmerge corrected_OWN.bed corrected_OWN.bim corrected_OWN.fam --make-bed --out data

# QC
plink --bfile data --geno .05 --make-bed --out dataclean

# MDS
plink --bfile data --genome --out genome
plink --bfile data --read-genome genome.genome --cluster --mds-plot 4 --silent --out mds

# figure
Rscript plot.R


