#!/bin/bash

# Files
FILE_PATH="/home/jupyter/DATA/Blood_ATAC/processed/PU1/qtltools/input/cqn"

# Index file
bgzip ${FILE_PATH}/PU1.norm_prop.txt && tabix -p bed ${FILE_PATH}/PU1.norm_prop.txt.gz

# Phenotype PCA
QTLtools pca --bed ${FILE_PATH}/PU1.norm_prop.txt.gz --center --scale --out ${FILE_PATH}/PU1.pheno_pca

# Genotype PCA
QTLtools pca --vcf /home/jupyter/DATA/Blood_ATAC/landerlab-vcf/1000_genomes_vcfs/Waszak_47_samples.chr1.vcf.gz --center --scale --out /home/jupyter/DATA/Blood_ATAC/processed/PU1/qtltools/input/cqn/PU1.geno_pca

# Return covariate corrected count matrix
QTLtools correct --bed quantifications.bed.gz --out quantifications_corrected.bed --cov technical_covars.txt --normal

QTLtools correct --bed processed/PU1/qtltools/input/cqn/PU1.norm_prop.txt.gz --out processed/PU1/qtltools/input/cqn/corrected.bed --cov processed/PU1/qtltools/input/cqn/PU1.covariates_prop.txt --normal
