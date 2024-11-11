#!/bin/bash

# Run shapeit5 on a multisample vcf file, phasing common SNPs one chromosome at a time, without chunking (for WES).
# When launched saves the results in the ROOT_DIR/ in a new directory named phased_genotypes

# Root directory where results will be saved
ROOT_DIR=/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing
# vcf file containing all samples to be phased
VCF_FILE=/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/vcf/TCGA-BRCA1_all_sample_nochr.bcf.gz

# threads to be used by shapeit5
THREADS=60

mkdir ${ROOT_DIR}/phased_genotypes
mkdir ${ROOT_DIR}/phased_genotypes/phase_common

# COMMON SNP PHASING

cd ${ROOT_DIR}/phased_genotypes/phase_common

for CHR in {1..22}; do

    MAP=${ROOT_DIR}/resources/maps/chr${CHR}.b38.gmap.gz

    OUT=TCGA-BRCA1-chr${CHR}-phase_common.bcf

    ${ROOT_DIR}/tools/shapeit5/phase_common_static \
        --input $VCF_FILE \
        --region $CHR \
        --map $MAP \
        --filter-maf 0.001 \
        --output $OUT \
        --thread $THREADS && \
    bcftools index -f $OUT

done