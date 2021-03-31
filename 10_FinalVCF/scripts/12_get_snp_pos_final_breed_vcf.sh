#!/bin/bash

pop=$1

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_$pop.vcf

echo "# get SNP for population"
$bcftools/bcftools query -f '%CHROM  %POS\n' $vcf > ${vcf/.vcf/.SNPids}
