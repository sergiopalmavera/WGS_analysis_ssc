#!/bin/bash

pop=$1

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

in_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf

pop_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_${pop}_all_sites.vcf

samps=../../sample_info/SOSFERT_${pop}_samples.txt

echo "# subset main vcf for population (exclude ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -o $pop_vcf
printf "\n\n"

echo "# Number of SNPs in pop:"
grep -v '^#' $pop_vcf | wc -l
printf "\n\n"


