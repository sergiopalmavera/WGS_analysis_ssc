#!/bin/bash

pop=$1

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

in_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_20Kb_thinned.recode.vcf
pop_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_20Kb_thinned.recode.${pop}.vcf
samps=../../sample_info/${pop}_samples

echo "# subset main vcf for population (exclude ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -o $pop_vcf
printf "\n\n"
