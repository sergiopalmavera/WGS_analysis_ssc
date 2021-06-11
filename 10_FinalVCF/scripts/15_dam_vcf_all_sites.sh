#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

in_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf

pop_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_dam_lines_all_sites.vcf

echo "# subset main vcf for population (exclude ref alleles)"
$bcftools/bcftools view --samples-file ../../sample_info/dam_lines.txt $in_vcf -o $pop_vcf
printf "\n\n"

echo "# Number of SNPs in pop:"
grep -v '^#' $pop_vcf | wc -l
printf "\n\n"


