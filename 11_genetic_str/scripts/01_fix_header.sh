#!/bin/bash

# Line 17 in original header had a problem parsing the date. Changing month to Maerz solved the problem.

# New header was created by making an temp vcf (first 1000 lines), then line was fixed as above and header was extracted as the new header for replacing.


# Picard changes dots to commas in QUAL!!!!
#picard=~/FBN_HOME/Tools/picard_2.18.11/
#java -jar $picard/picard.jar FixVcfHeader \
#	I=../../10_FinalVCF/output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf \
#	O=../../10_FinalVCF/output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.vcf \
#	HEADER=../../10_FinalVCF/output/fixed_header_line_17_fixed.vcf

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin/
$bcftools/bcftools reheader -h ../../10_FinalVCF/output/fixed_header_line_17_fixed.vcf ../../10_FinalVCF/output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf > ../../10_FinalVCF/output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.vcf
