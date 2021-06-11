#!/bin/bash

# Starting from the previously flagged VCF file
# input vcf has been set to no-call if DP<4, GQ<20 and DP > mean DP + 3*SD 
in_vcf=cohort_biallelicSNPs_HardFiltered_WithMissingness.vcf
out_vcf=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

# keep samples sequenced for the SOSFERT project
$bcftools/bcftools view --samples-file ../../sample_info/SOSFERT_samples.txt ../output/$in_vcf | \
	# remove any site with at least one no-call (./.)
	$bcftools/bcftools view -e 'GT[*]="./."' | \
	# keep sites having at least one alternative allele
	# in other words, remove any site where all samples are hom-ref
	$bcftools/bcftools view -i 'GT[*]="alt"' > ../output/$out_vcf

# extract DP and GQ table to corroborate
$bcftools/bcftools query -f '[%GT ]\n' ../output/$out_vcf > ../output/${out_vcf/.vcf/_GT.tab}
$bcftools/bcftools query -f '[%DP ]\n' ../output/$out_vcf > ../output/${out_vcf/.vcf/_DP.tab}
$bcftools/bcftools query -f '[%GQ ]\n' ../output/$out_vcf > ../output/${out_vcf/.vcf/_GQ.tab}
