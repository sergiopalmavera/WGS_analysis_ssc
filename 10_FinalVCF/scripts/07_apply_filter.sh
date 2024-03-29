#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../ref_ensemble95

out_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf

$GATK/gatk SelectVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness.vcf \
	-L ../output/keep_snps_min_12_nonmiss_per_pop.intervals \
	-O ${out_vcf/.vcf/.TMP.vcf}

# some sites could loose the ALT allele after adding missingness. Keep sites having at least one ALT allele count
$bcftools/bcftools view -i 'COUNT(GT="AA")>=1 || COUNT(GT="het")>=1' ${out_vcf/.vcf/.TMP.vcf} -o $out_vcf 

rm ${out_vcf/.vcf/.TMP.vcf}
