#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
in_vcf=cohort.vcf.gz
out_vcf=cohort_biallelicSNPs.vcf.gz
REF=../../ref_ensemble95

$GATK/gatk SelectVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V ../output/$in_vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	-O ../output/$out_vcf
