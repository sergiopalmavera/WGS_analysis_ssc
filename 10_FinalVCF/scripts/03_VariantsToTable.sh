#!/bin/bash

chr=$1

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../ref_ensemble95

$GATK/gatk SelectVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness.vcf \
	-L $chr \
	-O ../TMP/tmp_chr$1.vcf 

$GATK/gatk VariantsToTable \
	-V ../TMP/tmp_chr$1.vcf \
	-O ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_chr$1.tsv \
	-F CHROM -F POS -GF GT

rm ../TMP/tmp_chr$1.vcf
