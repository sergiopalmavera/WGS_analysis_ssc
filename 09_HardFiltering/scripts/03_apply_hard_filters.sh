#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../ref_ensemble95
path_vcf=../../08_GenotypeGVCFs/output
in_vcf=cohort_biallelicSNPs.vcf.gz

$GATK/gatk SelectVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V $path_vcf/$in_vcf \
	-L ../output/hard_filters_final.intervals \
	-O ../output/${in_vcf/.vcf.gz/_HardFiltered.vcf}
