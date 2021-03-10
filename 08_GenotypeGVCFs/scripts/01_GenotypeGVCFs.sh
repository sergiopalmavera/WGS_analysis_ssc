#!/bin/bash

chr=$1

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
#GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

REF=../../ref_ensemble95

input=../../07_CombineGVCFs/output

$GATK/gatk --java-options "-Xmx200G" GenotypeGVCFs \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V $input/cohort_${chr}.g.vcf.gz \
	-O ../output/cohort_${chr}.vcf.gz
