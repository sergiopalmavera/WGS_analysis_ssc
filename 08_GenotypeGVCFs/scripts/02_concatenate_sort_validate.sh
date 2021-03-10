#!/bin/bash

BCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=~/FBN_HOME/Tools/gatk-4.1.5.0
REF=../../ref_ensemble95

out_nm=cohort.vcf.gz

echo "# concatenate and sort ..."
time $BCFTOOLS/bcftools concat ../output/*.vcf.gz -Ou | $BCFTOOLS/bcftools sort -Oz -o ../output/$out_nm
printf "\n"

echo "# Indexing ..."
$GATK/gatk IndexFeatureFile -I ../output/$out_nm
printff "\n"

echo "# Validate vcf format only ..."
$GATK/gatk ValidateVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V ../output/$out_nm \
	--validation-type-to-exclude ALL
