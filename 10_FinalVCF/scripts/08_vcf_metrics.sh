#!/bin/bash

PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

vcf_fl=cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf

$GATK/gatk IndexFeatureFile -F ../output/$vcf_fl

echo "# Metrics for final vcf"
java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
	INPUT=../output/$vcf_fl \
	OUTPUT=../output/${vcf_fl/.vcf/} \
	DBSNP=../../ref_ensemble95/variation/sus_scrofa_NoWhiteSpacesInINFO.vcf
