#!/bin/bash

pop=$1

pop_vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_$pop.vcf

PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

echo "# Making index"
$GATK/gatk IndexFeatureFile -F $pop_vcf
printf "\n\n"

echo "# get metrics for pop"
java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
        INPUT=$pop_vcf \
        OUTPUT=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_$pop \
        DBSNP=../../ref_ensemble95/variation/sus_scrofa_NoWhiteSpacesInINFO.vcf

