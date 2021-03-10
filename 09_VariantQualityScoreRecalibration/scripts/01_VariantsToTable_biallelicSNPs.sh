#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
path_vcf=../../08_GenotypeGVCFs/output
in_vcf=cohort_biallelicSNPs.vcf.gz
REF=../../ref_ensemble95

$GATK/gatk IndexFeatureFile -I $path_vcf/$in_vcf
 
 $GATK/gatk VariantsToTable \
 	-V $path_vcf/$in_vcf \
 	-F CHROM -F POS -F FILTER -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR \
 	-O ../output/${in_vcf/.vcf.gz/.table}
