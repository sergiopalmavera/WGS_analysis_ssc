#!/bin/bash

SnpEff=/home/fb4/palma-vera/FBN_HOME/Tools/SnpEff_5.0e/snpEff
VARS=../../10_FinalVCF/output/
IN=cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf
OUT=${IN/.vcf/.ann.vcf}

echo "# Stating SnpEff..."
java -Xmx4g -jar $SnpEff/snpEff.jar -v -csvStats ../output/${IN/.vcf/.summary.csv} Sscrofa11.1.99 $VARS/$IN > ../output/$OUT 
echo "#... SnpEff done"
