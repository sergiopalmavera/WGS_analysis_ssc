#!/bin/bash

SnpEff=/home/fb4/palma-vera/FBN_HOME/Tools/SnpEff_5.0e/snpEff
VARS=../../10_FinalVCF/output/
IN=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_private_LW.vcf
OUT=${IN/.vcf/.ann.vcf}

echo "# Stating SnpEff..."
java -Xmx4g -jar $SnpEff/snpEff.jar -v -stats ../output/${IN/.vcf/.summary.html} -csvStats ../output/${IN/.vcf/.summary.csv} Sscrofa11.1.99 $VARS/$IN > ../output/$OUT 
echo "#... SnpEff done"

echo "# Extracting annotation fields..."
cat ../output/$OUT | $SnpEff/scripts/vcfEffOnePerLine.pl | java -Xmx4g -jar $SnpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" > ../output/${OUT/.vcf/.tab}
echo "#... done"


