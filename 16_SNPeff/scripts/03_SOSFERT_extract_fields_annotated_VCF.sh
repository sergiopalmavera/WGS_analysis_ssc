#!/bin/bash

SnpEff=/home/fb4/palma-vera/FBN_HOME/Tools/SnpEff/snpEff
vcfEffOnePerLine=~/FBN_HOME/Tools/SnpEff/snpEff/scripts
IN=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.ann.vcf

cat ../output/$IN | $vcfEffOnePerLine/vcfEffOnePerLine.pl | java -Xmx4g -jar $SnpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" > ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.ann.tab




