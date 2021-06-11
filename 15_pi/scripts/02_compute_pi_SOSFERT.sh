#!/bin/bash

pop=$1

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

vcf_path=../../10_FinalVCF/output

vcf=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_${pop}_all_sites.vcf

$vcftools/vcftools --vcf $vcf_path/$vcf --window-pi 50000 --window-pi-step 25000 --out ../output/${vcf/.vcf/}
