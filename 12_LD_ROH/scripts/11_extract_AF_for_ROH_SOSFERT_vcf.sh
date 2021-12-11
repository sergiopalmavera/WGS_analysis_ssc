#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

vcf_path=../../10_FinalVCF/output

vcf_fl=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf

$bcftools/bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%AF\n' $vcf_path/$vcf_fl | bgzip -c > ../output/${vcf_fl/.vcf/.AF.tab.gz} && tabix -s1 -b2 -e2 ../output/${vcf_fl/.vcf/.AF.tab.gz}
