#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

vcf_path=../../10_FinalVCF/output

vcf_fl=cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf

$vcftools/vcftools --vcf $vcf_path/$vcf_fl --recode --recode-INFO-all --thin 20000 --out ../output/${vcf_fl/.vcf/_20Kb_thinned} 

# Reference post: https://www.biostars.org/p/347796/	
# From the official manual page:
# "--thin <integer> Thin sites so that no two sites are within the specified distance from one another."

