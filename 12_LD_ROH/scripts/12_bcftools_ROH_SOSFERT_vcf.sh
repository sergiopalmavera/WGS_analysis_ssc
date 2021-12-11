#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
vcf_path=../../10_FinalVCF/output
vcf_fl=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf
out_fl=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.roh

time $bcftools/bcftools roh --rec-rate 0.7e-8 --AF-file ../output/${vcf_fl/.vcf/.AF.tab.gz} $vcf_path/$vcf_fl -o ../output/$out_fl -S ../../sample_info/SOSFERT_samples.txt

# About recombination rate
# "we assumed a constant recombination rate of 0.7cM/Mb along the chromosomes. " 
# in: Characterization of a haplotype-reference panel for genotyping by low-pass sequencing in Swiss Large White pigs.
# Nosková A, Bhati M, Kadri NK, Crysnanto D, Neuenschwander S, Hofer A, Pausch H.
# BMC Genomics. 2021 Apr 21;22(1):290. doi: 10.1186/s12864-021-07610-5.

# Transform rec rate into physical units, based on: https://www.biostars.org/p/285449/
# 0.7cM => 0.7% chances of recombination in average, per Mb => 0.7/100 = 0.007
# 0.007/Mb => 0.007/1e6 = 7e-9 = 0.7e-8
