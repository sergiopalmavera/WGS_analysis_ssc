#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp
vcf_dir=../../10_FinalVCF/output
vcf_nm=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf 
lines_dir=../../sample_info

time $vcftools/vcftools --vcf $vcf_dir/$vcf_nm --weir-fst-pop $lines_dir/SOSFERT_Pietrain_samples.txt --weir-fst-pop $lines_dir/SOSFERT_LargeWhite_samples.txt --out ../output/fst_win_SOSFERT_PT_vs_LW
