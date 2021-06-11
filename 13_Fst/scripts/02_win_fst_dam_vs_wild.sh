#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp
vcf_dir=../../10_FinalVCF/output
vcf_nm=cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered.vcf 
lines_dir=../../sample_info

time $vcftools/vcftools --vcf $vcf_dir/$vcf_nm --weir-fst-pop $lines_dir/dam_lines.txt --weir-fst-pop $lines_dir/euro_wild_boar_samples --fst-window-size 50000 --fst-window-step 25000 --out ../output/fst_win_dam_vs_EWB
