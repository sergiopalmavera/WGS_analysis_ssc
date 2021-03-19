#!/bin/bash

fl_dir=../../10_FinalVCF/output

plink --vcf $fl_dir/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.vcf --make-bed --chr 1-18 --out ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed


