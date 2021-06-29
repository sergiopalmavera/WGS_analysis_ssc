#!/bin/bash

fl_dir=../../10_FinalVCF/output

vcf=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf

plink --vcf $fl_dir/$vcf --make-bed --chr 1-18 --out ../output/${vcf/.vcf/}


