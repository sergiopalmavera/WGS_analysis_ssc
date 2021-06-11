#!/bin/bash

grep -v '^#' ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_Pietrain.vcf | awk '{print $1 "\t" $2}' > ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_Pietrain_chr_pos.tsv

grep -v '^#' ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_LargeWhite.vcf | awk '{print $1 "\t" $2}' > ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_LargeWhite_chr_pos.tsv
