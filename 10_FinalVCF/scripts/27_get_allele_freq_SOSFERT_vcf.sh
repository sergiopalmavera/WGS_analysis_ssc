#!/bin/bash

vcftools=~/FBN_HOME/Tools/vcftools_0.1.13/cpp/
vcf=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf

$vcftools/vcftools --vcf $vcf --freq --out ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered
