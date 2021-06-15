#!/bin/bash

gatk=~/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../ref_ensemble95
path_to_intervals=../../00_analysis/data/sosfert_common_private_snps

$gatk/gatk SelectVariants \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-V ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.vcf \
	-L $path_to_intervals/sosfert_private_snps_LW.intervals \
	-O ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_private_LW.vcf
