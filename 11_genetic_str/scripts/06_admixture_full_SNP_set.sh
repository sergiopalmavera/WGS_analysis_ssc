#!/bin/bash

admixture=/home/fb4/palma-vera/FBN_HOME/Tools/ADMIXTURE/admixture_linux-1.3.0

FL=../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.bed

$admixture/admixture $FL 2
$admixture/admixture $FL 3
$admixture/admixture $FL 4
$admixture/admixture $FL 5
$admixture/admixture $FL 6

# move results to output
mv *.P *.Q ../output
