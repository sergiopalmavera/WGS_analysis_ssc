#!/bin/bash

vcf=cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered.roh

cat ../output/$vcf | grep -E "^# RG|^RG" > ../output/${vcf}.table
