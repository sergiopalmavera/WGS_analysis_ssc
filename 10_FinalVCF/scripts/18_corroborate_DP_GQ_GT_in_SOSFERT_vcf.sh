#!/bin/bash

echo "# min DP"
cat ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_DP.tab | tr ' ' \\n | grep . | sort -h | head -1

echo "# max DP"
cat ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_DP.tab | tr ' ' \\n | grep . | sort -h | tail -1

echo "# min GQ"
cat ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_GQ.tab | tr ' ' \\n | grep . | sort -h | head -1

echo "# max GQ"
cat ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_GQ.tab | tr ' ' \\n | grep . | sort -h | tail -1

echo "# unique genotypes"
cat ../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_SOSFERTsamples_GTfiltered_GT.tab | tr ' ' \\n | grep . | sort | uniq
