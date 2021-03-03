#!/bin/bash

# Sample id & FASTQ links 
sample_id=$1

echo  "# Moving prepared bam file:"
ls -1 ../output/TMP/tmp_$sample_id/*.sorted.dedup.bam
mv ../output/TMP/tmp_$sample_id/*.sorted.dedup.bam ../output/FINAL
printf "\n"


