#!/bin/bash

# This scripts sorts alignments, adds read groups and marks duplicates

# Sample id 
sample_id=$1

echo "# Final part bam prep (merge->sort->dedup) for sample $sample_id"
printf "\n"

# bam directories
TMP=../output/TMP
tmp_sub=tmp_${sample_id} #temp dir for this sample
mkdir $TMP/${tmp_sub}_picard

# absolute paths modify accordingly
PICARD=~/FBN_HOME/Tools/picard_2.18.11
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin

# Get all bams for this sample
bams=$(ls -1 $TMP/tmp_${sample_id}/*.RG.bam)
echo "## BAMs for this sample"
for bam in $bams; do echo $bam; done
printf "\n"

# Merge bams
bams_proc=$(for bam in $bams; do echo "I=$bam"; done)
#
echo "## Merging and sorting BAMs with RG added"
java -jar $PICARD/picard.jar MergeSamFiles \
	$(echo $bams_proc) \
	O=$TMP/tmp_${sample_id}/$sample_id.RG.merged.sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/${tmp_sub}_picard

echo "## Marking duplicates"
java -jar $PICARD/picard.jar MarkDuplicates \
       I=$TMP/tmp_${sample_id}/$sample_id.RG.merged.sorted.bam \
       O=$TMP/tmp_${sample_id}/$sample_id.RG.merged.sorted.dedup.bam \
       M=../output/FINAL/$sample_id.RG.merged.sorted.dedup.metrics \
       TMP_DIR=$TMP/${tmp_sub}_picard

