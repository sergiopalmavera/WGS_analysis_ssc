#!/bin/bash

# This scripts sorts alignments, adds read groups and marks duplicates

# Sample id & FASTQ links 
sample_id=$1

# Define bam name
TMP=../output/TMP
tmp_sub=tmp_${sample_id} #temp dir for this sample
bam=$(echo $TMP/$tmp_sub/*.bam)

# Define relative paths 
REF=../../ref_ensemble95/
FINAL=../output/FINAL

# absolute paths
PICARD=~/FBN_HOME/Tools/picard_2.18.11
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin

echo "# Adding read groups to bam files for sample $sample_id"
printf "\n"

echo "## Processing bam: $bam"
printf "\n"

# set read groups
PU=$(basename $bam | cut -d'_' -f2 | sed 's/.bam//g')
LB=lib_$sample_id
SM=$sample_id
PL=ILLUMINA

# Print read groups
echo "## Read groups:"
printf "\n"
echo "PU: $PU"
echo "LB: $LB"
echo "SM: $SM"
echo "PL: $PL"
printf "\n"

############
## START ! #
############

echo "# Hostname (aka node) and date"
hostname
date
printf "\n"

# BAM prep -----------------------------------------------------------

bam_nm=$(basename $bam) # set basename of bam

bam_rg=${bam_nm/.bam/.RG.bam}

echo "## Adding Read Group ..."
time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
	I=$TMP/$tmp_sub/$bam_nm \
	O=$TMP/$tmp_sub/$bam_rg \
	RGLB=$LB \
	RGPL=$PL \
	RGSM=$SM \
	RGPU=$PU
printf "\n"

tmp_sub_picard=${tmp_sub}_picard # tmp file for picard
mkdir $TMP/$tmp_sub/$tmp_sub_picard

bam_rg_sort=${bam_nm/.bam/.RG.sorted.bam}

echo "### Sorting BAM file"
time java -jar $PICARD/picard.jar SortSam \
	I=$TMP/$tmp_sub/$bam_rg \
	O=$TMP/$tmp_sub/$bam_rg_sort \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/$tmp_sub/$tmp_sub_picard
printf "\n"

echo "### Indexing BAM file ..."
samtools index $TMP/$tmp_sub/$bam_rg_sort
printf "\n"

bam_rg_sort_dedup=${bam_nm/.bam/.RG.sorted.dedup.bam}

echo "### Marking duplicates ..."
time java -jar $PICARD/picard.jar MarkDuplicates \
	I=$TMP/$tmp_sub/$bam_rg_sort \
	O=$TMP/$tmp_sub/$bam_rg_sort_dedup \
	M=$FINAL/${bam_rg_sort_dedup/.bam/.metrics.txt} \
	TMP_DIR=$TMP/$tmp_sub/$tmp_sub_picard
	printf "\n"
