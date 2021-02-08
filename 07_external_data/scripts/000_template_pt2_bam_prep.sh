#!/bin/bash

# This scripts sorts alignments, adds read groups and marks duplicates

# Sample id & FASTQ links 
sample_id=

# Manually define bam 
base_nm=
bam=${base_nm}.bam

# Define read groups
# Read groups added according to: https://www.biostars.org/p/47487/
# Platform unit as: instrument_RunNumber_flowwcellID_lane 
# According to https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.html
PU=YYYYYY
LB=lib_$sample_id
SM=$sample_id
PL=ILLUMINA

# Define relative paths 
REF=../../ref_ensemble95/
TMP=../output/TMP
FINAL=../output/FINAL
tmp_sub=tmp_${sample_id} #temp dir for this sample

# absolute paths modify accordingly
PICARD=~/FBN_HOME/Tools/picard_2.18.11
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin

###########
# START ! #
###########

echo "# Hostname (aka node) and date"
hostname
date
printf "\n"

# BAM prep -----------------------------------------------------------
tmp_sub_picard=${tmp_sub}_picard # tmp file for picard
mkdir $TMP/$tmp_sub/$tmp_sub_picard

echo "### Sorting BAM file"
time java -jar $PICARD/picard.jar SortSam \
	I=$TMP/$tmp_sub/$bam \
	O=$TMP/$tmp_sub/${bam/.bam/.sorted.bam} \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP/$tmp_sub/$tmp_sub_picard
printf "\n"

echo "### Indexing BAM file ..."
samtools index $TMP/$tmp_sub/${bam/.bam/.sorted.bam}
printf "\n"

echo "### Adding Read Group ..."
time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
	I=$TMP/$tmp_sub/${bam/.bam/.sorted.bam} \
	O=$TMP/$tmp_sub/${bam/.bam/.sorted.RG.bam} \
	RGLB=$LB \
	RGPL=$PLILLUMINA \
	RGSM=$SM \
	RGPU=$PU
printf "\n"

echo "### Marking duplicates ..."
time java -jar $PICARD/picard.jar MarkDuplicates \
	I=$TMP/$tmp_sub/${bam/.bam/.sorted.RG.bam} \
	O=$TMP/$tmp_sub/${bam/.bam/.sorted.RG.dedup.bam} \
	M=$FINAL/${bam/.bam/.sorted.RG.dedup.metrics.txt} \
	TMP_DIR=$TMP/$tmp_sub/$tmp_sub_picard
	printf "\n"


