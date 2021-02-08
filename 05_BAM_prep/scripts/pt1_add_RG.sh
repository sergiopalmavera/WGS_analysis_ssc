#!/bin/bash

# This scripts sorts alignments, adds read groups and marks duplicates

# Sample id 
sample_id=$1

echo "# Adding read groups to bam files for sample $sample_id"
printf "\n"

# set bam directories
bam_dir=../../04_Alignments_Raw
TMP=../output/TMP
tmp_sub=tmp_${sample_id} #temp dir for this sample
mkdir $TMP/$tmp_sub # make a tmp dir for this sample

# absolute paths modify accordingly
PICARD=~/FBN_HOME/Tools/picard_2.18.11
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin

# Get all bams for this sample
bams=$(ls -1 $bam_dir | grep $sample_id )
echo "## BAMs for this sample"
for bam in $bams; do echo $bam; done
printf "\n"

for bam in $bams
do 
	echo "## Processing bam: $bam"
	printf "\n"

	# set read groups
	PU=$(echo $bam | cut -d'_' -f3 | sed 's/.bam//g')
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

	echo "## Adding Read Group ..."
	time java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
		I=$bam_dir/$bam \
		O=$TMP/$tmp_sub/${bam/.bam/.RG.bam} \
		RGLB=$LB \
		RGPL=$PL \
		RGSM=$SM \
		RGPU=$PU
	printf "\n"
done


