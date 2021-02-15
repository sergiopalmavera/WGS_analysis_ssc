#!/bin/bash

# This script runs BQSR and HaplotypeCaller on prepared BAMs
# BAMs used as input have been sorted, RG added, merged if applicable and marked for duplicates

# Sample id & FASTQ links 
sample_id=$1
 
# Define relative paths 
REF=../../ref_ensemble95/
TMP=../output/TMP
FINAL=../output/FINAL
tmp_sub=tmp_${sample_id} #temp dir for this sample

# define bam
tmp_sub=tmp_${sample_id} #temp dir for this sample
bam=$(echo $TMP/$tmp_sub/*.sorted.dedup.bam) # find bam

echo "# Processing sample $sample_id prepared bam file $bam"

# absolute paths modify accordingly
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

###########
# START ! #
###########

echo "# Hostname (aka node) and date"
hostname
date
printf "\n"

bam_nm=$(basename $bam) # keep only the bam name

echo "### BaseRecalibrator ..."
time $GATK/gatk BaseRecalibrator \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $bam \
	--known-sites $REF/variation/sus_scrofa_NoWhiteSpacesInINFO_sorted.vcf \
	-O $FINAL/${bam_nm/.bam/.bqsr.table} 
printf "\n"

echo "### ApplyBQSR ..."
time $GATK/gatk ApplyBQSR \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $bam \
	--bqsr-recal-file $FINAL/${bam_nm/.bam/.bqsr.table} \
	-O $FINAL/${bam_nm/.bam/.bqsr.bam}
printf "\n"

echo "### HaplotypeCaller on bam $bam_nm ..."
time $GATK/gatk HaplotypeCaller \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $FINAL/${bam_nm/.bam/.bqsr.bam} \
	-O $FINAL/${bam_nm/.bam/.bqsr.g.vcf} \
	-ERC GVCF
printf "\n"







