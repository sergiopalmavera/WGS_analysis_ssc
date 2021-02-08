#!/bin/bash

# This script runs BQSR and HaplotypeCaller on prepared BAMs
# BAMs used as input have been sorted, RG added, merged if applicable and marked for duplicates

# Sample id & FASTQ links 
sample_id=

# Manually define bam 
base_nm=
bam=${base_nm}.bam

# Define relative paths 
REF=../../ref_ensemble95/
TMP=../output/TMP
FINAL=../output/FINAL
tmp_sub=tmp_${sample_id} #temp dir for this sample

# absolute paths modify accordingly
GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0

###########
# START ! #
###########

echo "# Hostname (aka node) and date"
hostname
date
printf "\n"

echo "### BaseRecalibrator ..."
time $GATK/gatk BaseRecalibrator \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $TMP/$tmp_sub/${bam/.bam/.sorted.RG.dedup.bam} \
	--known-sites $REF/variation/sus_scrofa_NoWhiteSpacesInINFO_sorted.vcf \
	-O $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.table} 
printf "\n"

echo "### ApplyBQSR ..."
time $GATK/gatk ApplyBQSR \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $TMP/$tmp_sub/${bam/.bam/.sorted.RG.dedup.bam} \
	--bqsr-recal-file $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.table} \
	-O $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam}
printf "\n"

echo "### HaplotypeCaller on bam $b ..."
time $GATK/gatk HaplotypeCaller \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-I $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} \
	-O $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.g.vcf} \
	-ERC GVCF
printf "\n"







