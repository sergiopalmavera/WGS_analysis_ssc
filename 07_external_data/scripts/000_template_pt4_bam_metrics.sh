#!/bin/bash

# This scripts collects metrics of the final BAM in the pipeline
# This BAM has been sorted, read group added, merged (if applicable, if different lanes for same samples), marked for duplicated and BQSR'd

# Sample id & FASTQ links 
sample_id=

# Manually define bam 
base_nm=
bam=${base_nm}.bam

# Define relative paths 
REF=../../ref_ensemble95/
FINAL=../output/FINAL

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

echo "### Collect Metrics after BAM preparation (Q20 MQ20)"
time java -jar $PICARD/picard.jar CollectWgsMetrics \
	I=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} \
	O=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.CollectWgsMetrics_q20mq20.txt} \
	R=$REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	Q=20 \
	MQ=20
printf "\n"

echo "### Collect Metrics after BAM preparation (Q0 MQ0)"
time java -jar $PICARD/picard.jar CollectWgsMetrics \
	I=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} \
	O=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.CollectWgsMetrics_q0mq0.txt} \
	R=$REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	Q=0 \
	MQ=0
printf "\n"

echo "### Running Picard CollectInsertSizeMetrics"
time java -jar $PICARD/picard.jar CollectInsertSizeMetrics \
	I=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} \
	O=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.CollectInsertSizeMetrics.tab} \
	H=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.insert_size_hist.pdf}
printf "\n"

echo "### starting picard tools BamIndexStats"
time java -jar $PICARD/picard.jar BamIndexStats	I=$FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} > $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.BAMIndexStats.tab}
printf "\n"

echo "### Indexing final BAM"
time $SAMTOOLS/samtools index $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam}
printf "\n"

echo "### Running flagstat on final BAM"
time $SAMTOOLS/samtools flagstat $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.bam} > $FINAL/${bam/.bam/.sorted.RG.dedup.bqsr.flagstat.tab}
printf "\n"


