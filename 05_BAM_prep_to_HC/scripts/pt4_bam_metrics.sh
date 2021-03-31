#!/bin/bash

# This scripts collects metrics of the final BAM in the pipeline
# This BAM has been sorted, read group added, merged (if applicable, if different lanes for same samples), marked for duplicated and BQSR'd

sample_nm=$1

# Define relative paths 
REF=../../ref_ensemble95/

# absolute paths modify accordingly
PICARD=~/FBN_HOME/Tools/picard_2.18.11
SAMTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/samtools-1.9_installed/bin


echo "Processing BAMs"
ls -l ../output/*${sample_nm}*.bam

for bam in ../output/*${sample_nm}*.bam
do
	echo "# Processing BAM: $bam"

	echo "### Collect BAM Metrics (Q20 MQ20)"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$bam \
		O=${bam/.bam/.CollectWgsMetrics_q20mq20.txt} \
		R=$REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
		Q=20 \
		MQ=20
	printf "\n"

	echo "### Collect BAM Metrics (Q0 MQ0)"
	time java -jar $PICARD/picard.jar CollectWgsMetrics \
		I=$bam \
		O=${bam/.bam/.CollectWgsMetrics_q0mq0.txt} \
		R=$REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
		Q=0 \
		MQ=0
	printf "\n"

	echo "### Running Picard CollectInsertSizeMetrics"
	time java -jar $PICARD/picard.jar CollectInsertSizeMetrics \
		I=$bam \
		O=${bam/.bam/.CollectInsertSizeMetrics.tab} \
		H=${bam/.bam/.insert_size_hist.pdf}
	printf "\n"

	echo "### Indexing final BAM"
	time $SAMTOOLS/samtools index $bam
	printf "\n"

	echo "### Running flagstat on final BAM"
	time $SAMTOOLS/samtools flagstat $bam > ${bam/.bam/.flagstat.tab}
	printf "\n"
done

