#!/bin/bash

# Define sample and corresponding project file
sample_id=SAMEA3376937
tsv_file=PRJEB9326.tsv
multi_fastq_pairs=$(grep $sample_id $tsv_file | awk 'BEGIN { FS = "\t" } ; {print $7}')

# Define relative paths 
REF=../../ref_ensemble95/
TMP=../output/TMP
FINAL=../output/FINAL

# absolute paths modify accordingly
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC
FASTP=~/FBN_HOME/Tools/fastp


echo "# Processing multiple FASTQs for sample $sample_id"
for i in $multi_fastq_pairs; do echo $i; done
printf "\n\n"

for fq in $multi_fastq_pairs
do
	# Separate tsv entry into links for fq1 and fq2 
	fq1_link=$(echo $fq | awk 'BEGIN { FS = ";"} ; { print $1 }')
	fq2_link=$(echo $fq | awk 'BEGIN { FS = ";"} ; { print $2 }')

	# Print a message to show which read pairs are being processed"
	echo "# Processing read pairs:"
	echo $fq1_link
	echo $fq2_link
	printf "\n"

	run_accession=$(echo $fq1_link | awk 'BEGIN { FS = "/" } ; { print $5 }')  # get accesstion name, same for both links
	base_nm=${sample_id}_${run_accession}
	echo "# Base name:" # inform the base name the files will have
	echo $base_nm
	printf "\n"

	# Define fastq and bam file names
	new_R1=${base_nm}_R1.fastq.gz
	new_R2=${base_nm}_R2.fastq.gz
	fastp_html=${base_nm}_R1_R2_fastp.html
	fastp_json=${base_nm}_R1_R2_fastp.json
	new_R1_clean=${new_R1/.fastq/.clean.fastq.gz}
	new_R2_clean=${new_R2/.fastq/.clean.fastq.gz}
	bam=${base_nm}.bam

	###########
	# START ! #
	###########

	echo "# Hostname (aka node) and date"
	hostname
	date
	printf "\n"

	# Get files from ENA ---------------------------------------------------
	echo "### Getting fastq files sample_run: $base_nm"
	tmp_sub=tmp_${base_nm} #temp dir for this sample

	wget --no-verbose $fq1_link -P $TMP/$tmp_sub
	if [ $? -eq 0 ]; then echo "download succesfull"; else echo "download failed"; fi #check if success
	printf "\n"

	wget --no-verbose $fq2_link -P $TMP/$tmp_sub
	if [ $? -eq 0 ]; then echo "download succesfull"; else echo "download failed"; fi #check if success
	printf "\n"

	# change name fastq files adding sample name to file ------------------
	echo "### changing raw fastq file names by adding sample name"
	old_R1=$(basename $fq1_link)
	mv $TMP/$tmp_sub/$old_R1 $TMP/$tmp_sub/$new_R1
	printf "\n"

	old_R2=$(basename $fq2_link)
	mv $TMP/$tmp_sub/$old_R2 $TMP/$tmp_sub/$new_R2 
	printf "\n"

	# Get FastQC reports from raw fastqs ---------------------------------
	echo "### FastQC on raw fastq files"
	time $FASTQC/fastqc -t 2 -out $FINAL $TMP/$tmp_sub/$new_R1 $TMP/$tmp_sub/$new_R2 
	printf "\n"

	# Clean with fastp ---------------------------------------------------
	echo "### Cleaning reads with fastp"
	time $FASTP/fastp -h $FINAL/$fastp_html -j $FINAL/$fastp_json -i $TMP/$tmp_sub/$new_R1 -I $TMP/$tmp_sub/$new_R2 -o $TMP/$tmp_sub/$new_R1_clean -O $TMP/$tmp_sub/$new_R2_clean
	printf "\n"

	# Get FastQC reports from clean fastqs -------------------------------
	echo "### FastQC on clean fastq files"
	time $FASTQC/fastqc -t 2 -out $FINAL $TMP/$tmp_sub/$new_R1_clean $TMP/$tmp_sub/$new_R2_clean
	printf "\n"	

	# Make alignments ----------------------------------------------------
	echo "### bwa-mem"
	bwa mem -t 20 -M $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa $TMP/$tmp_sub/$new_R1_clean $TMP/$tmp_sub/$new_R2_clean | samtools view -@ 20 -bS - > $TMP/$tmp_sub/$bam
	if [ $? -eq 0 ]; then echo "bwa-mem succesfull"; else echo "bwa-mem failed"; fi #check if success
	printf "\n"
done

