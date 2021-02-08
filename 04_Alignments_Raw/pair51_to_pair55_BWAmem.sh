FROM=51
TO=55 #max80

THR=20

REF=../ref_ensemble95/

FASTQ=../02_qc_trimming_filtering_fastp/01_single_fastq_files

R1s=$(ls -1 $FASTQ | grep "fastq.gz" | grep "R1" | sed -n $FROM,${TO}p )

for R1 in $R1s; do echo "# Processing pairs ${R1/R1/[R1,R2]}"; done

for R1 in $R1s
do
	echo "## Processing pair ${R1/R1/[R1,R2]}"
	R2=$(echo ${R1/R1/R2})
	echo "### R1: $R1"
	echo "### R2: $R2"
	FLNM=$(echo $R1 | cut -d'_' -f1,2,3 | cut -d'.' -f1)

	echo "### Starting bwa mem"
	#time bwa mem -t $THR -M $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz $FASTQ/$R1 $FASTQ/$R2 > ./$FLNM.sam
	time bwa mem -t $THR -M $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa $FASTQ/$R1 $FASTQ/$R2 | samtools view -@ $THR -bS - > ./$FLNM.bam 
	printf "\n"

	#echo "### Starting samtools view sam to bam"
	#time samtools view -@ $THR -bS ./$FLNM.sam > ./$FLNM.bam
	#if [ -e ./$FLNM.bam ]; then rm ./$FLNM.sam; fi
	#printf "\n"

	echo "### Process completed, this is no confirmation that all went well, you need to inspect the BAM files!"
done
