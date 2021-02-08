#FROM
#TO

FASTQ=/projekte/I2-SOS-FERT-pig/04_qc_trimming_filtering_fastp/01_single_fastq_files
OUT=.
TMP=$OUT/TMP

FLS=$(ls -1 $FASTQ | grep ".fastq.gz" | sed -n $FROM,${TO}p)
echo "# Processing files"
for fastq in $FLS; do echo $fastq; done

time for fastq in $FLS
do
	echo "## Calculating GB/read in  $fastq"
	time zcat $FASTQ/$fastq | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $OUT/${fastq%.fastq.gz}.NBP
	time zcat $FASTQ/$fastq | grep "@" | wc -l > $OUT/${fastq%.fastq.gz}.NR
done 
