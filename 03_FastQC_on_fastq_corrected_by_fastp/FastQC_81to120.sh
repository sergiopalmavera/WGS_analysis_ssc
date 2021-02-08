# Define Chunk
FROM=81
TO=120

CORRECTED=../02_qc_trimming_filtering_fastp/01_single_fastq_files/

FILES=$(ls -1 $CORRECTED | grep "fastq.gz" | grep "corrected.fastq.gz" | sed -n $FROM,${TO}p | for i in $(cat); do echo $CORRECTED/$i; done)
echo "# Processing files"
for i in $FILES; do echo $i; done

FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

OUT=./

$FASTQC/fastqc -t 40 -out $OUT $FILES
