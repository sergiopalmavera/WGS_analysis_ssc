FROM=$1
TO=$2
FASTP=~/FBN_HOME/Tools/fastp
IN=/projekte/I2-SOS-FERT-pig/Original
OUT=.

R1s=$( ls -1 $IN | grep "fastq.gz" | grep "R1" | sed -n $FROM,${TO}p )

echo "# Processing pairs" 
for R1_fastq in $R1s; do echo ${R1_fastq/R1/[R1,R2]}; done
printf "\n"

for R1_fastq in $R1s
do
	echo "## Processing pair ${R1_fastq/R1/[R1,R2]}"
	printf "\n"
	R2_fastq=${R1_fastq/R1/R2}
	NAME_REPORT=$( echo ${R1_fastq/R1/[R1,R2]} | sed 's/.fastq.gz//g')
	#echo $NAME_REPORT
	time $FASTP/fastp -h ${NAME_REPORT}.html -j ${NAME_REPORT}.json -i $IN/$R1_fastq -I $IN/$R2_fastq -o $OUT/${R1_fastq/.fastq/.corrected.fastq} -O $OUT/${R2_fastq/.fastq/.corrected.fastq} 
	printf "\n"
	echo "### Pair ${R1_fastq/R1/[R1,R2]} completed"
	printf "\n\n"
done
