# Define samples (there are 20 samples)
SAMPLES=$(ls -1 *.NBP | cut -d'-' -f1 | sort | uniq)

# Collect all .NBP files for each sample (there are 8 files per sample)
for sample in $SAMPLES
do
	echo "# Processing sample $sample"
	files=$(ls -1 *.NBP | grep $sample)
	cat $files > TOTAL_NBP_$sample.txt
	echo "## Done processing sample $sample"
	printf "\n"
done
