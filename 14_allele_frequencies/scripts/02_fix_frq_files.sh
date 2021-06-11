#!/bin/bash

for fl in ../output/*.frq
do
	# Prepare header
	echo -e "CHROM\tPOS\tN_ALLELES\tN_CHR\tREF\tREF_frq\tALT\tALT_frq" > tmp1
	# Prepare main columns
	sed -n '2,$p' $fl | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > tmp2
	# Prepare allele freq cols (ref)
	sed -n '2,$p' $fl | awk '{print $5}' | awk -F ":" '{print $1 "\t" $2}' > tmp3
	# Prepare allele freq cols (alt)
	sed -n '2,$p' $fl | awk '{print $6}' | awk -F ":" '{print $1 "\t" $2}' > tmp4
	# combine columns
	paste --delimiters='\t' tmp2 tmp3 tmp4 > tmp5
	# add header
	cat tmp1 tmp5 > ../output/${fl/.frq/.frq2}
	# remove tmp files
	rm tmp1 tmp2 tmp3 tmp4 tmp5
done
