#!/bin/bash

vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

ls -1 ../../10_FinalVCF/output/*all_sites*.vcf | for in_fl in $(cat)
do
	out_fl=$(basename $in_fl | sed 's/.vcf//')
	$vcftools/vcftools --freq --vcf $in_fl --out ../output/$out_fl
done
