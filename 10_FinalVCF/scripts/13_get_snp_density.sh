#!/bin/bash

vcftools=~/FBN_HOME/Tools/vcftools_0.1.13/cpp/


for vcf in ../output/*filtered*.vcf 
do
	echo "Processing vcf: $vcf"
	time $vcftools/vcftools --vcf $vcf --SNPdensity 1000 --out ${vcf/.vcf/}
	printf "\n\n"
done
