#!/bin/bash

chr=$1 # go one chr at a time

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../ref_ensemble95

# collect sosfert samples gvcfs 
sosfert_gvcf=$(ls -1 ../../05_BAM_prep_to_HC/output/*.g.vcf | for i in $(cat); do echo "--variant $i"; done)

# collect external samples gvcfs
external_gvcf=$(ls -1 ../../06_external_data/output/FINAL/*.g.vcf | for i in $(cat); do echo --variant $i; done)

# Combine all gvcfs
all_gvcfs=$(echo $sosfert_gvcf $external_gvcf)

# Run tool
$GATK/gatk --java-options "-Xmx200g" CombineGVCFs \
	-R $REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	$(echo $all_gvcfs) \
	-O ../output/cohort_${chr}.g.vcf.gz \
	-L $chr
