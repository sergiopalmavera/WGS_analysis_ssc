#!/bin/bash

# Define Chunk
FROM=41
TO=80

# Define paths
FASTQ=/projekte/I2-SOS-FERT-pig/Original
FASTQC=/home/fb4/palma-vera/FBN_HOME/Tools/FastQC

# Define samples to process
FILES=$(ls -1 $FASTQ | grep "fastq.gz" | sed -n $FROM,${TO}p )
cd $FASTQ
$FASTQC/fastqc -t 40 -out . $FILES
