#!/bin/bash

bedtools=~/FBN_HOME/Tools/bedtools_version_2.29.2/bedtools.static.binary
bed_dir=../../00_analysis/data/bed_files_ROH

fls=$(ls $bed_dir/*.bed)
fl_id=$(for f in $fls; do basename ${f/.bed/}; done)

res=$(for i in $fl_id; do echo "\"$i\","; done)

echo $res
