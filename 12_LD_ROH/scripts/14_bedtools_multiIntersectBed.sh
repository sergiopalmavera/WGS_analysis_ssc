#!/bin/bash

bedtools=~/FBN_HOME/Tools/bedtools_version_2.29.2/bedtools.static.binary
bed_dir=../../00_analysis/data/bed_files_ROH

fls=$(ls $bed_dir/*.bed)
fl_id=$(for f in $fls; do basename ${f/.bed/}; done)
$bedtools multiinter -i $fls -names $fl_id > ../output/SOSFERT_ROH_intersections.bed
