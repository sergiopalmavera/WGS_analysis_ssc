#!/bin/bash

awk '{ if ($2 == "Large") print $1 }' internal_external_sample_info.tsv > large_white_samples
awk '{ if ($2 == "Landrace") print $1 }' internal_external_sample_info.tsv > landrace_samples

awk '{ if ($2 == "Pietrain") print $1 }' internal_external_sample_info.tsv > pietrain_samples
awk '{ if ($2 == "Duroc") print $1 }' internal_external_sample_info.tsv > duroc_samples

awk '{ if ($2 == "European") print $1 }' internal_external_sample_info.tsv > euro_wild_boar_samples
