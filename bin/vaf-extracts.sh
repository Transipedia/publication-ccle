#!/bin/bash

# Author: Benoit Guibert <benoit.guibert@free.fr>
# Version 0.1.0

# DESCRIPTION
# vaf-extract.sh

input_dir=$1
output_dir=$2


### Handle headers
fp_head="gene\tmutation\tsample\trdeer_AF\trdeer_VAF"
tp_head="gene\tmutation\tsample\trdeer_AF\trdeer_VAF\tdepmap_AF\tdepmap_VAF\trdeer_alt_AF/depmap_alt_AF\trdeer_VAF/depmap_VAF"

### write headers
echo -e "${fp_head}" > ${output_dir}/false-pos.tsv
echo -e "${tp_head}" > ${output_dir}/true-pos.tsv

### grep
for file in ${input_dir}/*
do
   grep 'False+' $file | cut -f1-5 >> ${output_dir}/false-pos.tsv
   grep 'True+'  $file | cut -f1-9  >> ${output_dir}/true-pos.tsv
done
