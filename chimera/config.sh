#!/bin/bash

# Initialize our variables :

## annotation files
genome="/scratch/projects/chloe/CCLE/hg38.fa" #reference genome fasta file
exons="/scratch/projects/chloe/CCLE/gencode.v42.annotation_exons.bed" #bed file of the exons from gencode
prime3="./gencode/gencode.v42.annotation_exons-3prime.bed" #bedfile of the 3' end of exons generated from the exons bed file
prime5="./gencode/gencode.v42.annotation_exons-5prime.bed" #bedfile of the 5' end of exons generated from the exons bed file

## rdeer configuration
server="ella"
index="CCLE-1019-cut-disk"

## other parameters
min_pos_count=3
