#!/bin/bash

# transform the depmap fusion annotation file into a fusion kmer file

#================================================================
# DESCRIPTION
#================================================================

SCRIPT_NAME=$( basename "$0" )
VERSION="1.0.0"

__usage="
 SYNOPSIS
    ${SCRIPT_NAME} [-hv] [-a annotfile] [-c cclefile]

 DESCRIPTION
    This is a script to create bed files from
    the depmap fusion file and generate associated 
    count tables. 

 OPTIONS
    -a [file]		   	Give SRR annotation file
    -c [file]   		Give CCLE fusion input file
    -t [INT]            Threshold (default is > 0)
    -h, --help      	Print this help
    -v, --version       Print script information 

 EXAMPLES
    ${SCRIPT_NAME} -a annot.csv -c CCLE_fusion.csv 
"

#================================================================
# ARGS
#================================================================

OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Returns the count of arguments that are in short or long options

if [ "$#" -eq 0 ]; then
  echo "$__usage"
  exit 0
fi

threshold=0

while getopts "h?va:c:t:" opt; do
    case "$opt" in
    h|\?)
        echo "$__usage"
        exit 0
        ;;
    v)  vers=$VERSION
		echo "$SCRIPT_NAME : version $VERSION"
		exit 0
        ;;
    a)  annot=${OPTARG}
        ;;
    c)  ccle=${OPTARG}
        ;;
    t)  threshold=${OPTARG}
        ;;
    esac
done


shift "$(( OPTIND - 1 ))"

if [ -z "$annot" ] || [ -z "$ccle" ]; then
        echo 'Missing -a or -c' >&2
        exit 1
fi


#================================================================
# MAIN
#================================================================

source config.sh

# All our variables 

echo -e "\nCCLE fusion file : "$ccle
echo -e "\nAnnotation file : "$annot
echo -e "\nReference genome file : "$genome
echo -e "\nReference exons file : "$exons
echo -e "\nThreshold : "$threshold

base=$(basename "$ccle");
out=${base%.*}


# STEP 1 : Convert the sample names into SRR
echo -e "\n\n STEP 1 : DepMap sample IDs to SRR"

awk -v threshold=$threshold 'BEGIN{FS=","} $3>threshold' $ccle | sed 's/\",\"/_/g' | cut -d"," -f1,2,3,6,7,8,9,17 | sort -u | grep -v FusionName | sed 's/,/\t/g' | sed 's/ /_/g' | awk 'BEGIN {FS=OFS="\t"}{print $1,$2"_"$5"_"$7,$5"_"$7,$3,$8}' | sort -u > $out"_uniq_pairs_sample-chimeras_with_counts_occ_"$threshold".tsv"
join -13 -21 <(sort -k3,3 $annot) <(sort -k1,1 $out"_uniq_pairs_sample-chimeras_with_counts_occ_"$threshold".tsv") | sed 's/\s/\t/g' > $out"_uniq_pairs_sample-chimeras_with_counts_occ_"$threshold"_SRR.tsv"
rm $out"_uniq_pairs_sample-chimeras_with_counts_occ_"$threshold".tsv"


# STEP 2 : keep the names
echo -e "\n STEP 2 : DepMap keeping ID / sample pairs"

cut -f4 $out"_uniq_pairs_sample-chimeras_with_counts_occ_"$threshold"_SRR.tsv"| sort -u | sed 's/_/\t/g' > $out"_filter_uniq-names-and-coord_"$threshold"_SRR.tsv"


# STEP 3 : generate the bed files for the left and right parts of the chimera
echo -e "\n STEP 3 : bedfiles and k-mers fasta file"

cut -f2,3 $out"_filter_uniq-names-and-coord_"$threshold"_SRR.tsv" | sort -u > tmp

get_fusions_seq () {
    left=$1
    right=$2
    length=$((left+right))
    cut -f1 tmp | awk -v left=$left 'BEGIN {FS=":";OFS="\t"}{if ($3=="+") print $1,$2-left,$2,$1":"$2":"$3,"0",$3; else print $1,$2-1,$2+left-1,$1":"$2":"$3,"0",$3}' > $out"_left-part-"$left"-bases_"$threshold".bed"
	cut -f2 tmp | awk  -v right=$right 'BEGIN {FS=":";OFS="\t"}{if ($3=="+") print $1,$2-1,$2+right-1,$1":"$2":"$3,"0",$3; else print $1,$2-right,$2,$1":"$2":"$3,"0",$3}' > $out"_right-part-"$right"-bases_"$threshold".bed"
    rm tmp
    paste -d "" <(bedtools getfasta -fi $genome -bed $out"_left-part-"$left"-bases_"$threshold".bed" -name -s) <(bedtools getfasta -fi $genome -bed $out"_right-part-"$right"-bases_"$threshold".bed" -name -s) | sed 's/\(>.*\)::.*>\(.*\)::.*/\1_\2/g' > $out"_filter_uniq-names-and-coord_"$length"_"$threshold".fa"
}

get_fusions_seq 26 25


# STEP 4 : keeping the chimeras that intersect with exon borders
echo -e "\n STEP 4 : select exonic fusions"
grep -A1 -f <(bedtools intersect -a $out"_left-part-"$left"-bases_"$threshold".bed" -b $exons -wa -s -u -f 1 | cut -f4 | awk '$0=$0"_"') $out"_filter_uniq-names-and-coord_"$length"_"$threshold".fa" --no-group-separator > $out"_filter_left-part-in-exon_"$length"_"$threshold".fa"
grep -A1 -f <(bedtools intersect -a $out"_right-part-"$right"-bases_"$threshold".bed" -b $exons -wa -s -u -f 1 | cut -f4 | sed 's/.*\(_.*\)/\1/g'  | awk '$0="_"$0') $out"_filter_left-part-in-exon_"$length"_"$threshold".fa" --no-group-separator > $out"_filter_left-and-right-part-in-exon_"$length"_"$threshold".fa"
rm $out"_filter_left-part-in-exon_"$length"_"$threshold".fa"

grep -A1 -f <(bedtools intersect -a <(awk 'BEGIN {FS=OFS="\t"}{print $7,$8,$9,$10,$11,$12,$4}' <(bedtools intersect -a <(paste $out"_left-part-"$left"-bases_"$threshold".bed" $out"_right-part-"$right"-bases_"$threshold".bed") -b $prime3 -s -wa -u)) -b $prime5 -s -wa -u | awk 'BEGIN {FS=OFS="\t"}{print $7,$4}' | sed 's/\t/_/g' | cut -d'_' -f1,2,4 ) $out"_filter_left-and-right-part-in-exon_"$length"_"$threshold".fa" --no-group-separator | head

grep -A1 -f <(bedtools intersect -a <(awk 'BEGIN {FS=OFS="\t"}{print $7,$8,$9,$10,$11,$12,$4}' <(bedtools intersect -a <(paste $out"_left-part-"$left"-bases_"$threshold".bed" $out"_right-part-"$right"-bases_"$threshold".bed") -b $prime3 -s -wa -u)) -b $prime5 -s -wa -u | awk 'BEGIN {FS=OFS="\t"}{print $7,$4}' | sed 's/\t/_/g' | cut -d'_' -f1,2,4 ) $out"_filter_left-and-right-part-in-exon_"$length"_"$threshold".fa" --no-group-separator > $out"_filter_edges_"$length"_"$threshold".fa"


# STEP 5 : kmerator to keep only specific kmers
echo -e "\n STEP 5 : kmerator on exonic fusions"
if [ -e "kmerator_"$threshold"_"$length"_edges/kmers.fa" ]
then
    echo -e "\t\tKmerator already done."
else
    kmerator -f $out"_filter_edges_"$length"_"$threshold".fa" --chimera -r 108 -o "kmerator_"$threshold"_"$length"_edges" -D
fi
sed 's/\.kmer[0-9]*//g' "kmerator_"$threshold"_"$length"_edges/kmers.fa" > $out"_filter_edges_"$length"_"$threshold"_kmerator.fa"
grep -f <(grep -v ">" $out"_filter_edges_"$length"_"$threshold"_kmerator.fa" | sort | uniq -c | awk '{print $2"\t"$1}' | awk '$2==1' | cut -f1) -B1 $out"_filter_edges_"$length"_"$threshold"_kmerator.fa" --no-group-separator > $out"_filter_edges_"$length"_"$threshold"_kmerator_spe.fa"
rm $out"_filter_edges_"$length"_"$threshold"_kmerator.fa"

# STEP 5 bis : remove low complexity kmers
echo -e "\n STEP 5' : remove low complexity kmers"
Rscript ../bin/complexity.R $out"_filter_edges_"$length"_"$threshold"_kmerator_spe.fa" $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_complex.tsv"  $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_complex.fa"
if [ ! $(awk 'BEGIN{RS=">"}{print $1"\t"$2;}' $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_complex.fa" | sed '1d' |  grep -v '\-.*AAAAAA\|\-.*TTTTTT\|\-.*GGGGGG\|\-.*CCCCCC' | awk -F'\t' -v OFS='\n' '{$1 = ">" $1} 1' | wc -l ) -eq 0 ] ; then awk 'BEGIN{RS=">"}{print $1"\t"$2;}' $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_complex.fa" | sed '1d' | grep -v '\-.*AAAAAA\|\-.*TTTTTT\|\-.*GGGGGG\|\-.*CCCCCC' | awk -F'\t' -v OFS='\n' '{$1 = ">" $1} 1' > $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_comp.fa" ; else touch $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_comp.fa" ; fi ;

# STEP 6 : 
echo -e "\n STEP 6 : get reindeer counts for all fusion k-mers"

rdeer-client query -s $server $index -q $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_comp.fa" -o $out"_filter_edges_"$length"_"$threshold"_kmerator_spe_comp_on_"$index".tsv"
