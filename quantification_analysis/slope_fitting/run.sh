#!/bin/sh

WK_DIR="~/Projects/working-bench/ReindeerAppli/" # General working directory

CMP_TAB=$WK_DIR"res_cmp2Kallisto/reindeer-kmerator2kallisto/kmerator-reindeer.tsv"
OUT_DIR=$WK_DIR"fitting_slope/"

Rscript a_show_by_gene.R $CMP_TAB $OUT_DIR
Rscript b_summarised_fitting.R $CMP_TAB $OUT_DIR
