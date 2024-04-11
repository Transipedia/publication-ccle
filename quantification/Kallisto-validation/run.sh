#!/bin/sh

WK_DIR="~/Projects/working-bench/ReindeerAppli/" # General working directory

GTF_PATH=$WK_DIR"data/Homo_sapiens.GRCh38.108.gtf" # GTF file downloaded from Gencode website
ENST2SYMBOL_PATH=$WKDIR"data/geneSymbol_SEQC_to_ENST.tsv" # provided in the current directory
MISSING_GENG_PATH=$WKDIR"data/missing_genes.tsv" # provided in the current directory

KALLISTO_DIR=$WK_DIR"data/kallisto-ensembl108/" # Kallisto folder with results produced by a_tximport.R
REINDEER_ONLY_DIR=$WK_DIR"data/reindeeronly/" # Reindeer query results without Kmerator
KMERATOR_REINDEER_DIR=$WK_DIR"data/kmerator-reindeer/" # Reindeer query results with Kmerator

REINDEER_ONLY_OUT_DIR=$WK_DIR"res_cmp2Kallisto/reindeeronly/"
KMERATOR_REINDEER_OUT_DIR=$WK_DIR"res_cmp2Kallisto/reindeer-kmerator2kallisto/"

# Tximport to convert transcript quantifications to gene level
Rscript a_tximport.R $GTF_PATH $KALLISTO_DIR

# Compare Reindeer (without Kmerator) to Kallisto-tximport results
Rscript b1_ReindeerOnly_vs_Kallisto.R $ENST2SYMBOL_PATH $GTF_PATH $MISSING_GENG_PATH \
                                      $KALLISTO_DIR $REINDEER_ONLY_DIR \
                                      $REINDEER_ONLY_OUT_DIR

# Compare Reindeer (with Kmerator) to Kallisto-tximport results
Rscript b2_Kmerator-Reindeer_vs_Kallisto.R $MISSING_GENG_PATH $GTF_PATH \
                                           $KALLISTO_DIR $KMERATOR_REINDEER_DIR \
                                           $KMERATOR_REINDEER_OUT_DIR
