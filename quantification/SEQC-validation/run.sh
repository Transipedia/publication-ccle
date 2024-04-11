#!/bin/sh

WK_DIR="~/Projects/working-bench/ReindeerAppli/" # General working directory

TAQMAN_TAB=$WK_DIR"data/merged-Taqman-raw_reindeer-IDs_tab.tsv" # Taqman table
ENSG2SYMBOL_PATH="geneSymbol_SEQC_to_ENST.tsv" # provided in the metadata folder

REINDEER_ONLY_TAB=$WK_DIR"data/reindeerOnly/SEQC-genes-whole-canonical-seq_on_SEQC_raw_counts_k31.out"
REINDEER_ONLY_OUT_DIR=$WK_DIR"res_cmp2SEQC/reindeeronly/"

KMERATOR_REINDEER_TAB=$WK_DIR"data/kmerator-reindeer/SEQC-genes_on_SEQC_raw_counts_k31_kmerator.out"
KMERATOR_REINDEER_OUT_DIR=$WK_DIR"res_cmp2SEQC/reindeer-kmerator2taqman/"

# Reindeer only comparison
Rscript a1_ReindeerOnly_vs_GroundTruth.R $TAQMAN_TAB $ENSG2SYMBOL_PATH \
                                         $REINDEER_ONLY_TAB $REINDEER_ONLY_OUT_DIR

# Kmerator-Reindeer comparison
Rscript a2_Kmerator-Reindeer_vs_GroundTruth.R $TAQMAN_TAB $KMERATOR_REINDEER_TAB $KMERATOR_REINDEER_OUT_DIR
