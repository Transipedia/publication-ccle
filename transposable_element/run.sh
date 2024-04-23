#!/bin/sh

WK_DIR="~/Projects/working-bench/ReindeerAppli/" # General working directory

SMP_LIST="ERV_column_names.tsv" # sample names of the table, provided in the metadata folder

TELESCOPE_TAB=$WK_DIR"data/1000ERV-56samples/Telescope_rmsk_colon_RAW.tsv" # Telescope query results
ERVS2CONSIDER="TELESCOPE_ERVlist.1000.tsv" # 1000 ERVs to consider, provided in the metadata folder

REDISCOVERTE_CPM_TAB=$WK_DIR"data/1000ERV-56samples/REdiscoverTE_colon_CPM.tsv"
REDISCOVERTE_RAW_TAB=$WK_DIR"data/1000ERV-56samples/REdiscoverTE_raw.tsv"

REINDEER_TELESCOPE_RES=$WK_DIR"data/1000ERV-56samples/contigs_HERV_on_CCLE-56-cut.out"
REINDEER_REDISCOVERTE_RES=$WK_DIR"data/1000ERV-56samples/REdiscoverERVs_max_tx_100_contig_spe_on_CCLE-1019-cut-disk_56-2.out"

OUT_DIR=$WK_DIR"res_HERV1000/"

# Compare to Telescope
Rscript b1_cmp1000ERV.R $TELESCOPE_TAB $SMP_LIST $ERVS2CONSIDER \
                        $REINDEER_TELESCOPE_RES \
                        $OUT_DIR

# Compare to REdiscoverTE, both for CPM and raw counts
Rscript b2_cmp2REdiscoverTE.R $REDISCOVERTE_CPM_TAB $REDISCOVERTE_RAW_TAB $SMP_LIST \
                              $REINDEER_REDISCOVERTE_RES \
                              $OUT_DIR
