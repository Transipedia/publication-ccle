rm(list = ls())

library(stringr)
library(tidyr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
INDIR <- paste0(WKDIR, "data/")
OUTDIR <- paste0(WKDIR, "res_cmp2SEQC/")

# Load and parse Reindeer results - with Kmerator
tab.with <- read.table(paste0(INDIR, "kmerator-reindeer/SEQC-genes_on_SEQC_raw_counts_k31_kmerator.out"),
                       header = TRUE, sep = "\t") %>%
    dplyr::mutate(seq_name = str_replace(seq_name, "C10orf93", "C10orf92")) %>%
    unique() %>%
    pivot_longer(cols = -seq_name, values_to = "Query.reindeer", names_to = "Sample.assay") %>%
    dplyr::mutate(Symbol = str_extract(seq_name, pattern = "ENST[0-9]*")) %>%
    dplyr::rename(Query.with = Query.reindeer)
tab.with <- aggregate(tab.with$Query.with,
                      by = list(tab.with$Symbol, tab.with$Sample.assay),
                      FUN = function(qs) paste0(qs, collapse = ",")) %>%
    dplyr::rename(Symbol = Group.1, Sample.assay = Group.2, Query.with = x) %>%
    dplyr::mutate(Count.reindeer = lapply(Query.with, 
                                          FUN = function(qstr) {
                                              qs.parsed <- strsplit(qstr, split = ",") %>% unlist()
                                              count.vect <- NULL
                                              for (s in qs.parsed) {
                                                  if (nchar(s) > 1) { # if the query is not a single "*"
                                                      pos.info <- strsplit(s, split = "-|:") %>% unlist()
                                                      b <- as.integer(pos.info[1])
                                                      e <- as.integer(pos.info[2])
                                                      q <- pos.info[3]
                                                      count.vect <- c(count.vect,
                                                                      rep(ifelse(q != "*",
                                                                                 yes = as.integer(q),
                                                                                 no = 0),
                                                                          e - b + 1))
                                                  }
                                              }
                                              return(count.vect)
                                          })) %>%
    dplyr::mutate(with.Kmerator = lapply(Count.reindeer,
                                         FUN = function(x) ifelse(!is.null(x), yes = mean(x), no = 0)) %>% unlist()) %>%
    dplyr::select(-Count.reindeer)
    
# Load and parse Reindeer results - no Kmerator
tab.without <- read.table(paste0(INDIR, "reindeerOnly/SEQC-genes-whole-canonical-seq_on_SEQC_raw_counts_k31.out"),
                          header = TRUE, sep = "\t") %>%
    dplyr::mutate(seq_name = str_remove_all(seq_name, " ")) %>%
    dplyr::rename(Symbol = seq_name) %>%
    pivot_longer(cols = -Symbol, values_to = "Query.reindeer", names_to = "Sample.assay") %>%
    dplyr::rename(Query.without = Query.reindeer)
tab.without <- aggregate(tab.without$Query.without,
                         by = list(tab.without$Symbol, tab.without$Sample.assay),
                         FUN = function(qs) paste0(qs, collapse = ",")) %>%
    dplyr::rename(Symbol = Group.1, Sample.assay = Group.2, Query.without = x) %>%
    dplyr::mutate(Count.reindeer = lapply(Query.without, 
                                          FUN = function(qstr) {
                                              qs.parsed <- strsplit(qstr, split = ",") %>% unlist()
                                              count.vect <- NULL
                                              for (s in qs.parsed) {
                                                  if (nchar(s) > 1) { # if the query is not a single "*"
                                                      pos.info <- strsplit(s, split = "-|:") %>% unlist()
                                                      b <- as.integer(pos.info[1])
                                                      e <- as.integer(pos.info[2])
                                                      q <- pos.info[3]
                                                      count.vect <- c(count.vect,
                                                                      rep(ifelse(q != "*",
                                                                                 yes = as.integer(q),
                                                                                 no = 0),
                                                                          e - b + 1))
                                                  }
                                              }
                                              return(count.vect)
                                          })) %>%
    dplyr::mutate(without.Kmerator = lapply(Count.reindeer,
                                            FUN = function(x) ifelse(!is.null(x), yes = mean(x), no = 0)) %>% unlist()) %>%
    dplyr::select(-Count.reindeer)

# Merged table
tab.merged <- merge(tab.with, tab.without)
write.table(tab.merged, paste0(OUTDIR, "REINDEER_query_mean.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
