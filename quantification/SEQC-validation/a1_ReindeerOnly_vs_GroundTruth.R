rm(list = ls())

library(tidyr)
library(stringr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
INDIR <- paste0(WKDIR, "data/")
OUTDIR <- paste0(WKDIR, "res_cmp2SEQC/reindeeronly/")

if (!dir.exists(OUTDIR)) {
    dir.create(OUTDIR, recursive = TRUE)
    cat("Created output folder:", OUTDIR, "\n")
}

# Load and parse Taqman table
tab.taqman <- read.table(paste0(INDIR, "merged-Taqman-raw_reindeer-IDs_tab.tsv"),
                         header = TRUE, sep = "\t") %>%
    pivot_longer(cols = -Symbol, values_to = "Counts.true", names_to = "Sample.assay")

# Load and parse Reindeer table
ensg2symbol <- read.table(paste0(INDIR, "geneSymbol_SEQC_to_ENST.tsv"), sep = "\t",
                          header = FALSE) %>%
    dplyr::filter(V1 != "C10orf93") %>%
    tibble::column_to_rownames("V2")
tab.reindeer <- read.table(paste0(INDIR, "reindeerOnly/SEQC-genes-whole-canonical-seq_on_SEQC_raw_counts_k31.out"),
                           header = TRUE, sep = "\t") %>%
    dplyr::mutate(seq_name = str_remove_all(seq_name, " ")) %>%
    dplyr::mutate(Symbol = ensg2symbol[seq_name, "V1"]) %>%
    dplyr::select(-seq_name) %>%
    pivot_longer(cols = -Symbol, values_to = "Query.reindeer", names_to = "Sample.assay") 
tab.reindeer <- aggregate(tab.reindeer$Query.reindeer,
                          by = list(tab.reindeer$Symbol, tab.reindeer$Sample.assay),
                          FUN = function(qs) paste0(qs, collapse = ",")) %>%
    dplyr::rename(Symbol = Group.1, Sample.assay = Group.2, Query.reindeer = x) %>%
    dplyr::mutate(Count.reindeer = lapply(Query.reindeer, 
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
    dplyr::mutate(Mean.reindeer = lapply(Count.reindeer,
                                         FUN = function(x) ifelse(!is.null(x), yes = mean(x), no = 0)) %>% unlist(),
                  Median.reindeer = lapply(Count.reindeer,
                                           FUN = function(x) ifelse(!is.null(x), yes = median(x), no = 0)) %>% unlist(),
                  Max.reindeer = lapply(Count.reindeer,
                                        FUN = function(x) ifelse(!is.null(x), yes = max(x), no = 0)) %>% unlist(),
                  Sum.reindeer = lapply(Count.reindeer,
                                        FUN = function(x) ifelse(!is.null(x), yes = sum(x), no = 0)) %>% unlist()) %>%
    dplyr::select(-Count.reindeer)

# Merge tables to compare
cmp.tab <- merge(tab.taqman, tab.reindeer, by = c("Symbol", "Sample.assay"))
# Compute correlations
cor.tab <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                  FUN = function(mthd) {
                      data.frame(row.names = mthd,
                                 corr.pearson = cor(cmp.tab[, "Counts.true"], cmp.tab[, mthd], method = "pearson"),
                                 corr.spearman = cor(cmp.tab[, "Counts.true"], cmp.tab[, mthd], method = "spearman"))
                  }) %>%
    do.call(what = rbind)
# Plot scatter density
plt_lst <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                  FUN = function(mthd) {
                      ggplot() +
                          geom_pointdensity(aes(x = cmp.tab[, "Counts.true"] + 1,
                                                y = cmp.tab[, mthd] + 1),
                                            size = 0.1) +
                          geom_text(aes(x = Inf, y = 1, 
                                        label = paste0("Pearson: ", round(cor.tab[mthd, "corr.pearson"], 2), "\n",
                                                       "Spearman: ", round(cor.tab[mthd, "corr.spearman"], 2))),
                                    size = 5, hjust = 1, vjust = 0) +
                          xlab("SEQC/MAQC-III qPCR") +
                          ylab(paste0(mthd)) +
                          scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
                          scale_x_log10(breaks = c(1, 6, 11, 21, 31),
                                        labels = c(0, 5, 10, 20, 30)) +
                          scale_y_log10(breaks = c(1, 6, 11, 21, 31),
                                        labels = c(0, 5, 10, 20, 30)) +
                          scale_x_log10() +
                          scale_y_log10() +
                          theme_bw() +
                          theme(text = element_text(size = 15), legend.position = "none")
                  })
plt_lst[[1]]

ggsave((plt_lst[[1]] | plt_lst[[2]]) / (plt_lst[[3]] | plt_lst[[4]]),
       filename = paste0(OUTDIR, "SupplFig_reindeerOnly.k31.4metrics.png"),
       width = 9, height = 5, dpi = 600)
