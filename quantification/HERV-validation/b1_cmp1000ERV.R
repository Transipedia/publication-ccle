rm(list = ls())

library(tidyr)
library(stringr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
INDIR <- paste0(WKDIR, "data/1000ERV-56samples/")
OUTDIR <- paste0(WKDIR, "res_HERV1000/")

# Load and parse Taqman table
old.tab <- read.table(paste0(INDIR, "Telescope_rmsk_colon_CPM_1000.tsv"), header = TRUE)
tab.TE <- read.table(paste0(INDIR, "Telescope_rmsk_colon_RAW.tsv"),
                     header = TRUE, sep = "\t")
colnames(tab.TE) <- colnames(old.tab)
tab.TE <- tab.TE[tab.TE$Transcript %in% old.tab$Transcript, ] %>%    
    pivot_longer(cols = -Transcript, values_to = "Counts.Telescope", names_to = "Sample.assay") %>%
    dplyr::rename(Symbol = Transcript)

# Load and parse Reindeer table
tab.reindeer <- read.table(paste0(INDIR, "contigs_HERV_on_CCLE-56-cut.out"),
                           header = TRUE, sep = "\t") %>%
    pivot_longer(cols = -seq_name, values_to = "Query.reindeer", names_to = "Sample.assay") %>%
    dplyr::mutate(Symbol = str_remove(seq_name, "\\.contig[0-9]+")) %>%
    dplyr::select(-seq_name)

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
cmp.tab <- merge(tab.TE, tab.reindeer, by = c("Symbol", "Sample.assay"))
write.table(cmp.tab, file = paste0(OUTDIR, "TE_vs_Kmerator-Reindeer.4metrics.csv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cmp.tab[cmp.tab$Counts.Telescope == 0 & cmp.tab$Max.reindeer > 0, ],
            file = paste0(OUTDIR, "TE_vs_Kmerator-Reindeer.questionable.4metrics.csv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
# Compute correlations
cor.tab <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                  FUN = function(mthd) {
                      data.frame(row.names = mthd,
                                 corr.pearson = cor(cmp.tab[, "Counts.Telescope"], cmp.tab[, mthd], method = "pearson"),
                                 corr.spearman = cor(cmp.tab[, "Counts.Telescope"], cmp.tab[, mthd], method = "spearman"))
                  }) %>%
    do.call(what = rbind)
# Plot scatter density
plt_lst <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                  FUN = function(mthd) {
                      ggplot() +
                          geom_pointdensity(aes(x = cmp.tab[, "Counts.Telescope"] + 1,
                                                y = cmp.tab[, mthd] + 1),
                                            size = 0.1) +
                          geom_text(aes(x = Inf, y = 1, 
                                        label = paste0("Pearson: ", round(cor.tab[mthd, "corr.pearson"], 2), "\n",
                                                       "Spearman: ", round(cor.tab[mthd, "corr.spearman"], 2))),
                                    size = 5, hjust = 1, vjust = 0) +
                          xlab("Telescope Counts") +
                          ylab(paste0(mthd)) +
                          scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
                          scale_x_log10(breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001),
                                        labels = c(0, 10, 100, 1000, 10000, 100000, 1000000)) +
                          scale_y_log10(breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001),
                                        labels = c(0, 10, 100, 1000, 10000, 100000, 1000000)) +
                          theme_bw() +
                          theme(text = element_text(size = 15), legend.position = "none")
                  })

ggsave((plt_lst[[1]] | plt_lst[[2]]) / (plt_lst[[3]] | plt_lst[[4]]),
       filename = paste0(OUTDIR, "TE_vs_Kmerator-Reindeer.4metrics.pdf"),
       width = 9, height = 5, dpi = 600)
