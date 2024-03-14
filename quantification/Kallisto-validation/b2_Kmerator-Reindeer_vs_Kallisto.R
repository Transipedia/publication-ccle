rm(list = ls())

library(tidyr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
INDIR <- paste0(WKDIR, "data/")
OUTDIR <- paste0(WKDIR, "res_cmp2Kallisto/reindeer-kmerator2kallisto/")

if (!dir.exists(OUTDIR)) {
    dir.create(OUTDIR, recursive = TRUE)
    cat("Created output folder:", OUTDIR, "\n")
}

# Transcript ID to gene symbol
syno2symbol <- read.table(paste0(INDIR, "missing_genes.tsv"), header = FALSE, row.names = 1)

# Gene ID to gene symbol
id2symbol <- read.table(paste0(INDIR, "Homo_sapiens.GRCh38.108.gtf"), sep = "\t") %>%
    dplyr::filter(V3 == "gene")
id2symbol <- lapply(id2symbol$V9,
                    FUN = function(x) data.frame("SYMBOL" = strsplit(x, ";")[[1]][3],
                                                 "ID" = strsplit(x, ";")[[1]][1])) %>%
    do.call(what = rbind) %>%
    dplyr::mutate(SYMBOL = str_remove_all(SYMBOL, "gene_name| "),
                  ID = str_remove_all(ID, "gene_id| ")) %>%
    tibble::column_to_rownames("ID")

# Load and parse Taqman table
tab.count <- read.table(paste0(INDIR, "kallisto-ensembl108/gene-counts-tximport.tsv"),
                        header = TRUE, sep = "\t") %>%
    tibble::rownames_to_column("ID") %>%
    pivot_longer(cols = -ID, values_to = "Kallisto.counts", names_to = "Sample.assay") %>%
    dplyr::mutate(Symbol = id2symbol[ID, "SYMBOL"])
tab.abundance <- read.table(paste0(INDIR, "kallisto-ensembl108/gene-abundance-tximport.tsv"),
                            header = TRUE, sep = "\t") %>%
    tibble::rownames_to_column("ID") %>%
    pivot_longer(cols = -ID, values_to = "Kallisto.TPM", names_to = "Sample.assay") %>%
    dplyr::mutate(Symbol = id2symbol[ID, "SYMBOL"])

# Load and parse Reindeer table
tab.reindeer <- read.table(paste0(INDIR, "kmerator-reindeer/SEQC-genes_on_SEQC_raw_counts_k31_kmerator.out"),
                           header = TRUE, sep = "\t") %>%
    dplyr::mutate(seq_name = str_replace(seq_name, "C10orf93", "C10orf92")) %>%
    unique() %>%
    pivot_longer(cols = -seq_name, values_to = "Query.reindeer", names_to = "Sample.assay") %>%
    dplyr::mutate(Symbol = lapply(seq_name,
                                  FUN = function(x) strsplit(x, split = ":")[[1]][1]) %>% unlist())

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

to_consider <- tab.reindeer$Symbol %in% rownames(syno2symbol)
tab.reindeer[to_consider, "Symbol"] <- syno2symbol[tab.reindeer$Symbol[to_consider], "V2"]

# Merge tables to compare
cmp.tab <- merge(tab.count, tab.abundance, by = c("ID", "Symbol", "Sample.assay"))
cmp.tab <- merge(cmp.tab, tab.reindeer, by = c("Symbol", "Sample.assay"))

# Scatter plot and correlations
for (col2cor in c("Kallisto.counts", "Kallisto.TPM")) {
    ## Compute correlations
    cor.tab <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                      FUN = function(mthd) {
                          data.frame(row.names = mthd,
                                     corr.pearson = cor(cmp.tab[, col2cor],
                                                        cmp.tab[, mthd],
                                                        method = "pearson"),
                                     corr.spearman = cor(cmp.tab[, col2cor],
                                                         cmp.tab[, mthd],
                                                         method = "spearman"))
                      }) %>%
        do.call(what = rbind)
    ## Plot scatter density
    plt_lst <- lapply(paste0(c("Mean", "Median", "Max", "Sum"), ".reindeer"),
                      FUN = function(mthd) {
                          ggplot() +
                              geom_pointdensity(aes(x = cmp.tab[, col2cor] + 1,
                                                    y = cmp.tab[, mthd] + 1),
                                                size = 0.1) +
                              geom_text(aes(x = Inf, y = 1, 
                                            label = paste0("Pearson: ", round(cor.tab[mthd, "corr.pearson"], 2), "\n",
                                                           "Spearman: ", round(cor.tab[mthd, "corr.spearman"], 2))),
                                        size = 5, hjust = 1, vjust = 0) +
                              xlab(col2cor) +
                              ylab(mthd) +
                              scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
                              scale_x_log10(breaks = c(1, 11, 101, 1001, 10001),
                                            labels = c(0, 10, 100, 1000, 10000)) +
                              scale_y_log10(breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001, 10000001),
                                            labels = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
                              theme_bw() +
                              theme(text = element_text(size = 15), legend.position = "none")
                      })
    ggsave((plt_lst[[1]] | plt_lst[[2]]) / (plt_lst[[3]] | plt_lst[[4]]),
           filename = paste0(OUTDIR, "SupplFig_kmerator-reindeer.", col2cor, ".k31.4metrics.png"),
           width = 9, height = 5, dpi = 600)
}

# Dropout by Reindeer
all((cmp.tab$Kallisto.counts > 0) == (cmp.tab$Kallisto.TPM > 0))
all((cmp.tab$Mean.reindeer == 0) == (cmp.tab$Max.reindeer == 0))
write.table(cmp.tab, paste0(OUTDIR, "kmerator-reindeer.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cmp.tab[cmp.tab$Mean.reindeer == 0 & cmp.tab$Kallisto.counts > 0, ],
            paste0(OUTDIR, "kmerator-reindeer_Dropped.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
