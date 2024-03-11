rm(list = ls())

library("stringr")

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"

nkmer.ctg <- read.table(paste0(WKDIR,
                               "data/masked_kmer_ratios/contigs_uniq_k31_nb-kmers.tsv"),
                        header = FALSE)
colnames(nkmer.ctg) <- c("contig", "nb.kmer")
nkmer.ctg$gene <- lapply(nkmer.ctg$contig,
                         FUN = function(x) str_split(x, ":")[[1]][1]) %>%
    unlist()
nkmer.ctg$transcript <- str_extract(nkmer.ctg$contig, pattern = "ENST[0-9]+")
nkmer.ctg <- nkmer.ctg[nkmer.ctg$gene != "C10orf93", ]

nkmer <- aggregate(nkmer.ctg$nb.kmer, 
                   by = list(nkmer.ctg$gene, nkmer.ctg$transcript),
                   FUN = sum)
colnames(nkmer) <- c("Symbol", "ID", "nb.kmer.after")

nkmer.tx <- read.table(paste0(WKDIR, 
                              "data/masked_kmer_ratios/SEQC_canonical_transcript_sequences_nb-kmers.tsv"),
                       header = FALSE)
colnames(nkmer.tx) <- c("ID", "nb.kmer.before")
nkmer.tx$ID <- str_remove(nkmer.tx$ID, pattern = "\\.[0-9]*")

nkmer <- merge(nkmer, nkmer.tx)
nkmer$ratio.masked <- nkmer$nb.kmer.after / nkmer$nb.kmer.before

cmp.tab <- read.table(paste0(WKDIR, "res_cmp2Kallisto/reindeer-kmerator2kallisto/kmerator-reindeer.tsv"),
                      header = TRUE, sep = "\t")
fit.tab <- cmp.tab %>%
    dplyr::filter((Symbol != "C10orf93") & (Sum.reindeer > 0) & (Kallisto.counts > 0)) %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(mean.Sum.reindeer = mean(Sum.reindeer),
                     sd.Sum.reindeer = sd(Sum.reindeer),
                     mean.Kallisto.counts = mean(Kallisto.counts),
                     sd.Kallisto.counts = sd(Kallisto.counts),
                     slope = lm(formula = Sum.reindeer ~ Kallisto.counts)$coefficients[2],
                     R2 = summary(lm(formula = Sum.reindeer ~ Kallisto.counts))$adj.r.squared,
                     n_smp = dplyr::n())
fit.tab <- merge(fit.tab, nkmer, all = TRUE)
