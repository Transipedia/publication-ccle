rm(list = ls())

library(tximport)
library(stringr)
library(GenomicFeatures)


WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
GTFPATH <- paste0(WKDIR, "data/Homo_sapiens.GRCh38.108.gtf")
KALLISTODIR <- paste0(WKDIR, "data/kallisto-ensembl108/")

samples <- list.files(KALLISTODIR) %>%
    str_subset(pattern = "seqc")
print(paste("Sample number: ", length(samples)))

abundance.files <- paste0(paste0(KALLISTODIR), samples, "/abundance.h5")
names(abundance.files) <- samples

txdb <- makeTxDbFromGFF(file = GTFPATH)
tx2gene <- select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
    
gene.kallisto <- tximport(files = abundance.files, type = "kallisto", 
                          txIn = TRUE, txOut = FALSE, tx2gene = tx2gene, ignoreTxVersion = T)
gene.counts <- as.data.frame(gene.kallisto$counts)
write.table(gene.counts, file = paste0(KALLISTODIR, "gene-counts-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
gene.abund <- as.data.frame(gene.kallisto$abundance)
write.table(gene.abund, file = paste0(KALLISTODIR, "gene-abundance-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
