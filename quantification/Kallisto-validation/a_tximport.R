rm(list = ls())

library(tximport)
library(stringr)
library(GenomicFeatures)

# WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
# GTF_PATH <- paste0(WKDIR, "data/Homo_sapiens.GRCh38.108.gtf")
# KALLISTO_DIR <- paste0(WKDIR, "data/kallisto-ensembl108/")

cmd.Args <- commandArgs(trailingOnly = TRUE)
GTF_PATH <- cmd.Args[1]
KALLISTO_DIR <- cmd.Args[2]

samples <- list.files(KALLISTO_DIR) %>%
    str_subset(pattern = "seqc")
print(paste("Sample number: ", length(samples)))

abundance.files <- paste0(paste0(KALLISTO_DIR), samples, "/abundance.h5")
names(abundance.files) <- samples

txdb <- makeTxDbFromGFF(file = GTF_PATH)
tx2gene <- select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
    
gene.kallisto <- tximport(files = abundance.files, type = "kallisto", 
                          txIn = TRUE, txOut = FALSE, tx2gene = tx2gene, ignoreTxVersion = T)
gene.counts <- as.data.frame(gene.kallisto$counts)
write.table(gene.counts, file = paste0(KALLISTO_DIR, "gene-counts-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
gene.abund <- as.data.frame(gene.kallisto$abundance)
write.table(gene.abund, file = paste0(KALLISTO_DIR, "gene-abundance-tximport.tsv"), 
            col.names = T, row.names = T, sep = "\t", quote = F)    
