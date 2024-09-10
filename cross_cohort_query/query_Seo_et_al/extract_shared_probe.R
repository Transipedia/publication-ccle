rm (list = ls())

library(readxl)
library(BSgenome.Hsapiens.UCSC.hg19)

WKDIR <- "../../../../../working-bench/ReindeerAppli/res_query_luadSEO/"

tab.probe <- read_excel(paste0(WKDIR, "SupplTabs_Transipedia.xlsx"), sheet = 1, skip = 1,
                  col_names = TRUE, col_types = "text") %>%
    dplyr::rename("Mutation" = `Mutation (HG19 coord.)`) %>%
    dplyr::filter(!str_detect(Mutation, pattern = "_-_|_-"))

tab.true <- read_excel(paste0(WKDIR, "SuppTable3-Seo.xls"), sheet = 1, skip = 1,
                       col_names = TRUE, col_types = "text") %>%
    dplyr::rename("wt" = `wt\nallele`, mut = `var\nallele`) %>%
    dplyr::mutate("Mutation" = paste0(chr, "_", pos, "_", wt, "_", mut),
                  "Probe" = getSeq(BSgenome.Hsapiens.UCSC.hg19, names = paste0("chr", chr),
                                   start = as.numeric(pos) - 30, end = as.numeric(pos) + 30, strand = "+",
                                   as.character = TRUE)) %>%
    dplyr::select(c(Probe, Mutation, chr, pos, wt, mut, annotation))
print(all(str_sub(tab.true$Probe, start = 31, end = -31) == tab.true$wt))
str_sub(tab.true$Probe, start = 31, end = -31) <- tab.true$mut
print(all(str_sub(tab.true$Probe, start = 31, end = -31) == tab.true$mut))

tab.probe %>%
    dplyr::filter(Probe %in% tab.true$Probe) %>%
    write.table(file = paste0(WKDIR, "shared_query/shared_probe.tsv"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
tab.probe %>%
    dplyr::filter(!(Probe %in% tab.true$Probe)) %>%
    write.table(file = paste0(WKDIR, "shared_query/transipedia_specific_probe.tsv"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
tab.true %>%
    dplyr::filter(!(Probe %in% tab.probe$Probe)) %>%
    write.table(file = paste0(WKDIR, "shared_query/Seo_specific_probe.tsv"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)