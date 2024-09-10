rm (list = ls())

WKDIR <- "../../../../../working-bench/ReindeerAppli/res_query_luadSEO/"

shared.tab <- read.table(paste0(WKDIR, "shared_query/shared_probe.tsv"),
                         header = TRUE, sep = "\t")
fa <- NULL
for (i in 1 : nrow(shared.tab)) {
    fa <- c(fa, paste0(">", shared.tab$Gene[i], ":", shared.tab$Mutation[i],
                       "\n", shared.tab$Probe[i]))
}
fp <- file(paste0(WKDIR, "shared_query/shared_probe.fa"))
writeLines(fa, fp, sep = "\n")
close(fp)