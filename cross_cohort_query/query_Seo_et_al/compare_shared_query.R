rm (list = ls())

library(readxl)
library(ggplot2)
library(stringr)
library(tidyr)
library(patchwork)
library(caret)

WKDIR <- "../../../../../working-bench/ReindeerAppli/res_query_luadSEO/"
FILT.THRES <- 3

# Query results from shared probe
qry.res <- read.table(paste0(WKDIR, "shared_query/query_results.tsv"),
                      header = TRUE, sep = "\t", row.names = 1)

# LUAD ground truth
true.res <- read_excel(paste0(WKDIR, "SuppTable3-Seo.xls"), sheet = 1, skip = 1,
                       col_names = TRUE, col_types = "text") %>%
    dplyr::mutate(seq_name = paste0(chr, "_", pos, "_", `wt\nallele`, "_", `var\nallele`)) %>%
    dplyr::filter(seq_name %in% str_split_i(rownames(qry.res), ":", 2)) %>%
    dplyr::select(c(-chr, -pos, -`wt\nallele`, -`var\nallele`, -annotation,
                    -`wt\nAminoAcid`, -`variant\nAminoAcid`, -Blosum,
                    -`is_in_COSMIC(v56)\n(LungCA)\n(ND=notDetected)`,
                    -`is_in_COSMIC(v56)\n(OtherCA)\n(ND=notDetected)`,
                    -`is_in\ndbSNP132common`, -`is_on\nsegmental\nduplication`,
                    -`#\nsamples\ninvolved`)) %>%
    pivot_longer(cols = -seq_name, names_to = "sample_name", values_to = "has.mut") %>%
    dplyr::mutate(has.mut = ifelse(has.mut == "-", yes = "FALSE", no = "TRUE"))
true.res[is.na(true.res)] <- "FALSE" # One case in Seo et al. table appear to be blank, and likely to be a '-'

# Seo et al. metadata
meta.data <- read.table(paste0(WKDIR, "filereport_read_run_PRJEB2784_tsv-2.txt"),
                        header = TRUE, sep = "\t", fill = TRUE) %>%
    dplyr::filter(run_accession %in% colnames(qry.res)) %>%
    dplyr::mutate(sample_name = str_extract(sample_title, "LC_.*"),
                  condition = str_extract(sample_title, "adenocarcinoma|adjacent_normal")) %>%
    dplyr::select(c(run_accession, sample_name, condition)) %>%
    dplyr::arrange(condition)

qry.res <- (qry.res >= FILT.THRES) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("seq_name") %>%
    pivot_longer(cols = -seq_name, names_to = "run_accession", values_to = "find.mut") %>%
    dplyr::mutate(find.mut = ifelse(find.mut, yes = "yes", no = "no"),
                  annot = seq_name,
                  seq_name = str_split_i(seq_name, ":", 2)) %>%
    merge(meta.data, by = "run_accession") %>%
    dplyr::filter(condition == "adenocarcinoma")

cmp.res <- merge(qry.res, true.res, by = c("seq_name", "sample_name")) %>%
    dplyr::mutate(`Transipedia | Ground Truth` = paste0(find.mut, " | ", has.mut))
cmp.res2plt <- pivot_wider(cmp.res, id_cols = annot,
                           names_from = run_accession,
                           values_from = `Transipedia | Ground Truth`) %>%
    dplyr::select_if(function(x) !all(x == "no | FALSE")) %>%
    pivot_longer(cols = -annot, names_to = "run_accession",
                 values_to = "Transipedia | Ground Truth")

plt <- ggplot(data = cmp.res2plt) +
    geom_tile(aes(y = annot, x = run_accession, fill = `Transipedia | Ground Truth`),
              lwd = 0.5, color = "gray") +
    scale_fill_manual(values = c("yes | TRUE" = "royalblue", "yes | FALSE" = "#e66101",
                                 "no | TRUE" = "#fdb863", "no | FALSE" = "white")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "top", text = element_text(size = 15))
ggsave(plt, filename = paste0(WKDIR, "shared_query/probe_hit_by_prob_thres", FILT.THRES, ".pdf"),
       width = 15, height = 7, dpi = 300)

print(confusionMatrix(factor(str_split_i(cmp.res2plt$`Transipedia | Ground Truth`, " \\| ", 1) == "yes"),
                      factor(str_split_i(cmp.res2plt$`Transipedia | Ground Truth`, " \\| ", 2) == "TRUE"),
                      positive = "TRUE", mode = "everything"))
