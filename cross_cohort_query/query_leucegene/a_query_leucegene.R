rm (list = ls())

library(stringr)
library(tidyr)
library(ggplot2)
library(caret)

WKDIR <- "../../../../../working-bench/ReindeerAppli/res_query_leucegene/"
fname <- "Leucegene-Fusion_Selection_InAML.csv"
# fname <- "query_results_BEAT-AML.csv"

cmp.tab <- read.csv(paste0(WKDIR, fname),
                    header = TRUE, colClasses = "character") %>%
    dplyr::select_if(function(x) !all(x == 0)) %>%
    dplyr::select(-SUM) %>%
    pivot_longer(cols = -seq_name, names_to = "sample", values_to = "query.res") %>%
    dplyr::mutate(`Transipedia | Ground Truth` = lapply(query.res,
                                                        FUN = function(x) {
                                                            if (x %in% c("0", "1", "2")) return ("no | FALSE")
                                                            else if (x == "FP") return ("yes | FALSE")
                                                            else if (x == "FN") return ("no | TRUE")
                                                            else return ("yes | TRUE")
                                                        }) %>% unlist())
plt <- ggplot(data = cmp.tab) +
    geom_tile(aes(y = seq_name, x = sample, fill = `Transipedia | Ground Truth`),
              lwd = 0.5, color = "gray") +
    scale_fill_manual(values = c("yes | TRUE" = "royalblue", "yes | FALSE" = "#e66101",
                                 "no | TRUE" = "#fdb863", "no | FALSE" = "white")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "top", text = element_text(size = 15))
ggsave(plt, filename = paste0(WKDIR, str_replace(fname, "csv", "pdf")),
       width = 15, height = 5, dpi = 300,
       limitsize = FALSE)

print(confusionMatrix(factor(str_split_i(cmp.tab$`Transipedia | Ground Truth`, " \\| ", 1) == "yes"),
                      factor(str_split_i(cmp.tab$`Transipedia | Ground Truth`, " \\| ", 2) == "TRUE"),
                      positive = "TRUE", mode = "everything"))
