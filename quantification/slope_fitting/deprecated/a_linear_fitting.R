rm(list = ls())

library(tidyr)
library(stringr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"
INDIR <- paste0(WKDIR, "data/1000ERV-56samples/")
OUTDIR <- paste0(WKDIR, "res_HERV1000/")

mthd <- "Sum.reindeer"

fit_plot <- function(filepath, count_column) {
    cmp.tab <- read.table(paste0(WKDIR, filepath), header = TRUE, sep = "\t")
    cor.tab <- data.frame(row.names = mthd,
                          corr.pearson.raw = cor(cmp.tab[, count_column], cmp.tab[, mthd], method = "pearson"),
                          corr.spearman.raw = cor(cmp.tab[, count_column], cmp.tab[, mthd], method = "spearman"))
    cmp.tab.pos <- cmp.tab[(cmp.tab[, mthd] > 0) & (cmp.tab[, count_column] > 0), ]
    fit.res <- lm(formula = cmp.tab.pos[, mthd] ~ cmp.tab.pos[, count_column])
    
    plt <- ggplot() +
        geom_smooth(aes(x = cmp.tab.pos[, count_column] + 1,
                        y = cmp.tab.pos[, mthd] + 1), method = "lm", color="orange") +
        geom_pointdensity(aes(x = cmp.tab[, count_column] + 1,
                              y = cmp.tab[, mthd] + 1),
                          size = 0.1) +
        geom_text(aes(x = Inf, y = 1, 
                      label = paste0("Pearson: ", round(cor.tab[mthd, "corr.pearson.raw"], 2), "\n",
                                     "Spearman: ", round(cor.tab[mthd, "corr.spearman.raw"], 2), "\n",
                                     "Slope: ", round(fit.res$coefficients[2], 2))),
                  size = 5, hjust = 1, vjust = 0) +
        # scale_x_log10(breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001, 10000001),
        #               labels = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
        # scale_y_log10(breaks = c(1, 11, 101, 1001, 10001, 100001, 1000001, 10000001),
        #               labels = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
        xlab(count_column) +
        ylab(paste0(mthd)) +
        scale_colour_gradient(low = "#56B1F7", high = "#132B43") +
        theme_bw() +
        theme(text = element_text(size = 15), legend.position = "none")
    return(plt)
}

plt1 <- fit_plot("res_cmp2Kallisto/reindeer-kmerator2kallisto/kmerator-reindeer.tsv", "Kallisto.counts")
plt2 <- fit_plot("res_HERV1000/TE_vs_Kmerator-Reindeer.4metrics.csv", "Counts.Telescope")
plt3 <- fit_plot("res_HERV1000/REdiscoverTE_vs_Kmerator-Reindeer.4metrics.csv", "count.REdiscoverTE")

ggsave(plt1 / plt2 / plt3,
       filename = paste0(WKDIR, "Fitted_sum2raw.kallisto-Telescope-REdiscoverTE.nolog.pdf"),
       width = 4.5, height = 7.5, dpi = 600)
