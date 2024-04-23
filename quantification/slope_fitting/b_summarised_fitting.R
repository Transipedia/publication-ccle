rm(list = ls())

library(stringr)
library(ggplot2)
library(ggpointdensity)
library(ggExtra)

# WKDIR <- "~/Projects/working-bench/ReindeerAppli/"

CMP_TAB_PATH <- cmd.Args[1] # res_cmp2Kallisto/reindeer-kmerator2kallisto/kmerator-reindeer.tsv
OUT_DIR <- cmd.Args[2] # fitting_slope/
if (!dir.exists(OUT_DIR)) {
    dir.create(OUT_DIR, recursive = TRUE)
    cat("Created output folder:", OUT_DIR, "\n")
}

mthd <- "Sum.reindeer"
count_column <- "Kallisto.counts"

cmp.tab <- read.table(CMP_TAB_PATH, header = TRUE, sep = "\t")
cmp.tab.pos <- cmp.tab[(cmp.tab[, mthd] > 0) & (cmp.tab[, count_column] > 0), ]
fit.tab <- cmp.tab.pos %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(mean.Sum.reindeer = mean(Sum.reindeer),
                     sd.Sum.reindeer = sd(Sum.reindeer),
                     mean.Kallisto.counts = mean(Kallisto.counts),
                     sd.Kallisto.counts = sd(Kallisto.counts),
                     slope = lm(formula = Sum.reindeer ~ Kallisto.counts)$coefficients[2],
                     R2 = summary(lm(formula = Sum.reindeer ~ Kallisto.counts))$adj.r.squared,
                     n_smp = dplyr::n())
write.table(fit.tab, paste0(OUT_DIR, "fit_tab2.tsv"),
            sep="\t", row.names = FALSE)

theme_set(theme_bw())
plot_glb <- ggplot(fit.tab) +
    geom_point(aes(x = slope, y = R2), size = 1, color = "royalblue") +
    geom_vline(xintercept = 140, color = "orange", size = 1) +
    theme(text = element_text(size = 15), legend.position = "none")
plot_glb <- ggMarginal(plot_glb, type="histogram", bins = 100,
                       fill = "royalblue", color = "royalblue")
ggsave(plot_glb, filename = paste0(OUT_DIR, "scatter_R2_slope.global.pdf"),
       width = 5, height = 5)

fit.tab.sub <- fit.tab %>% 
    dplyr::filter(R2 >= 0.95)
plot_lcl <- ggplot(fit.tab.sub) +
    geom_point(aes(x = slope, y = R2), size = 1, color = "royalblue") +
    geom_vline(xintercept = 140, color = "orange", size = 1) +
    theme(text = element_text(size = 15), legend.position = "none")
plot_lcl <- ggMarginal(plot_lcl, type="histogram", bins = 100,
                       fill = "royalblue", color = "royalblue")
ggsave(plot_lcl, filename = paste0(OUT_DIR, "scatter_R2_slope.zoomed.pdf"),
       width = 5, height = 5)
