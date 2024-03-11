rm(list = ls())

library(stringr)
library(ggplot2)
library(ggpointdensity)
library(patchwork)

WKDIR <- "~/Projects/working-bench/ReindeerAppli/"

mthd <- "Sum.reindeer"
count_column <- "Kallisto.counts"

# Load Kalisto counts 
cmp.tab <- read.table(paste0(WKDIR, "res_cmp2Kallisto/reindeer-kmerator2kallisto/kmerator-reindeer.tsv"),
                      header = TRUE, sep = "\t")
cmp.tab.pos <- cmp.tab[(cmp.tab[, mthd] > 0) & (cmp.tab[, count_column] > 0), ]

# Fit slopes by gene ID
fit.tab <- cmp.tab.pos %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise("slope" = lm(formula = Sum.reindeer ~ Kallisto.counts)$coefficients[2],
                     "R2" = summary(lm(formula = Sum.reindeer ~ Kallisto.counts))$adj.r.squared,
                     "n_smp" = dplyr::n())
write.table(fit.tab, paste0(WKDIR, "fitting_slope/fit_tab.tsv"), sep="\t", row.names = FALSE)

# Figure 1: histogram of adjusted R2
plt1 <- ggplot(na.omit(fit.tab)) +
    geom_histogram(aes(x = R2), binwidth = 0.05, fill = "royalblue", center=0) +
    geom_vline(xintercept = 0.75, color = "navy", linetype = "dashed") +
    geom_text(aes(x = 0.75, y = 600, label = "adj.R2 = 0.75"), color = "navy") +
    xlab("adjusted R-squared") +
    theme_bw() +
    theme(text = element_text(size=15))
ggsave(plt1, filename = paste0(WKDIR, "fitting_slope/histoplot_R2_distribution.pdf"),
       width = 9, height = 3)

# Remove fitting results by less than 5 samples
fit.tab.2plot <- fit.tab %>%
    dplyr::filter(n_smp >= 5) # from 922 genes to 885 genes

# Plot best fitting cases
top_genes <- fit.tab.2plot$ID[order(fit.tab.2plot$R2, decreasing = TRUE)[1:25]]
plt_lst <- lapply(top_genes, FUN = function(gene2sel) {
    cmp.tab.pos.sub <- cmp.tab.pos[cmp.tab.pos$ID == gene2sel, ]
    lm.res <- lm(formula = cmp.tab.pos.sub$Sum.reindeer ~ cmp.tab.pos.sub$Kallisto.counts)
    ggplot(cmp.tab.pos.sub) +
        geom_point(aes(x = Kallisto.counts, y = Sum.reindeer)) +
        geom_text(aes(x = Inf, y = 0,
                      label = paste0("slope = ", round(lm.res$coefficients[2], 2), "\n",
                                     "adj.R2 = ", round(summary(lm.res)$adj.r.squared, 2))),
                  hjust = 1.1, vjust = 0.1) +
        ggtitle(gene2sel) +
        theme_bw() +
        theme(text = element_text(size = 12))
})
ggsave(wrap_plots(plt_lst, ncol = 5),
       filename = paste0(WKDIR, "fitting_slope/top25_genes_fitted.pdf"),
       width = 15, height = 15)

# Plot OK fitting cases
ok_genes <- fit.tab.2plot$ID[order(fit.tab.2plot$R2, decreasing = TRUE)[701:725]]
plt_lst <- lapply(ok_genes, FUN = function(gene2sel) {
    cmp.tab.pos.sub <- cmp.tab.pos[cmp.tab.pos$ID == gene2sel, ]
    lm.res <- lm(formula = cmp.tab.pos.sub$Sum.reindeer ~ cmp.tab.pos.sub$Kallisto.counts)
    ggplot(cmp.tab.pos.sub) +
        geom_point(aes(x = Kallisto.counts, y = Sum.reindeer)) +
        geom_text(aes(x = Inf, y = 0,
                      label = paste0("slope = ", round(lm.res$coefficients[2], 2), "\n",
                                     "adj.R2 = ", round(summary(lm.res)$adj.r.squared, 2))),
                  hjust = 1.1, vjust = 0.1) +
        ggtitle(gene2sel) +
        theme_bw() +
        theme(text = element_text(size = 12))
})
ggsave(wrap_plots(plt_lst, ncol = 5),
       filename = paste0(WKDIR, "fitting_slope/rank701to725_genes_fitted.pdf"),
       width = 15, height = 15)

# Plot worst fitting cases
bottom_genes <- fit.tab.2plot$ID[order(fit.tab.2plot$R2, decreasing = FALSE)[1:25]]
plt_lst <- lapply(bottom_genes, FUN = function(gene2sel) {
    cmp.tab.pos.sub <- cmp.tab.pos[cmp.tab.pos$ID == gene2sel, ]
    lm.res <- lm(formula = cmp.tab.pos.sub$Sum.reindeer ~ cmp.tab.pos.sub$Kallisto.counts)
    ggplot(cmp.tab.pos.sub) +
        geom_point(aes(x = Kallisto.counts, y = Sum.reindeer)) +
        geom_text(aes(x = Inf, y = 0,
                      label = paste0("slope = ", round(lm.res$coefficients[2], 2), "\n",
                                     "adj.R2 = ", round(summary(lm.res)$adj.r.squared, 2))),
                  hjust = 1.1, vjust = 0.1) +
        ggtitle(gene2sel) +
        theme_bw() +
        theme(text = element_text(size = 12))
})
ggsave(wrap_plots(plt_lst, ncol = 5),
       filename = paste0(WKDIR, "fitting_slope/bottom25_genes_fitted.pdf"),
       width = 15, height = 15)