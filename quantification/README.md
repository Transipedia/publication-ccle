# Validation of `REINDEER` Quantification in Different Situations

This repository contains scripts to evaluate `REINDEER` quantification of arbitrary sequences in diverse situations, reproducing figures in preprint article [[C. Bessière, et. al. (2024)]](https://doi.org/10.1101/2024.02.27.581927).

## 1. Accuracy of RNA expression measure

### 1.1. `SEQC-validation/`: contrast to qRT-PCR measurement

The folder `SEQC-validation/` contains scripts comparing gene expression measured from REINDEER query vs qRT-PCR measurement. Two scripts are enclused:

- `a1_ReindeerOnly_vs_GroundTruth.R` and `a2_Kmerator-Reindeer_vs_GroundTruth.R` respectively use qRT-PCR measurement to evaluate `REINDEER`'s estimation without/with masking of non-specific k-mers.

These scripts reproduces Fig2 and FigS2 of the article.

### 1.2. `Kallisto-validation/`: contrast to `Kallisto-tximport` estimation

The folder `Kallisto-validation/` contains scripts comparing gene expression measured from REINDEER query vs Kallisto-tximport quantification. Three scripts are enclused:

- `a_tximport.R` estimates gene expression from `Kallisto`'s quantification results on the level of transcript isoforms, utilizing `tximport` package.
- `b1_ReindeerOnly_vs_Kallisto.R` and `b2_Kmerator-Reindeer_vs_Kallisto.R` respectively use `Kallisto-tximport` results (both estimated raw counts and TPM) to evaluate `REINDEER`'s estimation without/with masking of non-specific k-mers.

These scripts reproduces Fig2 and FigS2 of the article.

### 1.3. `slope_fitting/`: effect of read to k-mer conversion

The folder `slope_fitting/` encloses two scripts to investigate effect of conversion from sequence reads to k-mers, contrasting REINDEER estimation of gene expression to the `Kallisto-tximport` estimation of gene raw counts.

- `a_show_by_gene.R` studies cases of 25 best, 25 worst and 25 compromised (ranked from 701 to 725 out of 885 in total, but with acceptable $R^2$ coefficents) genes ranked by the $R^2$ coefficient.
- `b_summarised_fitting.R` summarizes all cases into a single scatter plot.

These scripts reproduces FigS9 of the article.

## 2. Query of transposable elements

The folder `HERV-validation/` includes two scripts `b1_cmp1000ERV.R` and `b2_cmp2REdiscoverTE.R` to compare `REINDEER`'s quantification of transcposable elements to either `Telescope` or `REdiscoverTE`, respectively.

These scripts reproduces Fig4 A-B, FigS6 and FigS7 of the article.
