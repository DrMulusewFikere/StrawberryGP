## Content
Data for Manuscript: "Accounting for Population Structure in Genomic Prediction of Strawberry Sweetness at a Global Scale"

This repository contains the refined genotype, map file, and phenotypic data used in our manuscript submitted to Sci Report.

Requirements
The analysis requires several R packages, most notably ASReml-R.

Data and Code
Data: The folder Data contains the phenotypic, genotypic, and SNP map data.
Code and Output: Within the Data folder, you'll find a tar-zipped file that includes the code and output (workspace).
Contact
For additional information, please contact: mulusewfikere@gmail.com

## Imputation and population
After initila quality control. The "strawberryImputation" code implements a genotype imputation pipeline for structured strawberry populations. Initially, it loads the required R libraries (e.g., snpReady, ggplot2, reshape2, dplyr) and imports population structure results along with a genetic map. Genotype data for two subpopulations are reformatted by separating sample IDs and SNPs, followed by transformation into long format to facilitate visualization. Missing SNP call distributions across the genome are then visualized using ggplot2, with missing data plotted by genomic intervals and chromosomes for both populations. The imputation is carried out separately for each subpopulation and also jointly across the full dataset using the raw.data() function from the snpReady package. The imputation strategy includes quality filtering based on minor allele frequency and call rate, and uses either mean or KNN-based methods. Concordance between the jointly-imputed dataset and population-specific imputations is assessed by comparing shared SNPs and calculating per-SNP percentage agreement, which is then visualized via density plots. Finally, the script performs a second mean-based imputation to generate a clean 012 matrix format and calculates concordance between the two imputation strategies (KNN vs mean). The workspace is saved for downstream use.
