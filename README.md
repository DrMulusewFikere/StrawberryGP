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
After initial quality control. To evaluate the genetic structure of the study population, a cluster analysis was conducted using agglomerative hierarchical clustering (Ward’s method with Euclidean distance) based on SNP markers that passed quality filtering. The clustering was performed on a matrix of reference allele counts per accession using R. The resulting dendrogram was visualized with testing locations mapped to each accession. The optimal number of genetic clusters (k) was determined using the silhouette method, which assesses the cohesiveness of individuals within clusters and their separation from other clusters by computing average silhouette widths across a range of k values. The value of k that maximized the average silhouette width was selected as optimal. To further investigate population structure, ancestry proportions were estimated using ADMIXTURE with 20-fold cross-validation, and the k value corresponding to the lowest cross-validation error was selected.

Genotype imputation was performed using a deterministic algorithm implemented in FImpute v3, both across the entire population and within the sub-populations defined through clustering, to assess the impact of population structure on imputation accuracy. Accuracy was evaluated by masking 2,000 known genotypes at random and comparing the imputed values to the original true genotypes. This process was repeated ten times, and imputation accuracy was summarized using the mean Pearson correlation between true and imputed genotypes. Additionally, the concordance rate—defined as the proportion of correctly imputed SNPs—was calculated to compare the consistency between population-wide and sub-population-specific imputation strategies.
