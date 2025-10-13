## Content
### Manuscript Title
"Accounting for Population Structure in Genomic Prediction of Strawberry Sweetness at a Global Scale"

This repository contains the input/output (I/O) data and pipeline used to replicate the study. The manuscript has been submitted to *Scientific Reports*.

### Requirements
The analysis requires several R packages, most notably *ASReml-R*, *FIMPUTE*, *ADMIXTURE*.

### Input and pipeline
Data: The folder Data contains the phenotypic, genotypic, and SNP map data.
For additional information, please contact: mulusewfikere@gmail.com

## 1. Imputation and population structure
Following initial quality control, the genetic structure of the study population was evaluated using agglomerative hierarchical clustering (Ward’s method with Euclidean distance) based on SNP markers that passed the filtering criteria. The clustering was performed on a matrix of reference allele counts per accession using R. The resulting dendrogram was visualized with testing locations mapped to each accession. The optimal number of genetic clusters (k) was determined using the silhouette method, which assesses the cohesiveness of individuals within clusters and their separation from other clusters by computing average silhouette widths across a range of k values. The value of k that maximized the average silhouette width was selected as optimal. To further investigate population structure, ancestry proportions were estimated using ADMIXTURE with 20-fold cross-validation, and the k value corresponding to the lowest cross-validation error was selected.

Genotype imputation was performed using a deterministic algorithm implemented in *FImpute v3*, both across the entire population and within the sub-populations defined through clustering, to assess the impact of population structure on imputation accuracy. Accuracy was evaluated by masking 2,000 known genotypes at random and comparing the imputed values to the original true genotypes. This process was repeated ten times, and imputation accuracy was summarized using the mean Pearson correlation between true and imputed genotypes. Additionally, the concordance rate—defined as the proportion of correctly imputed SNPs—was calculated to compare the consistency between population-wide and sub-population-specific imputation strategies.

## 2. Fitting statistical models
To investigate the impact of population structure on genomic prediction accuracy, three approaches were used to construct the additive genomic relationship matrix (GRM). These included the standard GBLUP method using average allele frequencies across all populations, a re-parameterized GBLUP model using eigenvalue decomposition of the GRM, and a structure-corrected GRM based on population-specific allele frequencies. Modeling was performed in R using ASReml-R v4.0, where within-location genomic environments were defined by trial-season combinations with homogeneous variances. Initial models were fitted separately by location and trial to identify the most parsimonious structure based on *Akaike Information Criterion (AIC)*, followed by combined multi-trial models. A final multi-trial model using a factor analytic (FA) approach was fitted to account for genotype-by-environment interactions, with FA1 to FA3 models compared to select the best-fitting and most parsimonious option.

## 3. Description of the Pipeline for Model Implementation <sub><sup>🚧 *Uploading in Progress*</sup></sub>
Below, we describe the pipeline implementation and input/output (I/O) structure.
#### 3.1. QC of Raw Data
#### 3.2. Genomic Relationship Matrix and Population Structure
#### 3.3. Individual Trial Model
#### 3.4. Multi-Trial and Factor Analytic (FA) Model
  ##### *3.4.1. Standard GBLUP Model (Gfa)*
  ##### *3.4.2. Incorporating PCA Eigenvalues and Re-Parameterization of the GBLUP Method (Pfa)*
  ##### *3.4.3. Multi-Population GBLUP Approach Using Sub-Population Genomic Relationship Matrix (Wfa)*
## 4. R Code and Command-Line Scripts (HPC Linux Environment) for Generating Figures <sub><sup>🚧 *Uploading in Progress*</sup></sub>
#### 4.1. Marker Density and Distribution of SNPs Across the Strawberry Genome
#### 4.2. Constructing Genomic Relationship Matrix
#### 4.3. Principal Coordinate Analysis (PCoA) and Selection of K Values
#### 4.4. Population Structure Using ADMIXTURE and Cross-Validation Error for K Selection
#### 4.5. Variance Component and Prediction Accuracy Plots
