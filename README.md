# scoup
Simulate Codon Sequences with Darwinian Selection Incorporated as an Ornstein-Uhlenbeck Process

## About
A R package intended for the Bioconductor platform that avails an opportunity to simulate molecular data primarily useful for investigating the effects of Darwinian natural selection. Concepts from the population genomics and the phylogenetics literature were merged so that components of observed selection signatures in genetic sequences may be interrogated more thoroughly. This was achieved by exploiting (a.) the explicit use of selection coefficients by Halpern-Bruno in their popular 1998 mutation-selection codon model and (b.) the versatility of the stochastic Ornstein-Uhlenbeck algorithm. A new selection discriminant statistic is also proposed. That is, the ratio of the non-synonymous to synonymous selection coefficient variance (vn/vS). Overall, the package represents a promising addition to the phylogenetics modelling literature in at least two ways. First, the capacity for data generation with respect to the validation of codon models of evolution is extended. Second, an exciting avenue is created for the development of sophisticated evolutionary models by the posited statistic.

## Contents
- DESCRIPTION: file describing the package more elaborately.
- inst: a folder containing an unfinalised citation file
- LICENSE: an appropriate license file obtain from Cran-R
- man: a folder where the R markdown files were saved.
- NAMESPACE: what functions and/or classes were imported and/or exported.
- NEWS: Progress and/or some to-do comments
- R: a folder of R scripts that contain the functions that make up the package
- test: a directory of R scripts where sample applications of the package are saved. The scripts were executed to generate the results intended for the accompanying journal article. They were also (minimally) used for compiling the vignette.
- vignettes: a folder where the markdown and citation files for building the vignette are saved.

## Installation
Upon acceptance on the Bioconductor platform, `scoup` can be installed from the platform with the following command:

``` r
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
 BiocManager::install("scoup")
```
