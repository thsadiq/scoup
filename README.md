# scoup
Simulate Codon Sequences with Darwinian Selection Incorporated as an Ornstein-Uhlenbeck Process

## About
A R package intended for the Bioconductor platform. It avails an opportunity to simulate molecular data that are primarily useful for investigating the effects of Darwinian natural selection. Concepts from the population genomics and the phylogenetics literature are merged so that components of observed selection signatures in genetic sequences may be interrogated more thoroughly. This is achieved by combining (a.) the explicit use of selection coefficients by Halpern-Bruno in their popular 1998 mutation-selection codon model and (b.) the versatility of the stochastic Ornstein-Uhlenbeck algorithm. Selection pressure is controlled with a new feature of the evolutionary biology research. That is, the ratio of the variance of non-synonymous to synonymous selection coefficients (vN/vS). Overall, the package represents a promising addition to the phylogenetics modelling literature in at least two ways. First, it extends the capacity to generate molecular data to validate codon models of evolution. It offers an uncommon opportunity to thoroughly examine the implications of different fitness landscape settings and the interactions among different evolutionary factors. Second, it creates an exciting avenue for the development of sophisticated evolutionary models that are functions of the introduced vN/vS statistic.

## Contents
- DESCRIPTION: file describing the package more elaborately.
- inst: a folder containing a temporary citation file
- man: a folder where the R markdown files are saved.
- NAMESPACE: what functions and/or classes are imported and/or exported.
- NEWS: Details of the edits, bug fixes and updates made to the package.
- R: a folder of R scripts that contains different (back-end) functions of the package.
- test: a directory of R scripts where sample applications of the package are saved. The scripts are intended to assist users get familiar with the package. They are also (minimally) used within the vignette.
- vignettes: a folder where the markdown and citation files for building the vignette are saved.

## Installation
`scoup` can be installed from the Bioconductor platform with the following command. Interested users should visit the `devel` section on the Bioconductor website for access to the developmental version of the package that include the most recent updates and bug fixes.

``` r
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
 BiocManager::install("scoup")
```
