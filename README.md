# scoup
Simulate Codon Sequences with Darwinian Selection Incorporated as an Ornstein-Uhlenbeck Process

## About
A R package intended for the Bioconductor platform that avails an opportunity to simulate molecular data primarily useful for investigating the effects of Darwinian natural selection. Concepts from the population genomics and the phylogenetics literature were merged so that components of observed selection signatures in genetic sequences may be interrogated more thoroughly. This was achieved by exploiting the explicit use of selection coefficients by Halpern-Bruno in their mutation-selection codon model and the versatility of the stochastic Ornstein-Uhlenbeck framework. A new selection discriminant statistic is also proposed. That is, the ratio of the non-synonymous to synonymous selection coefficient variance (vn/vS). Overall, the package represents a promising addition to the phylogenetics modelling literature in at least two ways. First, the capacity for data generation with respect to the validation of codon models of evolution is extended. Second, an exciting avenue is created for the development of sophisticated evolutionary models by the posited statistic.

## Contents
- DESCRIPTION: file describing the package more elaborately. This file is the only difference between the main branch and the cranR branch. For unknown reasons, BiocCheck complains about the Author and Maintainer fields. On the other hand, R CMD check fails without them!
- inst: a folder containing an unfinalised citation file
- man: a folder where the R markdown files were saved.
- NAMESPACE: what functions and/or classes were imported and/or exported.
- NEWS: Progress and/or some to-do comments
- R: a folder of R scripts that contain the functions that make up the package
- test: a directory R scripts where example application of the package are save. The scripts were executed to generate the results intended for the corresponding journal article. They were also (partially) used for compiling the vignette.
- vignettes: a folder where the markdown and citation files for building the vignette are saved.

## Installation
The folder can be cloned directly from git using `git clone https://github.com/thsadiq/scoup.git`. Then the `R CMD BUILD` and `R CMD INSTALL` commands could be used to build and instal it locally to R.

## Parting Words
Enjoy playing around with the package. Useful comments and/or criticisms are welcome!
