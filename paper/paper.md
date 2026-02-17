---
title: 'scoup: Simulate Codon Sequences with Darwinian
        Selection Incorporated as an Ornstein-Uhlenbeck Process'
tags:
  - R
  - Molecular Evolution
  - Ornstein-Uhlenbeck
  - Phylogenetics
  - Population Genetics
  - Simulation
  - Statistics
authors:
  - name: Hassan Sadiq
    corresponding: true
    orcid: 0000-0003-0192-7134
    affiliation: "1, 2"
  - name: Darren P. Martin
    orcid: 0000-0002-8785-0870
    affiliation: 2
affiliations:
 - name: Department of Statistics and Actuarial Science,
         Stellenbosch University, South Africa
   index: 1
 - name: Institute of Infectious Diseases and Molecular Medicine,
         Division of Computational Biology,
         Department of Integrative Biomedical Sciences,
         University of Cape Town, South Africa
   index: 2
date: 17 February 2026
bibliography: paper.bib
---

# Summary

Genetic analyses of natural selection within and between populations have
increasingly developed along separate paths. The two important genres of
evolutionary biology (i.e. phylogenetics and population genetics) borne
from the split can only benefit from research that seeks to bridge
the gap. Simulation algorithms that combine fundamental concepts from
both genres are important to achieve such unifying objective. 
We introduce `scoup`, a codon sequence simulator that is implemented
in R and hosted on the Bioconductor platform. There is hardly any other
simulator dedicated to genetic sequence generation for natural selection
analyses on the platform. Concepts from the Halpern-Bruno mutation-selection
model and the Ornstein-Uhlenbeck (OU) evolutionary algorithm were creatively
fused such that the end-product is a novelty with respect to computational
genetic simulation. Users are able to seamlessly adjust the model parameters
to mimic complex evolutionary procedures that may have been otherwise
infeasible. For example, it is possible to explicitly interrogate the
concepts of static and changing fitness landscapes with regards to Darwinian
natural selection in the context of codon sequences from multiple
populations.

# Statement of need

Statistical inference of the extent to which Darwinian natural selection has
impacted genetic data, commands a healthy portion of the phylogenetic
literature [@jacques2023]. Validation of these largely codon-based models
relies heavily on simulated data. Given the ever increasing diversity of
natural selection inference models that exist [@yang2007; @hyphy2020],
there is a need for more sophisticated simulators to match the
expanding model complexities.

Bioconductor [@bioc2004] is a leading platform where peer-reviewed
bioinformatic software useful for biological data analyses are hosted. A
search of the entries on the platform, in Version 3.19 on 29 October 2024,
with keywords including, `codon`, `mutation`, `selection`, `simulate`, and
`simulation` returned a total of 72 unique packages out of the 2300 available.
None of the retrieved entries was dedicated to codon data simulation for
natural selection analyses. Thus, `scoup` is designed on the basis of the
mutation-selection (MutSel) framework [@halpern1998] as an overdue
contribution to the void. Software and/or packages for simulating genetic
sequences are also rare in the scientific literature [@gearty2024]. 

# Algorithm

`scoup` is further unique for at least three reasons. First, it incorporates
Darwinian natural selection into the MutSel model in terms of variability of
selection coefficients, an extension of an idea from @spielman2015. Second,
it directly utilises the concept of fitness landscapes. Third, fitness
landscape updates can be executed in either a deterministic or a stochastic
format. The stochastic updates are implemented in terms of the more
biologically amenable, Ornstein-Uhlenbeck (OU) process
[@bartoszek2017; @uhlenbeck1930]. A crude summary of how substitution
events are executed in `scoup` is presented in \autoref{sfrrame}.


![\label{sfrrame}**Summarised `scoup` algorithm.**
  After each substitution event, the process returns to *STEP A*,
  until the input tree length ($\tau \in \mathbf{SEQ}$) is exhausted.
  $\sigma^{2}_{n}=$ variance of amino acid selection coefficients.
  $\sigma^{2}_{s}=$ variance of synonymous codon selection
  coefficients. $\Sigma^{2}_{}=$ OU asymptotic variance.
  $\theta=$ OU mean reversion rate. $\mathbf{SEQ}=$ sequence
  information. $x_{\star}^{}=$ codon. $\mathbf{s}_{\star}^{}=$
  codon selection coefficient vector.](FIG1.pdf)

We highlight two important design choices from \autoref{sfrrame}. First,
we assume that a static fitness landscape is obtained from a single set of
parameters ($\xi$) needed to sample a $20$-element numerical vector of amino
acid selection coefficients (that is, $s_{0}^{}$ in \autoref{sfrrame}). The
coefficients are subsequently used as inputs of the corresponding MutSel model.
This ensured that a seascape setting is then defined as a function of multiple
sets of parameters ($\xi_{1}^{}$, $\xi_{2}^{}$, \ldots, $\xi_{k}^{}$, where
$k \leq$ extant taxa size). Second, the coefficient update ($s_t$) step is
done after every substitution event. In addition, the Ornstein-Uhlenbeck
update process is discretised. In other words, the OU jump sizes are fixed
and pre-specified as an input to the simulation functions.

# Implementation

`scoup` may be installed directly from Bioconductor using the following
`R` code in \autoref{installcode}.

![\label{installcode}**R installation code for `scoup`**. Allows the most
  recently published version of the package to be installed from Bioconductor.
  Development version of the package may be installed by adding
  `BiocManager::install(version="devel")` before the final line.](FIG2.pdf)

A sample code for executing a simulation run with `scoup` is presented in
\autoref{pseudocode}. The code executes a stochastic OU framework on a
balanced phylogeny with $64$ extant taxa.

![\label{pseudocode}**An example R code for simulating a codon sequence
  alignment with `scoup`**. Default values were left unchanged. `Line01`:
  OU adaptation parameters where, $\mu=0$, $\Sigma^{2}_{}=0.01$ and
  $\theta=0.01$. `Line02`: evolution model input where,
  $\mathbf{s} \sim \text{Gamma}(1,\sigma_{n}^{-1})$,
  $\sigma^{2}_{n}=10^{-5}_{}$, $\sigma^{2}_{s}=10^{-5}_{}$ and
  effective population size, $N_{\texttt{e}}^{}=1000$. `Line03`:
  sequence information where, site count is $250$, extant taxa count
  is $64$ and branch length is $0.1$.](FIG3.pdf)
  
# Conclusions

We present `scoup`, a R package that allows for simulation of codon
sequences in a way that is capable of recapitulating the evolutionary
processes of biological systems more realistically than most existing
simulators. Our framework creatively incorporates the Ornstein-Uhlenbeck
process into the mutation-selection evolutionary model. This attribute
could potentially unlock exciting research avenues that will improve
existing knowledge about the complex interactions of different,
potentially interacting, molecular evolutionary processes. In another
unique contribution to the literature, the magnitude of the Darwinian
selection affect on the simulated sequences was controlled with the ratio
of the variances of selection coefficients.

# Code availability

`scoup` is published for free public use under the
GPL-2 license. It is available for download from the
[Bioconductor platform](doi.org/10.18129/B9.bioc.scoup),
along with detailed documentation and tutorial files.
Some additional sample code are accessible in
`tests/` and `vignettes/` folders of the package.

# Acknowledgements

We thank Ben Murrell for suggesting modelling varying selection coefficients
with an OU process. Computations were performed using the
[`HPC1`](http://www.sun.ac.za/hpc) facility at Stellenbosch University, South
Africa.


# References
