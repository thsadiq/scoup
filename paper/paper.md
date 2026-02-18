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

Statistical inference models that are used to analyse the degree of the impact
of Darwinian natural selection on observed genetic data, command a healthy
portion of the phylogenetic literature [@gupta2023]. Validation of these
largely codon-based models relies heavily on simulated data. Given the ever
increasing diversity of natural selection inference models that exist
[@yang2007; @arenas2015; @hyphy2020], there is a need for more sophisticated
simulators to match the expanding model complexities.

Bioconductor [@bioc2004] is a leading bioinformatics platform distributing
peer-reviewed R packages. A search of the entries on the platform, in
Version 3.22 on 18 February 2026, with keywords including, `codon`, `mutation`,
`selection`, `simulate`, and `simulation` returned a total of 70 packages
(excluding `scoup`) out of the 2361 available. None of the retrieved entries
was dedicated to codon data simulation for natural selection analyses. Thus,
`scoup` is designed on the basis of the mutation-selection (MutSel) framework
[@halpern1998] as an overdue contribution to the void.

Software and/or packages for simulating molecular protein sequences are
a few in the scientific literature [@peng2015]. Existing simulators tend
to be more suitable for quantitative character evolution. These include,
`ape` [@paradis2019], `ouch` [@cressler2015; @butler2004] and `geiger`
[@pennell2014]. Other extensively used DNA sequence simulators including,
`Seq-Gen` [@rambaut1997], `INDELible` [@fletcher2009], `PhyloSim`
[@sipos2011] and `phangorn` [@klaus2011] are parameterized in accordance
with $\omega$-based models [@goldman1994; @muse1994]. More recent sequence
simulators, such as, `phastSim` [@nicola2022] and `AliSim-HPC` [@nhan2023]
prioritized output capacity. Only few genetic simulators were built upon the
more elaborate MutSel evolutionary concept. These include, `Pyvolve`
[@wilke2015] and `SGWE` [@arenas2014]. To the best of our knowledge, these
existing MutSel friendly simulators are only able to generate data from
static landscapes. With our proposed simulator (`scoup`), it is possible
to generate codon sequences from landscapes that are static or those that
are changing (also known as *seascapes*) [@lassig2009].

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


![\label{sfrrame}**Summarised `scoup` algorithm.** The flowchart
  shows the process for a single substitution event. After each
  substitution event, the process returns to *STEP A*,
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
The seascape setting is then defined as a function of
multiple sets of parameters ($\xi_{1}^{}$, $\xi_{2}^{}$, \ldots, $\xi_{k}^{}$,
for $k \leq$ extant taxa size). Second, the coefficient update ($s_t$) step
is done after every substitution event. In addition, the Ornstein-Uhlenbeck
update process is discretised. In other words, the OU jump sizes are fixed
and pre-specified as an input to the simulation functions.

# Implementation

`scoup` is primarily designed using base functions in `R`. Some important
complementary functions are imported from the
[`Matrix`](doi.org/10.32614/CRAN.package.Matrix) and
[`Biostrings`](doi.org/10.18129/B9.bioc.Biostrings) packages.
  
# Conclusions

We present `scoup`, a R package for codon sequences simulation, where the
evolutionary processes are mirrored more realistically than most existing
simulators. Our framework creatively incorporates the Ornstein-Uhlenbeck
process into the mutation-selection evolutionary model. This attribute
could potentially unlock exciting research avenues that will improve
existing knowledge about the complex interactions of different,
potentially interacting, molecular evolutionary processes.

# Code availability

`scoup` is published for free public use under the GPL-2 license. It is
available for download from the
[Bioconductor platform](doi.org/10.18129/B9.bioc.scoup),
along with detailed documentation and tutorial files. Some additional sample
code are accessible in the
[`tests`](https://github.com/thsadiq/scoup/tree/main/tests) and
[`vignettes`](https://github.com/thsadiq/scoup/tree/main/vignettes)
folders of the package.

# Acknowledgements

We thank Ben Murrell for suggesting modelling varying selection coefficients
with an OU process. Computations were performed using the
[`HPC1`](http://www.sun.ac.za/hpc) facility at Stellenbosch University, South
Africa.


# References
