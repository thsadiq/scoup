% File src/library/scoup/man/..Rd
\name{scoup}
\alias{scoup}
\alias{scoup-package}
\docType{package}

\title{Simulate Codons with Darwinian Selection Added as an OU Process}

\description{
    The primary objective of this package is to facilitate more rigorous
    understanding of phylogenetic inferences of natural selection from codon
    sequences. Concepts from the Halpern-Bruno mutation-selection model and
    the Ornstein-Uhlenbeck stochastic process were creatively fused such that
    the end-product is a novelty with respect to computational genetic
    simulation. Users are able to seamlessly adjust the model parameters to
    mimic complex evolutionary procedures that may have been otherwise
    infeasible. For example, it is possible to explicitly interrogate the
    concepts of static and changing fitness landscapes with regards to
    Darwinian natural selection in the context of DNA sequences. The ratio of
    the variance in selection coefficients, vN/vS, is presented as a new
    measure of the net selection effect acting on genetic sequences. This
    package could be very useful for generating more appropriate test data
    sets for validation of likelihood-based (\eqn{\omega}) codon models of
    evolution.
}

\details{
    Three simulation algorithms are available. (a.) The Ornstein-Uhlenbeck
    simulation technique. This technique was built around the stochastic
    Brownian motion evolutionary paradigm. Explicit parameters exist to control
    the extent of drift, mutation and selection that are acting on the
    biological system. (b.) The frequency-dependent approach where every
    substitution event that corresponds to a shift in the fitness landscape.
    (c.) The deterministic method where the model parameters may be fixed
    for each internal node of the phylogeny.
}

\author{
    Hassan Sadiq
}

\references{
    Halpern, A. L. and Bruno, W. J. (1998). Evolutionary Distances for
    Protein-Coding Sequences: Modelling Site-Specific Residue Frequencies,
    \emph{Molecular Biology and Evolution} \bold{15}(7): 910-917.

    Sadiq, H. et al. (in progress) \code{scoup}: \emph{Simulate Codon Sequences
    with Darwinian Selection Incorporated as an Ornstein-Uhlenbeck Process}.

    Uhlenbeck, G. E. and Ornstein, L. S. (1930), On the Theory of the Brownian
    Motion, \emph{Physical Review} \bold{36}: 823-841.
}
