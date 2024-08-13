% File src/library/scoup/man/..Rd
\name{alignsim}
\alias{alignsim}
\alias{aminoSC}
\alias{aaSCupdate}
\alias{aaSCupdate.discrete}
\alias{aaSCupdate.ou}
\alias{alignsim.discrete}
\alias{alignsim.ou}
\alias{alignsim.omega}
\alias{codonvalues}
\alias{discrete}
\alias{omega}
\alias{ou}
\alias{print.aminoSC}
\alias{print.codonvalues}
\alias{sitesim}
\alias{sitesim.discrete}
\alias{sitesim.ou}
\alias{sitesim.omega}

\title{Simulate Codon Sequence Alignment}
\description{
    Obtain an alignment of codon sequences that have been artificially subjected
    to natural selection, imposed as changes in the fitness landscape along the
    branches of a symmetric evolutionary tree.
}
\usage{
    alignsim(adaptIn, seqIn, modelIn, filename=NA)
}
\arguments{
    \item{adaptIn}{A list of class \code{discrete}, \code{omega} or
        \code{ou}.}
    \item{seqIn}{A list of class \code{seqParameters}.}
    \item{modelIn}{A list of class \code{hbParameters}. Not applicable
        when \code{adaptIn} is an object of class \code{discrete}.}
    \item{filename}{A string that specifies the full path of the file
        that will contain the simulated alignment in NEXUS format. Say it is
        given as \code{"seq.nex"}, a file with that name will be printed in
        the working directory. When set as \code{NA} (default), no file will
        be saved.}
}
\details{
    This is the primary function of the package. Codon sequence alignment be
    simulated in terms of the population genetics paradigm. Fitness landscape
    may be kept static or set to be renewed along the branches of a balanced
    phylogeny based on any of the three available methods: Ornstein-Uhlenbeck,
    frequency-dependent or deterministic.
}
\value{
    A NEXUS format file is saved in the specified (or working) directory. In
    addition, a list that contains the following entries is returned.
    \describe{
        \item{\code{seqs}}{A matrix of integers between 1 and 61. The integers
            are the positions of the simulated codons within an ordered set of
            nucleotide triplets. The rows are the extant sequences and the
            columns are alignment sites.}
        \item{\code{DNDS}}{A matrix of the corresponding site-wise (or
            codon-wise) dN/dS value for all the fitness landscapes utilised in
            the simulation.}
        \item{\code{aInfo}}{A string of text that contains details of the
            parameter values that were used during simulation of the codon
            sequence alignment.}
        \item{\code{cseq}}{A dataframe of the simulated codon sequence
            alignment.}
        \item{\code{seqCOL}}{A DNAStringSet object with colorful sequences.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    Complementary functions that are useful for defining the simulation
    parameters needed to successfully utilise this function. These include,
    (a.) \code{\link{discreteInput}}, (b.) \code{\link{hbInput}},
    (c.) \code{\link{ouInput}}, (d.) \code{\link{seqDetails}} and (d.)
    \code{\link{wInput}}. See also \code{DNAStringSet} in the \code{Biostrings}
    package.
}

\references{
    Sadiq, H. et al. (in progress) \code{scoup}: \emph{Simulate Codon Sequences
    with Darwinian Selection Incorporated as an Ornstein-Uhlenbeck Process}.

    Pages H, Aboyoun P, Gentleman R, DebRoy S (2024). \emph{Biostrings:
    Efficient manipulation of biological strings. R package version 2.72.1},
    https://bioconductor.org/packages/Biostrings.
}

\examples{
    alignsim(wInput(), seqDetails(c(ntaxa=8,nsite=10)), NULL, NA)$seqCOL
    
    alignsim(ouInput(), seqDetails(c(ntaxa=8,nsite=10)), hbInput(), NA)$seqCOL
    
    alignsim(discreteInput(), seqDetails(c(ntaxa=8,nsite=10)), NULL, NA)$seqCOL
}
