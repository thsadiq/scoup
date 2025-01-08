% File src/library/scoup/man/..Rd
\name{discreteInput}
\alias{discreteInput}

\title{Populate (Deterministic) Episodic Seascape Model Parameters}
\description{
    Create an object that has a \code{discrete} class attribute. It is
    particularly useful for defining one of the possible inputs of the main
    simulation function \code{alignsim}, when interested in simulating codon
    sequences that undergo episodic evolution. That is, sequences that evolve
    such that fitness landscapes change at every internal node.
}
\usage{
    discreteInput(defList=list())
}
\arguments{
    \item{defList}{A list that may contain up to five named entries. See
        \bold{Details} for further information.}
}
\details{
    If fully specified, \code{defList} will be a list with five elements. The
    preferred list content include (a.) \code{p02xnodes}: a 2-row matrix with
    rows that are properly named as \dQuote{\code{vNvS}} and
    \dQuote{\code{nsynVar}}. Entries in the \dQuote{\code{vNvS}} row should be
    the ratio of the variance of the non-synonymous selection coefficients to
    the variance of the synonymous coefficients. Entries in the
    \dQuote{\code{nsynVar}} row should be the variance of the non-synonymous
    selection coefficients. The number of the matrix columns will be used to
    determine the number of internal nodes to assume for the phylogeny to use
    for the simulation. Each column of the matrix will be used to determine the
    parameters of the sampling distribution where the coefficient updates
    will be sampled for every node. Default is a \eqn{2 \times 4} matrix,
    wherein all the values in the \dQuote{\code{vNvS}} row are equal to 1
    and all the entries in the \dQuote{\code{nsynVar}} row are equal to
    \eqn{10^{-5}_{}}. (b.) \code{technique}: a binary integer that could be
    \code{1} for \code{Gaussian} or \code{2} for \code{Gamma} (default)
    distribution. It informs about the probability distribution to use
    for updating the coefficients. (c.) \code{pSize}: (default = 1000) is
    the effective population size. (d.) \code{nodeIndex}: a nuisance input
    that is best left unspecified. It is updated within the \code{alignsim}
    operation. (e.) \code{leafModel}: a text that may be used to suffix the
    names of the terminal nodes (default = \code{NA}). Note that this
    function was not designed to be used in isolation. Its purpose is to
    complement the \code{alignsim} simulation function.
}
\value{
    A \code{discrete} object that contains the following.
    \describe{
        \item{\code{lscape}}{Matrix containing landscape parameters.}
        \item{\code{sampler}}{Name of the sampling distribution used for
            selection coefficient updates.}
        \item{\code{effpop}}{Effective population size.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    The primary simulation function \code{\link{alignsim}} as well as the
    two in-built functions for generating amino acid selection coefficients,
    \code{\link{aaGamma}} and \code{\link{aaGauss}}. It may also be useful to
    check \code{\link{biTree}} to understand the tree generating process.
}
\examples{
    dtest <- discreteInput()
    effpop(dtest)
    lscape(dtest)
    sampler(dtest)
}
