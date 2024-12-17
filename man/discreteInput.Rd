% File src/library/scoup/man/..Rd
\name{discreteInput}
\alias{discreteInput}

\title{Populate Deterministic Seascape Model Parameters}
\description{
    Create an object that has a \code{discrete} class attribute. It is
    particularly useful for defining one of the possible inputs of the main
    simulation function \code{alignsim}, when interested in simulating codon
    sequences that evolve with fitness landscapes that change at every internal
    node.
}
\usage{
    discreteInput(defList=list())
}
\arguments{
    \item{defList}{A list that may contain up to seven named entries. See
        \bold{Details} for further information.}
}
\details{
    If fully specified, \code{defList} will be a list with seven elements. The
    preferred list content include (a.) \code{p02xnodes}: a 2-row matrix with
    rows that are properly named as \dQuote{\code{vNvS}} and
    \dQuote{\code{nsynVar}}. Entries in the \dQuote{\code{vNvS}} row should be
    the ratio of the variance of the non-synonymous selection coefficients to
    the variance of the synonymous coefficients. Entries in the
    \dQuote{\code{nsynVar}} row should be the variance of the non-synonymous
    selection coefficients. The number of the matrix columns will be used to
    determine the number of internal nodes to assume for the simulation
    phylogeny. Each column of the matrix will be used to determine the
    parameters of the sampling distribution where the coefficient updates
    will be sampled at every node. Default is a \eqn{2 \times 4} matrix,
    wherein all the values in the \dQuote{\code{vNvS}} row are equal to 1
    and all the entries in the \dQuote{\code{nsynVar}} row are equal to
    \eqn{10^{-5}_{}}. (b.) \code{technique}: a binary integer that could be
    \code{1} for \code{Gaussian} or \code{2} for \code{Gamma} (default)
    distribution. It informs about the probability distribution to be used
    for updating the coefficients. (c.) \code{pSize}: (default = 1000) is
    the effective population size. (d.) \code{nodeIndex}: a nuisance input
    that is best left unspecified. It is updated within the \code{alignsim}
    operation. (e.) \code{leafModel}: a text that may be used to suffix the
    names of the terminal nodes (default = \code{NA}). (f.) \code{kappa}:
    transition to transversion rate between nucleotide pairs (default = 4).
    (g.) \code{mrate}: average rate at which mutation occurs between DNA pairs
    (default = 0.25). Note that this function was not designed to be used in
    isolation. Its purpose is to complement the \code{alignsim} simulation
    function.
}
\value{
    A \code{discrete} object that contains the following.
    \describe{
        \item{\code{lscape}}{Matrix containing landscape parameters.}
        \item{\code{sampler}}{Name of the sampling distribution used for
            selection coefficient updates.}
        \item{\code{effpop}}{Effective population size.}
        \item{\code{kappa}}{Transition to transversion rate between DNA pairs.}
        \item{\code{hkyMu}}{Average DNA mutation rate.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    \code{\link{alignsim}},\verb{  } \code{\link{aaGamma}},\verb{  }
    \code{\link{aaGauss}},\verb{  } \code{\link{biTree}},\verb{  }
    \code{\link{codonMutation}}.
}
\examples{
    dtest <- discreteInput()
    effpop(dtest)
    lscape(dtest)
    sampler(dtest)
}
