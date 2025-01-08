% File src/library/scoup/man/..Rd
\name{seqDetails}
\alias{seqDetails}

\title{Populate Sequence Alignment Information}
\description{
    Create a \code{seqParameters} object that contains features of the sequence
    that needs to be simulated.
}
\usage{
    seqDetails(seqVector=0)
}
\arguments{
    \item{seqVector}{A named vector that provides characteristics of the
        required sequence alignment. See \bold{Details} for further
        information.}
}
\details{
    If fully specified, \code{seqVector} should be a four-element named vector.
    That is, \code{seqVector = c(ntaxa, nsite, blength, terModel)}. \code{ntaxa}
    should be of the form \eqn{2^m}, where \code{m} is an integer. It
    corresponds to the number of extant taxa, default is 64. \code{nsite}, also
    an integer (default = 250), is the number of codon sites. \code{blength} is
    the length of each branch on the balanced symmetric tree that will be used
    for the simulation (default = 0.10). \code{terModel} is a text that will be
    added as a suffix to the leaf names on the phylogeny (default = \code{NA}).
    It is meant to facilitate assignment of models to the terminal nodes for
    branch-wise selection analyses. The purpose of this function is to
    complement \code{alignsim}.
}
\value{
    A \code{seqParameters} object that contains the following.
    \describe{
        \item{\code{sites}}{Number of alignment sites.}
        \item{\code{taxa}}{Number of extant taxa.}
        \item{\code{nodes}}{Number of internal (bifurcating) stages on the
            evolutionary tree. A tree with \eqn{2^{\code{m}}} leaves will have
            \code{m} internal stages.}
        \item{\code{branchL}}{Length of the branches on the phylogeny.}
        \item{\code{phylogeny}}{Evolutionary tree in newick format.}
        \item{\code{details}}{Text that describes the evolutionary tree.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    The codon sequence simulator \code{\link{alignsim}} and
    \code{\link{biTree}}, the balanced evolutionary tree generator.
}
\examples{
    t0 <- seqDetails()
    sites(t0)

    t1 <- seqDetails(c(ntaxa=16, nsite=10, blength=0.20, terModel=" #1"))
    details(t1)
}
