% File src/library/scoup/man/..Rd
\name{biTree}
\alias{biTree}
\alias{tree}

\title{Generate a Balanced Bifurcating Evolutionary Tree}
\description{
    Obtain an evolutionary tree that is such that all its internal nodes have
    exactly two offspring and all the branches on the tree have equal length.
}
\usage{
    biTree(ntaxa, bLength, terModel=NA)
}
\arguments{
    \item{ntaxa}{Number of extant taxa. It must be an integer (\code{t})
        that may be expressed as \code{2^m}, where \code{m} is itself a
        positive integer.}
    \item{bLength}{Branch length. All the branches of the generated tree
        will have the same length that is equal to the specified value.}
    \item{terModel}{A text that would be added as suffix to the extant
        taxa names. If set as \code{NA} (default), no suffix will be added.
        This is useful for assigning a model to the leaves in branch-specific
        analyses.}
}
\value{
    \describe{
        \item{\code{tree}}{A bifurcating tree in newick format.}
    }
}
\author{
    Hassan Sadiq
}

\examples{
    biTree(16, 0.01, "{foreground}")
    
    biTree(16, 0.01, " #1")

    biTree(16, 0.01)
}
