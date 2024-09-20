\name{seqParameters-class}
\docType{class}
\alias{seqParameters-class}
\alias{branchL}
\alias{branchL,seqParameters-method}
\alias{details}
\alias{details,seqParameters-method}
\alias{nodes}
\alias{nodes,seqParameters-method}
\alias{phylogeny}
\alias{phylogeny,seqParameters-method}
\alias{show,seqParameters-method}
\alias{sites}
\alias{sites,seqParameters-method}
\alias{taxa}
\alias{taxa,seqParameters-method}

\title{Simulated Codon Sequence Structure}

\description{
    A S4 object that contains information about the structure (that is, size,
    length, etc) of the simulated genetic sequence.
}

\section{Objects from the Class}{
    This is the object class of the output from the \code{\link{seqDetails}}
    function. It is a core input of the \code{\link{alignsim}} function.
    Objects can be created by calls of the form \code{new("seqParameters",
    sites=??, taxa=??, nodes=??, branchL=??, phylogeny=??, details=??)}.
}

\section{Slots}{
    \describe{
        \item{\code{sites}:}{numeric value returned by the \code{sites}
            method.}
        \item{\code{taxa}:}{numeric value returned by the \code{taxa}
            method.}
        \item{\code{nodes}:}{numeric value returned by the \code{nodes}
            method.}
        \item{\code{branchL}:}{numeric value returned by the \code{branchL}
            method.}
        \item{\code{phylogeny}:}{character returned by the \code{phylogeny}
            method.}
        \item{\code{details}:}{character returned by the \code{details}
            method.}
    }
}

\section{Methods}{
    \describe{
        \item{alignsim}{\code{signature(adaptIn="discrete",
            seqIn="seqParameters")}: an option of the primary simulation
            function in the \code{\link{scoup}} package. This setting
            activates the deterministic framework.}
        \item{alignsim}{\code{signature(adaptIn="omega",
            seqIn="seqParameters")}: an option of the primary simulation
            function in the \code{\link{scoup}} package. This setting
            activates the frequency-dependent framework.}
        \item{alignsim}{\code{signature(adaptIn="ou",
            seqIn="seqParameters")}: an option of the primary simulation
            function in the \code{\link{scoup}} package. This setting
            activates the Ornstein-Uhlenbeck framework.}
        \item{branchL}{\code{signature(xo="seqParameters")}: branch length.
            Only balanced evolutionary trees are permitted. Therefore, all
            tree nodes have the same length.}
        \item{details}{\code{signature(xo="seqParameters")}: note that
            contain the important parameter settings that generated the
            corresponding data. It is added as comments to the saved output.}
        \item{nodes}{\code{signature(xo="seqParameters")}: number of internal
            (bifurcating) stages of the balanced phylogeny. An evolutionary
            tree with \eqn{2^\code{m}} extant taxa will have \code{m} nodes.}
        \item{phylogeny}{\code{signature(xo="seqParameters")}: newick string
            of the phylogeny utilised for the codon sequence simulation.}
        \item{show}{\code{signature(object="seqParameters")}: summary
            descriptive details about the corresponding sequence alignment.}
        \item{sites}{\code{signature(xo="seqParameters")}: number of codon
            sites that make up the sequence.}
        \item{taxa}{\code{signature(xo="seqParameters")}: number of leaves on
            the phylogeny.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    Codon sequence simulator \code{\link{alignsim}} and the sequence
    preparatory function \code{\link{seqDetails}}.
}

\examples{
    t0 <- seqDetails()
    sites(t0)
}

\keyword{classes}
