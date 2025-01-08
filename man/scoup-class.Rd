\name{scoup-class}
\docType{class}
\alias{scoup-class}
\alias{aInfo}
\alias{aInfo,scoup-method}
\alias{cseq}
\alias{cseq,scoup-method}
\alias{dNdS}
\alias{dNdS,scoup-method}
\alias{seqCOL}
\alias{seqCOL,scoup-method}
\alias{seqs}
\alias{seqs,scoup-method}
\alias{show,scoup-method}

\title{Output from the \code{scoup::alignsim} Genetic Sequence Simulator}

\description{
    Stores the results from a successful implementation of any of the
    simulation algorithms available in the \code{\link{scoup}} package.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("scoup", seqs=...,
    DNDS=..., aInfo=..., cseq=..., seqCOL=...)}.
}

\section{Slots}{
    \describe{
        \item{\code{seqs}:}{numerical matrix returned by the \code{seqs}
            method.}
        \item{\code{DNDS}:}{numerical matrix returned by the \code{dNdS}
            method.}
        \item{\code{aInfo}:}{character phrase returned by the \code{aInfo}
            method.}
        \item{\code{cseq}:}{data frame returned by the \code{cseq} method.}
        \item{\code{seqCOL}:}{\code{DNAStringSet} object returned by the
            \code{seqCOL} method.}
    }
}

\section{Methods}{
    \describe{
        \item{aInfo}{\code{signature(x="scoup")}: details of the parameters
            used to execute the simulation process. This includes, the branch
            length of all the nodes of the balanced phylogeny, the name of the
            probability distribution where the amino acid selection
            coefficients is sampled as well as the (\code{vNvS} &
            non-synonymous) selection parameter set used at each internal node
            (``generation'') stage.}
        \item{cseq}{\code{signature(x="scoup")}: data frame that contains the
            simulated genetic sequence.}
        \item{dNdS}{\code{signature(x="scoup")}: analytical estimates of the
            magnitude of the imposed selection effect. It is calculated
            node-wise or after every substitution event (depending on the
            simulation technique in use) as the ratio of the non-synonymous
            to synonymous substitutions.}
        \item{seqCOL}{\code{signature(x="scoup")}: a \code{DNAStringSet}
            version of the simulated genetic sequence alignment.}
        \item{seqs}{\code{signature(x="scoup")}: expression of the simulated
            sequence as a matrix of integers, where each entry corresponds to
            the position of the associated codon in an an alphabetically
            increasing ordered set of the DNA triplets of the 61 sense codons.
            This is intended to facilitate utility of \code{R}'s concatenation
            function to generate site-partitioned sequence alignment.}
        \item{show}{\code{signature(object="scoup")}: sentence that contains
            the number of codon sites and the number of extant taxa that make
            up the simulated genetic sequence alignment.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    Simulation function \code{\link{alignsim}}. The \code{\link{seqWriter}}
    function is a good complement to facilitate the generation of
    site-partitioned sequences.
}

\examples{
    alignEntry <- seqDetails(c(ntaxa=8,nsite=10))
    dsim <- alignsim(discreteInput(), alignEntry)
    aInfo(dsim)
    cseq(dsim)
}

\keyword{classes}
