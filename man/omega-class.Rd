\name{omega-class}
\docType{class}
\alias{omega-class}
\alias{alignsim,omega,seqParameters-method}
\alias{kappa,omega-method}
\alias{hky85mu,omega-method}
\alias{effpop,omega-method}
\alias{lscape,omega-method}
\alias{nsynVar,omega-method}
\alias{sampler,omega-method}
\alias{show,omega-method}
\alias{sitesim,omega,numeric-method}
\alias{vNvS,omega-method}

\title{Frequency-Dependent Evolutionary Model Specification}

\description{
    Creates an object that contains the inputs that are necessary to define
    a frequency-dependent evolutionary algorithm and subsequently simulate
    genetic sequence alignment based on the framework using the
    \code{\link{alignsim}} function in the \code{\link{scoup}} package.
}

\section{Objects from the Class}{
    Objects of this class can be created by calls of the form
    \code{new("omega", nsynVar=..., psize=..., sampler=..., aaPlus=...,
    vNvS=..., kappa=..., mrate=...)}. The object is an important input of
    the \code{\link{alignsim}} function when interested in simulating
    sequences with respect to the frequency-dependent framework. The
    \code{\link{wInput}} function in the \code{\link{scoup}} package
    returns this kind of object.
}

\section{Slots}{
    \describe{
        \item{\code{nsynVar}:}{numeric value returned by the \code{nsynVar}
            method.}
        \item{\code{psize}:}{numeric value returned by the \code{effpop}
            method.}
        \item{\code{sampler}:}{numeric value that can be set as 1 or 2. It
            indicates the probability distribution where the amino acid
            selection coefficients should be sampled.}
        \item{\code{aaPlus}:}{indices of the amino acids (after the
            corresponding one-letter IUPAC names are arranged in increasing
            alphabetical order) that should be assigned non-zero \code{vNvS}.}
        \item{\code{vNvS}:}{numeric value returned by the \code{vNvS}
            method.}
        \item{\code{kappa}:}{numeric value returned by the \code{kappa}
            method.}
        \item{\code{mrate}:}{numeric value returned by the \code{hky85mu}
            method.}
    }
}

\section{Methods}{
    \describe{
        \item{alignsim}{\code{signature(adaptIn="omega",
            seqIn="seqParameters")}: primary simulation function availed in
            the \code{\link{scoup}} package.}
        \item{effpop}{\code{signature(x="omega")}: effective population size.}
        \item{lscape}{\code{signature(x="omega")}: IUPAC one-letter notations
            of the amino acids that were assigned non-zero \code{vNvS} values.}
        \item{nsynVar}{\code{signature(x="omega")}: variance of the
            non-synonymous selection coefficients \eqn{\sigma^2_\code{n}}.}
        \item{sampler}{\code{signature(x="omega")}: probability
            distribution where the amino acid selection coefficients should
            be randomly retrieved.}
        \item{show}{\code{signature(object="omega")}: prints the \code{vNvS}
            and the count of amino acids that had positive \code{vNvS} values.}
        \item{sitesim}{\code{signature(parameters="omega",
            nodeLength="numeric")}: background function that is not available
            to end-user. It generates the DNA data at each site independently.}
        \item{vNvS}{\code{signature(x="omega")}: ratio of the variance of the
            non-synonymous to synonymous selection coefficients.}
        \item{kappa}{\code{signature(x="omega")}: ratio of the rate of the
            transition to transversion of nucleotides.}
        \item{hky85mu}{\code{signature(x="omega")}: average nucleotide mutation
            rate.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    \code{\link{alignsim}}, \verb{  } \code{\link{wInput}}
}

\examples{
    w1 <- wInput(list(aaPlus=c(4,2,11), nsynVar=10))
    lscape(w1)
    w1
}

\keyword{classes}
