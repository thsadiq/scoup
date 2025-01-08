\name{discrete-class}
\docType{class}
\alias{discrete-class}
\alias{aaSCupdate,discrete-method}
\alias{alignsim,discrete,seqParameters-method}
\alias{effpop,discrete-method}
\alias{lscape,discrete-method}
\alias{sampler,discrete-method}
\alias{show,discrete-method}
\alias{sitesim,discrete,numeric-method}

\title{(Deterministic) Episodic Simulation Model Input}

\description{
    Creates an object suitable for use when interested in generating an
    alignment of genetic sequences following the episodic (referenced as,
    ``deterministic'') simulation technique available in the
    \code{\link{scoup}} package.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("discrete",
    lscape=..., sampler=..., nodeIndex=..., psize=..., t3mdl=...)}. Objects
    can also be created straightforwardly with the \code{\link{discreteInput}}
    function.
}

\section{Slots}{
    \describe{
        \item{\code{lscape}:}{numeric matrix returned by the \code{lscape}
            method.}
        \item{\code{sampler}:}{numeric value that can be set as 1 or 2. It
            indicates the probability distribution where the amino acid
            selection coefficients is sampled.}
        \item{\code{nodeIndex}:}{numeric input only relevant for implicit
            execution of the simulation algorithm. It is of no practical
            utility to the end-user.}
        \item{\code{psize}:}{numeric value returned by the \code{effpop}
            method.}
        \item{\code{t3mdl}:}{character input that may be used to specify
            suffix for the leaves on the returned phylogeny. It is intended
            to facilitate inference analyses with external software such as
            \code{PAML} or \code{HyPhy}.}
    }
}

\section{Methods}{
    \describe{
        \item{aaSCupdate}{\code{signature(parameters="discrete")}: background
            function that is not intended for end-use. It updates the amino
            acid selection coefficients intermittently during the sequence
            simulation process.}
        \item{alignsim}{\code{signature(adaptIn="discrete",
            seqIn="seqParameters")}: primary simulation function available in
            the \code{\link{scoup}} package.}
        \item{effpop}{\code{signature(x="discrete")}: effective population
            size.}
        \item{lscape}{\code{signature(x="discrete")}: numerical matrix that
            contains parameters of the fitness landscape. The first row will
            contain the ratio of the variance of the non-synonymous to
            synonymous selection coefficients (\code{vNvS}) and the second row
            will contain the variance of the non-synonymous selection
            coefficients \eqn{\sigma^2_\code{n}}. The number of columns will
            be equal to the number of internal (bifurcating) stages. A
            phylogeny with \eqn{2^{\code{m}}} leaves will have \code{m}
            internal stages.}
        \item{sampler}{\code{signature(x="discrete")}: probability distribution
            where the amino acid selection coefficients should be obtained.}
        \item{show}{\code{signature(object="discrete")}: prints characteristics
            of the corresponding genetic sequence, including the population
            size and the number of extant taxa.}
        \item{sitesim}{\code{signature(parameters="discrete",
            nodeLength="numeric")}: background function that is not to be used
            by an end-user. It generates the DNA data at each site
            independently.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    Objects of this class are returned by the \code{\link{discreteInput}}
    function. This class of objects are also required as inputs of the
    \code{\link{alignsim} function}.
}

\examples{
    dtest <- discreteInput()
    effpop(dtest)
    lscape(dtest)
    sampler(dtest)
}

\keyword{classes}
