\name{hbParameters-class}
\docType{class}
\alias{hbParameters-class}
\alias{kappa,hbParameters-method}
\alias{hky85mu,hbParameters-method}
\alias{effpop,hbParameters-method}
\alias{nsynVar,hbParameters-method}
\alias{sampler,hbParameters-method}
\alias{show,hbParameters-method}
\alias{vNvS,hbParameters-method}

\title{Halpern-Bruno Mutation-Selection Evolutionary Model Input}

\description{
    Creates an object of class (\code{hbParameters}) that contains the
    principal entries necessary to construct a Halpern-Bruno mutation-selection
    evolutionary model.
}

\section{Objects from the Class}{
    Objects of this class can be created by calls of the form
    \code{new("hbParameters", psize=??, vNvS=??, sampler=??, nsynVar=??,
    kappa=??, mrate=??, words=??)}. The object is an important input of
    the \code{\link{alignsim}} function when interested in simulating
    sequences with respect to the Ornstein-Uhlenbeck framework.
}

\section{Slots}{
    \describe{
        \item{\code{psize}:}{numeric value returned by the \code{effpop}
            method.}
        \item{\code{vNvS}:}{numeric value returned by the \code{vNvS}
            method.}
        \item{\code{sampler}:}{numeric value that can be set as 1 or 2. It
            indicates the probability distribution where the amino acid
            selection coefficients should be sampled.}
        \item{\code{nsynVar}:}{numeric value returned by the \code{nsynVar}
            method.}
        \item{\code{kappa}:}{numeric value returned by the \code{kappa}
            method.}
        \item{\code{mrate}:}{numeric value returned by the \code{hky85mu}
            method.}
        \item{\code{words}:}{comments about the specified Halpern-Bruno model
            parameters. It is a \code{character} format string that is
            eventually added to the simulated sequence for reference.}
    }
}

\section{Methods}{
    \describe{
        \item{effpop}{\code{signature(x = "hbParameters")}: effective
            population size.}
        \item{nsynVar}{\code{signature(x = "hbParameters")}: variance of the
            non-synonymous selection coefficients \eqn{\sigma^2_\code{n}}.}
        \item{sampler}{\code{signature(x = "hbParameters")}: probability
            distribution where the amino acid selection coefficients should
            be randomly retrieved.}
        \item{show}{\code{signature(object = "hbParameters")}: prints
            characteristics of the defined model including the population size,
            the \code{vNvS} and the \eqn{\sigma^2_\code{n}}.}
        \item{vNvS}{\code{signature(x = "hbParameters")}: ratio of the variance
            of the non-synonymous to synonymous selection coefficients.}
        \item{kappa}{\code{signature(x = "hbParameters")}: ratio of the rate of
            the transition to transversion of nucleotides.}
        \item{hky85mu}{\code{signature(x = "hbParameters")}: average nucleotide
            mutation rate.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    \code{\link{alignsim}}, \verb{   }\code{\link{hbInput}}
}

\examples{
    h1 <- hbInput(c(Ne=100, meth=2, vNvS=1e-08, nsynVar=1e-08))
    sampler(h1)
    h1
}

\keyword{classes}
