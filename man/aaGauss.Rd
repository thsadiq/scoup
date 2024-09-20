% File src/library/scoup/man/..Rd
\name{aaGauss}
\alias{aaGauss}

\title{Obtain Gaussian Distributed Amino Acid Selection Coefficients}
\description{
    Obtain a vector of values from a Normal distribution that may be
    conveniently used as amino acid selection coefficients.
}
\usage{
    aaGauss(vNvS, nsynVar)
}
\arguments{
    \item{vNvS}{Ratio of the variance of the selection coefficients
        of the non-synonymous codons relative to the variance of the
        synonymous coefficients. It can be assigned a value equal to zero
        to eliminate synonymous selection.}
    \item{nsynVar}{A non-negative value that corresponds to the variance
        of the coefficients among the non-synonymous codons. That is, the
        between-amino-acid variance.}
}
\details{
    An observation is sampled from a \code{Normal(0,nsynVar)} distribution
    independently for each of the 20 amino acid residues. The variance of the
    synonymous selection coefficients (\code{synVar}) is calculated as
    \code{nsynVar/vNvS}. If \code{nsynVar} is less than \eqn{10^{-12}}, all the
    amino acid coefficients will be equal to a single random draw from a
    \code{Normal(0,synVar)} distribution.
}
\value{
    Returns an object of class \code{aminoSC}, a vector that at least contains
    the following component.
    \describe{
        \item{\code{coeffs}}{A vector of 20 numeric elements that represent the
            sampled amino acid coefficients. The coefficients are ordered in
            terms of the 1-letter amino acid IUPAC labels. That is, (A, C,
            . . . , W, Y).}
        \item{\code{synVar}}{Variance of the selection coefficients of the
            synonymous codons.}
        \item{\code{nsynVar}}{Variance of the selection coefficients of the
            non-synonymous codons.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    An alternative sampling function, \code{\link{aaGamma}} is also available.
    The \code{\link{codonCoeffs}} function requires the output from this
    function (or from \code{\link{aaGamma}}).
}

\examples{
    case0 <- aaGauss(0.50, 1e-04)
    nsynVar(case0)
    synVar(case0)
    coeffs(case0)
    
    case1 <- aaGauss(1e-02, 0)
    nsynVar(case1)
    synVar(case1)
    coeffs(case1)
}
