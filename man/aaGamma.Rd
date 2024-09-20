% File src/library/scoup/man/..Rd
\name{aaGamma}
\alias{aaGamma}

\title{Obtain Gamma Distributed Amino Acid Selection Coefficients}
\description{
    Obtain a vector of values from the Gamma distribution that may be
    conveniently used as amino acid selection coefficients.
}
\usage{
    aaGamma(vNvS, nsynVar)
}
\arguments{
    \item{vNvS}{Ratio of the variance of the coefficients of the
        non-synonymous codons relative to the variance of the synonymous
        selection coefficients. Can be assigned a value equal to zero to
        eliminate synonymous selection.}
    \item{nsynVar}{A non-negative value that corresponds to the variance
        of the coefficients of the non-synonymous codons. That is, the
        the between-amino-acid variance.}
}
\details{
    Twenty random observations from a \eqn{\code{Gamma}(1,
    1/\sqrt{\code{nsynVar}})} distribution are sampled as the amino acid
    selection coefficients. The associated variance of the synonymous
    coefficients (\code{synVar}) is calculated as \code{nsynVar/vNvS}. If
    \eqn{\code{nsynVar} < 10^{-12}}, all the amino acid coefficients will be
    equal to a single random draw from a \eqn{\code{Gamma}(1,1/\code{synVar})}
    distribution.
}
\value{
    Returns an object of class \code{aminoSC}, a vector that at least contains
    the following components.
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
    An alternative approach for generating amino acid selection coefficients
    from a normal distribution \code{\link{aaGauss}}. There is
    \code{\link{codonCoeffs}} as well, a function designed to convert amino
    acid to codon selection coefficients.
}

\examples{
    test1 <- aaGamma(0.50, 1e-04)
    coeffs(test1)
    synVar(test1)
    nsynVar(test1)
    
    test2 <- aaGamma(1e-02, 0)
    coeffs(test2)
    synVar(test2)
    nsynVar(test2)
}
