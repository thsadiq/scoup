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
  % \describe{
    \item{vNvS}{Ratio of the variance of the coefficients of the
      non-synonymous codons relative to the variance of the synonymous selection
      coefficients. Can be assigned a value equal to zero to eliminate synonymous
      selection.}
    \item{nsynVar}{A non-negative value that corresponds to the variance
      of the coefficients of the non-synonymous codons. That is, the
      the between-amino-acid variance.}
  % }
}
\details{
  Twenty random observations from a \code{Gamma(1,nsynVar^(-0.5))} distribution
  are sampled as the amino acid selection coefficients. The associated variance
  of the synonymous coefficients (\code{synVar}) is calculated as
  \code{nsynVar/vNvS}. If \code{nsynVar < 10^(-12)}, all the amino acid
  coefficients will be equal to a single random draw from a
  \code{Gamma(1,synVar^(-1))} distribution.
}
\value{
  Returns an object of class \code{aminoSC}, a vector that at least contains
  the following components.
  \describe{
    \item{\code{aacoeff}}{A vector of 22 elements where the first 20 elements
      represent the sampled amino acid coefficients. The 21st element
      is the associated variance of synonymous coefficients while the
      22nd element is the input \code{nsynVar} value. The coefficients are
      ordered in terms of the 1-letter amino acid IUPAC labels. That
      is, (A, C, . . . , W, Y).}
  }
}

\author{
  Hassan Sadiq
}

\seealso{
  An alternative approach for generating amino acid selection coefficients from
  a normal distribution \code{\link{aaGauss}}. There is \code{\link{codonCoeffs}}
  as well, a function designed to convert amino acid to codon selection
  coefficients.
}

\examples{
  aaGamma(0.50, 1e-04)
  
  aaGamma(1e-02, 0)
}
