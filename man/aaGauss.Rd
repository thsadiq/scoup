% File src/library/scoup/man/..Rd
\name{aaGauss}
\alias{aaGauss}

\title{Obtain Gaussian Distributed Amino Acid Selection Coefficients}
\description{
  Obtain a vector of values from a Normal distribution that may be conveniently
  used as amino acid selection coefficients.
}
\usage{
  aaGauss(vNvS, nsynVar)
}
\arguments{
  % \describe{
    \item{vNvS}{Ratio of the variance of the selection coefficients
      of the non-synonymous codons relative to the variance of the synonymous
      coefficients. It can be assigned a value equal to zero to eliminate
      synonymous selection.}
    \item{nsynVar}{A non-negative value that corresponds to the variance
      of the coefficients among the non-synonymous codons. That is, the
      between-amino-acid variance.}
  % }
}
\details{
  An observation is sampled from a \code{Normal(0,nsynVar)} distribution
  independently for each of the 20 amino acid residues. The variance of the
  synonymous selection coefficients (\code{synVar}) is calculated as
  \code{nsynVar/vNvS}. If \code{nsynVar} is less than \eqn{10^(-12)}, all the
  amino acid coefficients will be equal to a single random draw from a
  \code{Normal(0,synVar)} distribution.
}
\value{
  Returns an object of class \code{aminoSC}, a vector that at least contains
  the following component.
  \describe{
    \item{\code{aacoeff}}{A vector of 22 elements of which the first 20 elements
      represent the sampled amino acid coefficients. The 21st element is the
      variance of the synonymous coefficients and the 22nd element is the input
      \code{nsynVar} value. The coefficients are ordered in terms of the
      1-letter amino acid IUPAC labels. that is, (A, C, ..., W, Y).}
  }
}

\author{
  Hassan Sadiq
}

\seealso{
  An alternative sampling function, \code{\link{aaGamma}} is also available. The
  \code{\link{codonCoeffs}} function requires the output from this function
  (or from \code{\link{aaGamma}}).
}

\examples{
  aaGauss(0.50, 1e-04)
  
  aaGauss(1e-02, 0)
}
