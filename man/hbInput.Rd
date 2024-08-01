% File src/library/scoup/man/..Rd
\name{hbInput}
\alias{hbInput}

\title{Generate Halpern-Bruno Substitution Model Parameters}
\description{
  Create a \code{hbParameters} object that contains values necessary to construct
  a Halpern-Bruno mutation-selection codon substitution model.
}
\usage{
  hbInput(hbVector=0)
}
\arguments{
  % \describe{
    \item{hbVector}{A named vector that provides values of the parameters
      necessary to successfully create a codon substitution model for simulation
      purpose. See \bold{Details} for further information.}
  % }
}
\details{
  When fully specified, \code{hbVector} will be a four-element named vector. That
  is, \code{hbVector=c(Ne, meth, vNvS, nsynVar)}. \code{Ne} is an integer that
  represents the effective population size (default = 1000). \code{meth} is a
  binary integer that indicates the probability distribution from where the
  initial selection coefficients must be sampled. It can be \code{1} for a
  truncated Gaussian or \code{2} for a Gamma (default) distribution. \code{vNvS}
  is the ratio of the variance of the non-synonymous to synonymous selection
  coefficients. Default value is \code{1}. The user can set the variance of the
  non-synonymous selection coefficients with \code{nsynVar} (default is
  \eqn{10^{-5}}). This function was not intended for independent use. Rather
  as a complement to the \code{alignsim} simulation function.
}
\value{
  A \code{hbParameters} object that contains the following.
  \describe{
    \item{\code{psize}}{Effective population size.}
    \item{\code{vNvS}}{Ratio of the variance of the non-synonymous to the
        synonymous selection coefficients.}
    \item{\code{nsVar}}{Variance of the non-synonymous selection coefficients.}
    \item{\code{meth}}{Probability distribution for generating the initial
        selection coefficients.}
    \item{\code{words}}{Comments that contain details of the set parameters, for
        the output alignment file.}
  }
}
\author{
  Hassan Sadiq
}
\seealso{
  Selection coefficients sampling functions, \code{\link{aaGamma}} and
  \code{\link{aaGauss}} as well as the primary simulation function
  \code{\link{alignsim}}.
}
\examples{
  hbInput()
  hbInput(c(Ne=100, meth=2, vNvS=1e-08, nsynVar=1e-08))
}
