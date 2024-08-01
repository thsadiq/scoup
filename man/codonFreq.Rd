% File src/library/scoup/man/..Rd
\name{codonFreq}
\alias{codonFreq}

\title{Generate Codon Frequencies From Selection Coefficients}
\description{
  Obtain codon frequencies from specified selection coefficients in a way
  that accounts for the magnitude of the coefficients in the real number
  line.
}
\usage{
  codonFreq(sc01x61)
}
\arguments{
  % \describe{
    \item{sc01x61}{Vector of sense codon selection coefficients that
      are ordered alphabetically in terms of the IUPAC nucleotide triplets
      nomenclature.}
  % }
}
\details{
  This conversion to frequencies accommodates the magnitude and signs of the
  selection coefficients because the frequency for the \emph{i}th codon is
  estimated as:
  \deqn{\pi_{i}^{} =
    \begin{cases}
      \log(s_{i}^{}) \big/ \sum_{j=1}^{61} \log(s_{j}^{})
        & \text{if } s_{j}^{} > 0 \;\; \forall_{j}, \\[0.5cm]
      e^{s_{i}^{}}_{} \big/ \sum_{j=1}^{61} e^{s_{j}^{}}_{}
        & \text{otherwise,} \\
    \end{cases}
  }
  where \eqn{s_{i}^{} \in \code{sc01x61} } is the selection coefficient of the
  \emph{i}th codon.
}
\value{
  Returns a \code{codonvalues} object that contains the following.
  \describe{
    \item{\code{codonF}}{A vector of 61 fractional values that sum to one and
      represent the frequencies of sense codons that are ordered alphabetically
      in terms of the IUPAC nucleotide triplets nomenclature.}
  }
}
\author{
  Hassan Sadiq
}
\seealso{
  \code{\link{codonCoeffs}}, a function that produces codon selection coefficients
    that may be used as an input.
}
\examples{
  aaEG1 <- aaGamma(1e-03, 0)
  csc01 <- codonCoeffs(aaEG1, 4)
  codonFreq(csc01)
}
