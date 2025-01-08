% File src/library/scoup/man/..Rd
\name{codonFreq}
\alias{codonFreq}

\title{Generate Codon Frequencies From Selection Coefficients}
\description{
    Obtain codon frequencies from specified selection coefficients in a way
    that accounts for the magnitude of the coefficients on the real number
    line.
}
\usage{
    codonFreq(sc01x61)
}
\arguments{
    \item{sc01x61}{Vector of selection coefficients for the sense codons. They
        are expected to be ordered alphabetically in terms of the IUPAC
        nucleotide triplets nomenclature.}
}
\details{
    The conversion formula is based on a softmax identity such that the
    magnitude and the signs of the selection coefficients are accommodated.
    The frequency for the \emph{i}th codon is estimated as:
    \deqn{\pi_{i}^{} =
            \frac{e^{s_{i}^{}-c}_{}}{\sum_{j=1}^{61} e^{s_{j}^{}-c}_{}},
    }
    where \eqn{s_{i}^{} \in \code{sc01x61} } is the selection coefficient of
    the \emph{i}th codon and \eqn{c=\max_{j} s_{j}^{}}.
}
\value{
    Returns a \code{codonvalues} object that contains the following.
    \describe{
        \item{\code{freqs}}{A vector of 61 fractional values that sum to one
            and represent the frequencies of sense codons. The output will be
            ordered alphabetically in terms of the IUPAC nucleotide triplets
            nomenclature.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    \code{\link{codonCoeffs}}, a function that produces codon selection
    coefficients that may be used as an input.
}
\examples{
    aaEG1 <- aaGamma(1e-03, 0)
    csc01 <- codonCoeffs(aaEG1)
    cFq <- codonFreq(csc01)
    freqs(cFq)
}
