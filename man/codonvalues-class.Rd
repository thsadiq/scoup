\name{codonvalues-class}
\docType{class}
\alias{codonvalues-class}
\alias{coeffs,codonvalues-method}
\alias{freqs}
\alias{freqs,codonvalues-method}
\alias{show,codonvalues-method}

\title{Codon Frequencies and Selection Coefficients}

\description{
    A numerical vector of values that correspond to the selection coefficients
    of the sense codons.
}

\section{Objects from the Class}{
    Objects of this class (\code{codonvalues}) can be created by calls of the
    form \code{new("codonvalues", cdnums=...)}. Two codon-related
    transformation functions (that is, \code{\link{codonCoeffs}} and
    \code{\link{codonFreq}}) that are available in the \code{\link{scoup}}
    package return objects of this class.
}

\section{Slots}{
    \describe{
        \item{\code{cdnums}:}{vector of 61 values that could correspond to the
            selection coefficients or the frequencies of the sense codons
            depending on the method called.}
    }
}

\section{Methods}{
    \describe{
        \item{coeffs}{\code{signature(x = "codonvalues")}: vector of 61 values
            that correspond to the selection coefficients of the sense codons.
            The entries are ordered in increasing alphabetical order in terms
            of the IUPAC nucleotide triplets naming structure.}
        \item{freqs}{\code{signature(x = "codonvalues")}: vector of 61 values
            that correspond to the frequencies of the sense codons. The entries
            are ordered in increasing alphabetical order in terms of the IUPAC
            nucleotide triplets naming structure.}
        \item{show}{\code{signature(object = "codonvalues")}: prints the first
            six relevant (that is, coefficients or frequencies) codon related
            values.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    The pair of functions that return objects of this class include,
    the generating function for codon selection coefficients
    \code{\link{codonCoeffs}} and the \code{\link{codonFreq}} function that
    converts codon coefficients to frequencies.
}

\examples{
    aasc1 <- aaGamma(1e-10, 1e-04)
    ccfs0 <- codonCoeffs(aasc1)
    cFq <- codonFreq(ccfs0)
    coeffs(ccfs0)
    freqs(cFq)
}


\keyword{classes}
