\name{aminoSC-class}
\docType{class}
\alias{aminoSC-class}
\alias{coeffs}
\alias{coeffs,aminoSC-method}
\alias{nsynVar,aminoSC-method}
\alias{show,aminoSC-method}
\alias{synVar}
\alias{synVar,aminoSC-method}

\title{Amino Acid Selection Coefficients}

\description{
    A numerical vector of values that are associated with the amino acid
    selection coefficients.
}

\section{Objects from the Class}{
    Objects of this class (\code{aminoSC}) can be created by calls of the form
    \code{new("aminoSC", coeffs=..., synVar=..., nsynVar=...)}. The two amino
    acid sampling functions (that is, \code{\link{aaGamma}} and
    \code{\link{aaGauss}}) that are available in the \code{\link{scoup}}
    package return objects of this class.
    
}
\section{Slots}{
    \describe{
        \item{\code{coeffs}:}{numeric vector returned by the \code{coeffs}
            method.}
        \item{\code{synVar}:}{numeric value returned by the \code{synVar}
            method.}
        \item{\code{nsynVar}:}{numeric value returned by the \code{nsynVar}
            method.}
    }
}
\section{Methods}{
    \describe{
        \item{coeffs}{\code{signature(x = "aminoSC")}: vector of twenty values
            that correspond to the amino acid selection coefficients. The
            entries are ordered in increasing alphabetical order in terms
            of the one-letter IUPAC naming structure.}
        \item{nsynVar}{\code{signature(x = "aminoSC")}: variance of the
            probability distribution where the returned amino acid selection
            coefficients were sampled.}
        \item{show}{\code{signature(object = "aminoSC")}: summary of the
            contents of the \code{aminoSC} object including a snippet of the
            sampled coefficients as well as the values of the \code{synVar}
            \eqn{(\sigma_s^2)} and the \code{nsynVar} \eqn{(\sigma_n^2)}
            parameters.}
        \item{synVar}{\code{signature(x = "aminoSC")}: variance of the uniform
            distribution where the synonymous selection coefficients should be
            sampled.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    \code{\link{aaGamma}}, \verb{  } \code{\link{aaGauss}}
}

\examples{
    aasc1 <- aaGamma(1e-10, 1e-04)
    coeffs(aasc1)
    show(aasc1)
}

\keyword{classes}
