% File src/library/scoup/man/..Rd
\name{codonCoeffs}
\alias{codonCoeffs}

\title{Transform Amino Acid to Codon Selection Coefficients}
\description{
    Convert a 20-element vector of amino acid selection coefficients to a
    61-element vector of codon selection coefficients.
}
\usage{
    codonCoeffs(s01x22)
}
\arguments{
    \item{s01x22}{A 22-element vector of class \code{aminoSC}.}
}
\details{
    Consider a vector of amino acid selection coefficients, (\code{s_x:}
    \code{s_A}, \code{s_C}, \code{s_D}, ..., \code{s_W}, \code{s_Y}) that
    are subset of \code{s01x22}. All the synonymous codons that encode each
    amino acid are assigned independently sampled values from
    \code{Uniform(s_x - 3*synVar; s_x + 3*synVar)} distribution, where
    \code{synVar} is the synonymous variance and it is also retrieved from
    \code{s01x22}. For amino acids \code{M} and \code{W}, the corresponding
    codon coefficient is simply set equal to \code{s_M} and \code{s_W},
    respectively. The output from the function is of the order (\code{s_(AAA)},
    \code{s_(AAC)}, \code{s_(AAG)}, ..., \code{s_(TTC)}, \code{s_(TTG)},
    \code{s_(TTT)}), excluding the stop codons.
}
\value{
    Returns a \code{codonvalues} object that will contain the following.
    \describe{
        \item{\code{coeffs}}{A 61-element vector of codon selection
            coefficients ordered alphabetically with respect to the IUPAC
            nucleotide triplets nomenclature.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    \code{\link{aaGamma}} and \code{\link{aaGauss}}, functions useful for
    generating \code{aminoSC} objects.
}
\examples{
    # Example 1:
    aasc1 <- aaGamma(1e-10, 1e-04)
    ccfs0 <- codonCoeffs(aasc1)
    coeffs(ccfs0)
    
    # Example 2:
    aasc2 <- aaGauss(1e-10, 1e-04)
    ccfs1 <- codonCoeffs(aasc2)
    coeffs(ccfs1)

    # Example 3:
    aasc3 <- aaGauss(1e-03, 0)
    ccfs2 <- codonCoeffs(aasc3)
    coeffs(ccfs2)
}
