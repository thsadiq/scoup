% File src/library/scoup/man/..Rd
\name{wInput}
\alias{wInput}

\title{Populate Frequency-Dependent Simulation Model Parameters}
\description{
    Create an \code{omega} object. The utility is for defining the parameters
    that are necessary to simulate codon sequences that mimic the
    evolutionary process described by the frequency-dependent models.
}
\usage{
    wInput(wList=list())
}
\arguments{
    \item{wList}{A list that may contain up to five named entries. See
        \bold{Details} for further information.}
}
\details{
    In its full form, \code{wList} will contain five named elements. The
    elements include (a.) \code{pSize}: an integer that represents the effective
    population size (default = 1000). (b.) \code{vNvS}: a numerical value that
    corresponds to the ratio of the variance of the non-synonymous to the
    synonymous selection coefficients (default = 1). (c.) \code{aaPlus}: its
    default is a vector of integers between 1 and 20, inclusive. It gives the
    indices, if the one-letter IUPAC amino acid notations were ordered
    alphabetically, of the residues that should be assigned non-zero coefficient
    variances. It is useful for executing directional evolution. (d.)
    \code{technique}: it informs of the preferred probability distribution
    where the selection coefficients should be sampled. It could be set as 1
    for Gaussian or 2 for Gamma (default) distribution. (e.) \code{nsynVar}:
    variance of the non-synonymous selection coefficients. This is a
    complementary function to \code{alignsim}.
}
\value{
    An \code{omega} object that contains the following.
    \describe{
        \item{\code{nsynVar}}{Variance of the non-synonymous selection
            coefficients.}
        \item{\code{technique}}{Probability density function for sampling the
            amino acid selection coefficients.}
        \item{\code{aaPlus}}{Indices of the amino acids to be assigned non-zero
            coefficient variance values.}
        \item{\code{vNvS}}{Ratio of the variance of the non-synonymous to the
            synonymous selection coefficients.}
        \item{\code{psize}}{Effective population size.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    Sequence simulation function \code{\link{alignsim}} as well as selection
    coefficient conversion functions \code{\link{aaGamma}} and
    \code{\link{aaGauss}}.
}
\references{
    Goldman, N. and Yang, Z. (1994), A Codon-Based Model of Nucleotide
    Substitution for Protein-Coding DNA Sequences, \emph{Molecular Biology and
    Evolution} \bold{11}(5): 725-736.
    
    Muse, S. V. and Gaut, B. S. (1994), A A Likelihood Approach for Comparing
    Synonymous and Nonsynonymous Nucleotide Substitution Rates, with Application
    to the Chloroplast Genome, \emph{Molecular Biology and Evolution}
    \bold{11}(5): 715-724.
}
\examples{
    w0 <- wInput()
    vNvS(w0)
    w0

    w1 <- wInput(list(aaPlus=c(4,2,11), nsynVar=10))
    lscape(w1)
    w1
}
