% File src/library/scoup/man/..Rd
\name{dndsCalculator}
\alias{dndsCalculator}

\title{Estimate dN/dS Value Analytically}
\description{
    Obtain an analytical estimate for the ratio of non-synonymous to synonymous
    rate (dN/dS) of codon substitution.
}
\usage{
    dndsCalculator(pi01x61, q61x61, kappa, mrate)
}
\arguments{
    \item{pi01x61}{Vector of sense codon frequencies that are ordered
        alphabetically with respect to nucleotide triplets according to IUPAC
        nomenclature.}
    \item{q61x61}{A \eqn{61 \times 61} matrix of sense codon instantaneous
        substitution rates, where the rows and the columns are ordered in terms
        of IUPAC-lettered nucleotide triplets.}
    \item{kappa}{Transition to transversion rate between nucleotide pairs
        suitable for the HKY85 DNA mutation model.}
    \item{mrate}{Average rate at which mutation occurs between DNA pairs
        suitable for the HKY85 DNA mutation model.}
}
\details{
    The returned dN/dS estimate is obtained from the ratio of the following
    expressions.
    \deqn{
    dN = \frac{\sum_{j}^{}\sum_{i\neq j}^{} \pi_{i}^{}\,\text{A}_{ij}^{}
        \,I_{N}^{}}{\sum_{j}^{}\sum_{i\neq j}^{} \pi_{i}^{}\,\mu_{ij}^{}
        \,I_{N}^{}},
    \hspace*{1.00cm}
    dS = \frac{\sum_{j}^{}\sum_{i\neq j}^{} \pi_{i}^{}\,\text{A}_{ij}^{}
        \,I_{S}^{}} {\sum_{j}^{}\sum_{i\neq j}^{} \pi_{i}^{}\,\mu_{ij}^{}
        \,I_{S}^{}},
    }
    where \eqn{A} and \eqn{\pi} are the input codon frequency vector
    (\code{pi01x61}) and the instantaneous substitution rate matrix
    (\code{q61x61}), respectively. The notation \eqn{\mu} denotes codon
    mutation matrix (embedded as a function of the HKY85 nucleotide model)
    while \eqn{I_{S}^{}} and \eqn{I_{N}^{}} are Boolean matrices with
    ones at positions occupied by synonymous and non-synonymous codons,
    respectively.
}
\value{
    \describe{
        \item{\code{dnds}}{An estimate for the corresponding dN/dS}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    Codon frequencies generating function \code{\link{codonFreq}},
    instantaneous substitution rate matrix function \code{\link{subsMatrix}}
    and the nucleotide to codon mutation matrix converter
    \code{\link{codonMutation}}.
}
\references{
    Spielman, S. J. and Wilke, C. O. (2015). The Relationship between dN/dS and
    Scaled Selection Coefficients, \emph{Molecular Biology and Evolution}
    \bold{32}(4): 1097-1108.
}
\examples{
    aasc <- aaGauss(0.5, 1e-03)
    codonsc <- codonCoeffs(aasc)
    piFreq <- codonFreq(codonsc)
    smat <- subsMatrix(codonsc, 1000, 4, 0.25)
    dndsCalculator(piFreq, smat, 4, 0.25)
}
