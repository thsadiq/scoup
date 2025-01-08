% File src/library/scoup/man/..Rd
\name{fixMatrix}
\alias{fixMatrix}

\title{Construct Fixation Rate Matrix}
\description{
    Construct a 61 \eqn{\times} 61 matrix that contains the rate at which a
    mutant (column) codon gets fixed at a position previously occupied by a
    resident (row) codon for all the pairwise combinations of the 61 sense
    codons.
}
\usage{
    fixMatrix(sc01x61, effpopsize)
}
\arguments{
    \item{sc01x61}{A vector of the selection coefficients for the sense codons,
        ordered alphabetically in terms of the IUPAC nucleotide triplets
        nomenclature.}
    \item{effpopsize}{Effective population size.}
}
\details{
    If the additive fitness of a mutant codon (\emph{j}) relative to a
    wild-type codon (\emph{i}) is given by \eqn{s_{ij}^{} = (s_{j}^{}-
    s_{i}^{})}, then the rate at which codon \emph{j} gets fixed in a
    codon position occupied by codon \emph{i} can be expressed as follows.
    \deqn{
        f_{ij}^{} =
        \begin{cases}
            \frac{1 - e^{-2s_{ij}^{}}_{}}{1 - e^{-2N_{\texttt{e}}^{}
                s_{ij}^{}}_{}} & \text{if codons $i$ and $j$ differs only by
                one nucleotide,} \\[1ex]
            0 & \text{otherwise,}
        \end{cases}
    }
    where \eqn{s_{i}^{},\,s_{j}^{} \in} \code{sc01x61} and \eqn{N_{e}^{}}
    (\code{effpopsize}) is the effective population size. If
    \eqn{s_{i}^{}-s_{j}^{} =0}, then \eqn{f_{ij}^{}=1}.
}
\value{
    \describe{
        \item{\code{fixMtx}}{A \eqn{61 \times 61} matrix of the fixation rates
            of the sense codons. The rows and columns are ordered alphabetically
            in terms of the IUPAC nucleotide triplets nomenclature. That is,
            \code{AAA}, \code{AAC}, \code{AAG}, \ldots, \code{TTG}, \code{TTT}.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    You may be interested in also consulting the \code{\link{codonCoeffs}}
    function designed to generate selection coefficients for the sense codons
    or similar functions, \code{\link{aaGamma}} and \code{\link{aaGauss}},
    targeted at amino acid.
}
\references{
    Bershtein, S. and Serohijos, A. W. R. and Shakhnovich, E. I. (2017),
    Bridging the Physical Scales in Evolutionary Biology: From Protein
    Sequence Space to Fitness of Organisms and Populations, \emph{Current
    Opinion in Structural Biology} \bold{42}: 31-40.

    Halpern, A. L. and Bruno, W. J. (1998). Evolutionary Distances for
    Protein-Coding Sequences: Modelling Site-Specific Residue Frequencies,
    \emph{Molecular Biology and Evolution} \bold{15}(7): 910-917.

    McCandlish, D. M. (2011), Visualizing Fitness Landscapes, \emph{Evolution}
    \bold{65}(6): 1544-1558.
}
\examples{
    aasc <- aaGauss(0.5, 1e-03)
    codonsc <- codonCoeffs(aasc)
    fixMatrix(codonsc, 1000)
}
