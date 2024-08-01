% File src/library/scoup/man/..Rd
\name{subsMatrix}
\alias{subsMatrix}

\title{Build Mutation-Selection Codon Substitution Matrix}
\description{
  Construct an instantaneous codon substitution matrix based on the
  mutation-selection framework.
}
\usage{
  subsMatrix(sc01x61, effpopsize)
}
\arguments{
  % \describe{
    \item{sc01x61}{Vector of selection coefficients associated with the
      61 sense codons, ordered alphabetically according to the nucleotide
      triplets and the IUPAC naming structure.}
    \item{effpopsize}{Effective population size.}
  % }
}
\details{
  Given an arbitrary scaling constant (\code{k}), codon fixation rates
  (\eqn{f_{ij}}) and mutation rates (\eqn{m_{ij}}), the instantaneous
  rate by which codon \code{i} is substituted by another codon \code{j}
  may be expressed as follows.
  \deqn{
    q_{ij} =
    \begin{cases}
      k \times \text{m}_{ij}^{} \times \text{f}_{ij}^{}
        & \text{if $i$ and $j$ differs by only one nucleotide,} \\[1ex]
      0 & \text{if $i$ and $j$ differs by more than one nucleotide,}
    \end{cases}
  }
  and \eqn{q_{ii}^{}=-\sum_{j}q_{ij}^{}}. The \code{HKY85} nucleotide
  mutation matrix was embedded into the \code{alignsim} function (with
  transition to transversion rate, \eqn{\kappa} = 4, average rate,
  \eqn{\mu} = 1 and equal equilibrium frequencies).
}
\value{
  \describe{
    \item{\code{mainMatrix}}{Instantaneous codon substitution matrix such that
      the rows and the columns are arranged with respect to the IUPAC naming
      structure of nucleotide triplets in alphabetical order.}
  }
}
\author{
  Hassan Sadiq
}
\seealso{
  Selection coefficients conversion function \code{\link{codonCoeffs}} and
  fixation matrix generating function \code{\link{fixMatrix}}.
}
\references{
  Halpern, A. L. and Bruno, W. J. (1998). Evolutionary Distances for
  Protein-Coding Sequences: Modelling Site-Specific Residue Frequencies,
  \emph{Molecular Biology and Evolution} \bold{15}(7): 910-917.
  
  Hasegawa, M., Kishino, H. and Yano, T. (1985). Dating of the Human-Ape
  Splitting by a Molecular Clock of Mitochondria DNA, \emph{Journal of
  Molecular Evolution} \bold{22}: 160-174.
}
\examples{
  aacoeffs <- aaGauss(1e-03, 0)
  codonsc <- codonCoeffs(aacoeffs)
  subsMatrix(codonsc, 1000)
}
