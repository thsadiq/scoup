% File src/library/scoup/man/..Rd
\name{codonMutation}
\alias{codonMutation}
\alias{kappa}
\alias{hky85mu}

\title{Construct a Codon Mutation Matrix}
\description{
    Build an instantaneous codon mutation matrix from a HKY85 DNA mutation
    matrix.
}
\usage{
    codonMutation(kappa=4, mrate=0.25)
}
\arguments{
    \item{kappa}{Transition to transversion rate between nucleotide pairs.
        Default value is \code{4}.}
    \item{mrate}{Average rate at which mutation occurs between DNA pairs.
        Default value is \code{0.25}.}
}
\details{
    The mutation rate is set to zero for pairs of codons that differ by more
    than one nucleotide position or for two codons that do not differ. For a
    pair of codons that differ at one nucleotide position, the mutation rate
    is set to the corresponding mutation rate of the non-matching nucleotide
    pair. The nucleotide rate is obtained as proposed by Hasegawa, Kishino
    and Yano (1985).
}
\value{
    \describe{
        \item{\code{codonMutate}}{Codon mutation matrix such that the rows and
            the columns are arranged with respect to the IUPAC naming structure
            of nucleotide triplets in alphabetical order.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    Codon substution matrix generating function \code{\link{subsMatrix}} and
    the fixation matrix generating function \code{\link{fixMatrix}}.
}
\references{
    Hasegawa, M., Kishino, H. and Yano, T. (1985). Dating of the Human-Ape
    Splitting by a Molecular Clock of Mitochondria DNA, \emph{Journal of
    Molecular Evolution} \bold{22}: 160-174.
}
\examples{
    cm0 <- codonMutation()
    head(cm0)
    dim(cm0)
    cm1 <- codonMutation(1, 0.032)
    tail(cm1)
    dim(cm1)
}
