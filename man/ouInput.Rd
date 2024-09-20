% File src/library/scoup/man/..Rd
\name{ouInput}
\alias{ouInput}

\title{Populate Parameters of the Ornstein-Uhlenbeck Algorithm}
\description{
    Create an \code{ou} object that will contain the parameters necessary
    to simulate a codon sequence alignment that evolves according to an
    Ornstein-Uhlenbeck (OU) process.
}
\usage{
    ouInput(ouVector=0)
}
\arguments{
    \item{ouVector}{A vector that contains carefully named elements.
        Each element represents a parameter in an OU model. See \bold{Details}
        for more information.}
}
\details{
    In its full form, \code{ouVector} is a three-element vector. Its contents
    each represents part of the parameters required to implement an OU process.
    The vector contents include, \code{eMean}, \code{eVar} and \code{Theta}.
    Input \code{eMean} is the asymptotic mean (\eqn{\mu}) and zero is its
    default value. \code{eVar} denotes the asymptotic variance (\eqn{\Sigma^2}).
    It has a 0.01 default value. \code{Theta} (default = 0.01) represents the
    reversion rate (\eqn{\theta}). This function was aimed as a complement to
    \code{alignsim}, not for use in isolation.
}
\value{
    An \code{ou} object that contains the following.
    \describe{
        \item{\code{asymMean}}{Asymptotic mean of the OU process.}
        \item{\code{asymVar}}{Asymptotic variance of the OU process.}
        \item{\code{reversion}}{Reversion rate of the OU process.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    The Ornstein-Uhlenbeck state generating function \code{\link{ouEvolve}} and
    the \code{\link{alignsim}} simulation function.
}
\references{
    Uhlenbeck, G. E. and Ornstein, L. S. (1930), On the Theory of the Brownian
    Motion, \emph{Physical Review} \bold{36}: 823-841.
}
\examples{
    o0 <- ouInput()
    reversion(o0)
    o0
    
    o1 <- ouInput( c(eVar=1e-02, Theta=10))
    asymMean(o1)
    asymVar(o1)
}
