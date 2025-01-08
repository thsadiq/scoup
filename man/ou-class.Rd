\name{ou-class}
\docType{class}
\alias{ou-class}
\alias{aaSCupdate,ou-method}
\alias{alignsim,ou,seqParameters-method}
\alias{asymMean}
\alias{asymMean,ou-method}
\alias{asymVar}
\alias{asymVar,ou-method}
\alias{reversion}
\alias{reversion,ou-method}
\alias{show,ou-method}
\alias{sitesim,ou,numeric-method}

\title{Ornstein-Uhlenbeck Stochastic Simulation Model Object}
\description{
    Contains the inputs that are necessary to implement the Ornstein-Uhlenbeck
    (OU) evolutionary algorithm.
}

\section{Objects from the Class}{
    Class \code{ou} objects can be created by calls of the form \code{new("ou",
    var=??, theta=??, mu=??, words=??)}. This type of object is returned by
    the \code{\link{ouInput}} function in the \code{\link{scoup}} package. It
    is an important input of the \code{\link{alignsim}} function when
    interested in codon sequences that evolved following the OU framework.
}

\section{Slots}{
    \describe{
        \item{\code{var}:}{numeric value returned by the \code{asymVar}
            method.}
        \item{\code{theta}:}{numeric value returned by the \code{reversion}
            method.}
        \item{\code{mu}:}{numeric value returned by the \code{asymMean}
            method.}
        \item{\code{words}:}{descriptive text that contains details of the
            input parameter values. Useful as reference comments to be
            included in the generated sequence alignment.}
    }
}

\section{Methods}{
    \describe{
        \item{aaSCupdate}{\code{signature(parameters = "ou")}: background
            function that is not intended for end-use. It updates the amino
            acid selection coefficients intermittently during the sequence
            simulation process.}
        \item{alignsim}{\code{signature(adaptIn="ou", seqIn="seqParameters")}:
            primary simulation function available in the \code{\link{scoup}}
            package.}
        \item{asymMean}{\code{signature(x="ou")}: asymptotic mean, \eqn{\mu},
            of the OU evolutionary algorithm.}
        \item{asymVar}{\code{signature(x="ou")}: asymptotic variance,
            \eqn{\Sigma^2}, of the OU evolutionary framework.}
        \item{reversion}{\code{signature(x="ou")}: reversion parameter,
            \eqn{\theta}, that acts as a selective pull in the OU process.}
        \item{show}{\code{signature(object="ou")}: prints the input values of
            \eqn{\Sigma^2}, \eqn{\mu} and \eqn{\theta}.}
        \item{sitesim}{\code{signature(parameters="ou", nodeLength="numeric")}:
            background function that is not to be used by an end-user. It
            generates the DNA data at each site independently.}
    }
}

\author{
    Hassan Sadiq
}

\seealso{
    Further help may be sought by consulting the help pages of the main
    simulation function, \code{\link{alignsim}}, and/or the page of the
    object generating function, \code{\link{ouInput}}.
}

\examples{
    o1 <- ouInput( c(eVar=1e-02, Theta=10))
    asymMean(o1)
    asymVar(o1)
}

\keyword{classes}
