% File src/library/scoup/man/..Rd
\name{ouEvolve}
\alias{ouEvolve}

\title{Simulate the Trend of an Ornstein-Uhlenbeck Process}
\description{
    Simulate the next state of an Ornstein-Uhlenbeck (OU) process for a given
    value.
}
\usage{
    ouEvolve(xInit, deltaT, sysTheta, asymptoteVar, asymptoteMew)
}
\arguments{
    \item{xInit}{Starting point of the OU process.}
    \item{deltaT}{Jump size.}
    \item{sysTheta}{Reversion rate.}
    \item{asymptoteVar}{Asymptotic variance.}
    \item{asymptoteMew}{Asymptotic mean.}
}
\details{
    The state at time \code{k} (that is, \eqn{x_{t_{k}}}) of a process that
    evolves according to an OU algorithm can be expressed as an observation
    from a Gaussian distribution as follows.
    \deqn{
        x_{t_{k}^{}}^{} \sim \textsf{Normal}\bigg(\mu + \left(x_{t_{k-1}^{}}^{}-
        \mu\right)e^{-\theta \Delta t}_{};\,\frac{\sigma^{2}_{}}{2\theta}\left(1
        -e^{-2\theta \Delta t}_{}\right)\bigg)
    }
    Observe that when time interval (\code{deltaT})
    \eqn{\Delta t = t_{k}-t_{k-1}} approaches infinity, the asymptotic mean
    (\code{asymptoteMew}) and the asymptotic variance (\code{asymptoteVar})
    of the distribution are \eqn{\mu} and \eqn{\Sigma^2=\sigma^2/2\theta}
    respectively, where \eqn{\theta} is the reversion rate. The function
    may be used to simulate toy OU trends. See \bold{Example} below.
}
\value{
    \describe{
        \item{\code{xnew}}{A scalar that represents the updated state of the OU
            process.}
    }
}
\author{
    Hassan Sadiq
}
\references{
    Uhlenbeck, G. E. and Ornstein, L. S. (1930), On the Theory of the Brownian
    Motion, \emph{Physical Review} \bold{36}: 823-841.
}
\examples{
    x0 <- 0.015
    xvec <- c()
    xvec[1] <- x0
    for(k in seq(2,100)){
        xstate <- ouEvolve(x0, 0.002, 10, 0.001, 0)
        xvec[k] <- xstate
        x0 <- xstate
    }
    plot(xvec, type="l")
}
