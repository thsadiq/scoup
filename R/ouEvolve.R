# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Ornstein-Uhlenbeck Simulation.                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 29 Apr, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Simulate Ornstein-Uhlenbeck Trend
ouEvolve <- function(xInit, deltaT, sysTheta, asymptoteVar, asymptoteMew){
    sys_mean <- asymptoteMew + (xInit-asymptoteMew) * exp((-1)*sysTheta*deltaT)
    sys_var <- asymptoteVar * (1 - exp((-2)*sysTheta*deltaT))
    xnew <- rnorm(1, sys_mean, sqrt(sys_var))
    return(xnew)
}
## ><>< ## Example:
# x0 <- 0.015
# xvec <- c()
# xvec[1] <- x0
# for(k in 2:100){
#   xstate <- ouEvolve(x0, 0.002, 10, 0.001, 0)
#   xvec[k] <- xstate
#   x0 <- xstate
# }
# plot(xvec, type="l")

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #