# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Simulate Ornstein-Uhlenbeck Trend
ouEvolve <- function(xInit, deltaT, sysTheta, asymptoteVar, asymptoteMew){
    sys_mean <- asymptoteMew + (xInit-asymptoteMew) * exp((-1)*sysTheta*deltaT)
    sys_var <- asymptoteVar * (1 - exp((-2)*sysTheta*deltaT))
    xnew <- rnorm(1, sys_mean, sqrt(sys_var))
    return(xnew)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #