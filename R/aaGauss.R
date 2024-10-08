# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Amino Acid Selection Coefficients (reviewed)
aaGauss <- function(vNvS, nsynVar){
    if(is(nsynVar, "character")) stop("`nsynVar` should be non-negative!")
    vaRatio <- ifelse(vNvS <= 0, Inf, vNvS)
    if(nsynVar > 1e-12){
        synVar <- (1 / vaRatio) * nsynVar
        aacoefs <- abs( rnorm(20, 0, sqrt(nsynVar)))
    }else{
        synSD <- ifelse(vNvS>0, sqrt(vaRatio), 0)
        meansyn <- abs( rnorm(1, 0, synSD))
        aacoefs <- rep(meansyn, 20)
        synVar <- synSD^2
    }
    names(synVar) <- NULL
    names(nsynVar) <- NULL
    aacoeff <- new("aminoSC",
        coeffs=abs(aacoefs), synVar=synVar, nsynVar=nsynVar)
    return(aacoeff)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #