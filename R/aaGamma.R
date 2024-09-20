# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Amino Acid Selection Coefficients (reviewed)
aaGamma <- function(vNvS, nsynVar){
    if(is(nsynVar, "character")) stop("`nsynVar` should be non-negative!")
    vNvS <- ifelse(vNvS <= 0, Inf, vNvS)
    if(nsynVar > 1e-12){
        synVar <- (1 / vNvS) * nsynVar
        aacoefs <- rgamma(20, 1, sqrt(1/nsynVar))
    }else{
        synVar <- max((1 / vNvS), 0)  # Eliminate negative input
        meansyn <- rgamma(1, 1, sqrt(1/synVar))
        aacoefs <- rep(meansyn, 20)
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