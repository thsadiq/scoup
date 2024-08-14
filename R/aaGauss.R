# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                      <- aaGauss -> Function                      ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
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
    names(synVar) <- NULL; names(nsynVar) <- NULL
    aacoeff <- c(aacoefs, synVar=synVar, nsynVar=nsynVar)
    class(aacoeff) <- "aminoSC"
    return(aacoeff)
}
## ><>< ## Example ## ><>< ##
# anEG <- aaGauss(1, 1e-04)
# print(anEG)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #