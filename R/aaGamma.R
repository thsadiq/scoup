# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                      <- aaGamma -> Function                      ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Amino Acid Selection Coefficients (reviewed)
aaGamma <- function(vNvS, nsynVar){
    vNvS <- ifelse(vNvS <= 0, Inf, vNvS)
    if(nsynVar > 1e-12){
        synVar <- (1 / vNvS) * nsynVar
        aacoefs <- rgamma(20, 1, sqrt(1/nsynVar))
    }else{
        synVar <- max((1 / vNvS), 0)  # Eliminate negative input
        meansyn <- rgamma(1, 1, sqrt(1/synVar))
        aacoefs <- rep(meansyn, 20)
    }
    names(synVar) <- NULL; names(nsynVar) <- NULL
    aacoeff <- c(aacoefs, synVar=synVar, nsynVar=nsynVar)
    class(aacoeff) <- "aminoSC"
    return(aacoeff)
}
## ><>< ## Example ## ><>< ##
# agEG <- aaGamma(1, 2e-04)
# print(agEG)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #