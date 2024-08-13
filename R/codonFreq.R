# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                     <- codonFreq -> Function                     ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Codon Frequencies
codonFreq <- function(sc01x61){
    sweights <- exp(sc01x61)
    if(all(sc01x61 > 1)) sweights <- log(sc01x61)
    codonF <- sweights / sum(sweights)
    class(codonF) <- "codonvalues"
    return(codonF)
}
# ><>< ## Example:
# aaEG1 <- aaGamma(1e-03, 0)
# codonsc <- codonCoeffs(aaEG1, 4)
# piFreq <- codonFreq(codonsc)
# print(piFreq)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #