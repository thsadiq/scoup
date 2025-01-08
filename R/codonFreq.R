# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Codon Frequencies
codonFreq <- function(sc01x61){
    smax <- max( coeffs(sc01x61))
    sShifted <- coeffs(sc01x61) - smax
    sweights <- exp( sShifted )
    codonF <- sweights / sum(sweights)
    codonF <- new("codonvalues", cdnums=codonF)
    return(codonF)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #