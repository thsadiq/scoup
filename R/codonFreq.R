# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Codon Frequencies
codonFreq <- function(sc01x61){
    sweights <- exp( coeffs(sc01x61) )
    if(all( coeffs(sc01x61) > 1)) sweights <- log( coeffs(sc01x61) )
    codonF <- sweights / sum(sweights)
    codonF <- new("codonvalues", cdnums=codonF)
    return(codonF)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #