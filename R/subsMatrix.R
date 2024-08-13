# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                    <- subsMatrix -> Function.                    ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Codon Substitution Matrix
subsMatrix <- function(sc01x61, effpopsize){
    fmatrix <- fixMatrix(sc01x61, effpopsize)
    qmatrix <- codon_m_matrix * fmatrix
    return(qmatrix)
}
## ><>< ## Example:
# aaEG <- aaGauss(.5, 1e-03)
# codonsc <- codonCoeffs(aaEG)
# smatrix <- subsMatrix(codonsc, 1000)
# print( range(smatrix) )

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #