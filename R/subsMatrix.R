# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Codon Substitution Matrix
subsMatrix <- function(sc01x61, effpopsize){
    fmatrix <- fixMatrix(sc01x61, effpopsize)
    qmatrix <- codon_m_matrix * fmatrix
    return(qmatrix)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #