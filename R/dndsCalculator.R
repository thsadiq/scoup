# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  <- dndsCalculator -> Function.                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Obtain Analytical dN/dS Estimate
dndsCalculator <- function(pi01x61, q61x61){
    piMatrix <- matrix(pi01x61, 61) %*% matrix(rep(1,61), 1)
    dnTop <- sum(piMatrix * q61x61 * nonsynonymID)
    dnBot <- sum(piMatrix * codon_m_matrix * nonsynonymID)
    dn <- dnTop / dnBot
    #####
    dsTop <- sum(piMatrix * q61x61 * synonymID)
    dsBot <- sum(piMatrix * codon_m_matrix * synonymID)
    ds <- dsTop / dsBot
    dnds <- dn / ds
    return(dnds)
}
## ><>< ## Example:
# aaEG2 <- aaGauss(0.5, 1e-03)
# codonsc <- codonCoeffs(aaEG2)
# piFreq <- codonFreq(codonsc)
# smat <- subsMatrix(codonsc, effpopsize=1000)
# DNS <- dndsCalculator(piFreq, smat)
# print(DNS)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #