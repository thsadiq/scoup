# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                     <- fixMatrix -> Function                     ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Generate Fixation Matrix
fixMatrix <- function(sc01x61, effpopsize){
    fixMtx <- matrix(0, 61, 61)
    for(a1 in seq(1,61)){
        for(a2 in seq(1,61)){
            if(qmatID[a1,a2] != 0){
                s_num <- (-2) * (sc01x61[a2] - sc01x61[a1])
                s_den <- (-2) * effpopsize * (sc01x61[a2] - sc01x61[a1])
                if(s_den == 0){ fixMtx[a1,a2] <- 1 }else{
                    numer <- 1 - exp(s_num)
                    denom <- 1 - exp(s_den)
                    fixMtx[a1,a2] <- ifelse(abs(denom)==Inf, 1, numer/denom)
    } } } }
    return(fixMtx)
}
## ><>< ## Example:
# aaEG <- aaGamma(.5, 1e-03)
# codonsc <- codonCoeffs(aaEG)
# fMat <- fixMatrix(codonsc, 1000)
# print(head(fMat))

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #