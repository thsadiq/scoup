# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Transform Amino Acid to Codon Selection Coefficients (Reviewed)
codonCoeffs <- function(s01x22){ # deleted "fixed=NULL"
    logicMat <- matrix(FALSE, 20, 61)
    codonSC <- rep(0, 61)
    count <- 0
    for(a0 in seq(1,20)){
        newID <- which(amino2codon==a0)
        logicMat[a0,newID] <- TRUE
    }
    
    for(a1 in seq(1,20)){
        count <- count + 1
        syns <- logicMat[a1,]
        if(sum(syns) == 1){
            codonSC[syns] <- coeffs(s01x22)[count]
        }else{
            minU <- coeffs(s01x22)[count] - sqrt(3*synVar(s01x22))
            minU <- max(minU, 0)
            maxU <- coeffs(s01x22)[count] + sqrt(3*synVar(s01x22))
            codonSC[syns] <- runif(sum(syns), minU, maxU) } }
    codonSC <- new("codonvalues", cdnums=codonSC)
    return(codonSC)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #