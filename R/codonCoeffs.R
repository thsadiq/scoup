# ><>< ========================================================================= ><>< #
# ><><   scoup: Simulate Codon Sequences with Darwinian Selection Incorporated   ><>< #
# ><><                     as an Ornstein-Uhlenbeck Process.                     ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                             General Functions                             ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                               18 May, 2024.                               ><>< #
# ><>< ========================================================================= ><>< #

# ><>< # Transform Amino Acid to Codon Selection Coefficients (Reviewed)
codonCoeffs <- function(s01x22, fixed=NULL){
  logicMat <- matrix(F, 20, 61); codonSC <- rep(0, 61); count <- 0
  for(a0 in seq(1,20)){ newID <- which(amino2codon==a0); logicMat[a0,newID] <- T }
  synVar <- s01x22["synVar"]; nsynVar <- s01x22["nsynVar"]; aacoeff <- s01x22[seq(1,20)]
  if(nsynVar > 1e-12){ aaID <- seq(1,20) }else{ aaID <- fixed }
  
  for(a1 in aaID){ count <- count + 1; syns <- logicMat[a1,]
    if(sum(syns) == 1){ codonSC[syns] <- aacoeff[count] }else{
      minU <- aacoeff[count] - sqrt(3*synVar); minU <- max(minU, 0)
      maxU <- aacoeff[count] + sqrt(3*synVar)
      codonSC[syns] <- runif(sum(syns), minU, maxU) } }
  class(codonSC) <- "codonvalues"
  return(codonSC)
}
## ><>< ============ ><>< #
## ><>< ## Example 1:
## ><>< ============ ><>< #
# aaEG0 <- aaGamma(1e-10, 1e-04)
# csc01 <- codonCoeffs(aaEG0)
# print(csc01)
## ><>< ============ ><>< #
## ><>< ## Example 2:
## ><>< ============ ><>< #
# aaEG0 <- aaGauss(1e-10, 1e-04)
# csc02 <- codonCoeffs(aaEG0)
# print(csc02)
## ><>< ============ ><>< #
## ><>< ## Example 3:
## ><>< ============ ><>< #
# aaEG1 <- aaGauss(1e-03, 0)
# csc03 <- codonCoeffs(aaEG1, c(2,6))
# print(csc03)

# ><>< ========================================================================= ><>< #
# ><><                              CODE ENDS HERE.                              ><>< #
# ><>< ========================================================================= ><>< #