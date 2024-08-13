# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Ornstein-Uhlenbeck Simulation.                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 29 Apr, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< ====== ><>< # Halpern-Bruno Model Input Generator
hbInput <- function(hbVector=0){
    Ne <- hbVector["Ne"]; popn <- ifelse(is.na(Ne), 1000, Ne)
    mthd <- hbVector["meth"]; mRule <- is.na(mthd) | (mthd != 1)
    vaRatio <- hbVector["vNvS"]; vNvS <- ifelse(is.na(vaRatio), 1, vaRatio)
    nsVar <- hbVector["nsynVar"]; nonVar <- ifelse(is.na(nsVar), 1e-05, nsVar)
    if(mRule){methd <- 2; fname <- "Gamma"}else{methd <- 1; fname <- "Gauss"}
    cptext <- paste0(fname,"(vNvS=",vNvS,",nsynVar=",nonVar,")")
    mdtxt <- paste0("Model Parameters: population.size=",
                    popn,"  method=", cptext)
    names(popn) <- NULL; names(vNvS) <- NULL
    names(nonVar) <- NULL; names(methd) <- NULL
    modelEntry <- list( psize=popn, vNvS=vNvS,
                        nsynVar=nonVar, meth=methd, words=mdtxt)
    class(modelEntry) <- "hbParameters"
    return(modelEntry)
}
# ><>< ====== ><>< # Example 1:
# print( hbInput() )
# ><>< ====== ><>< # Example 2:
# invec <- c(meth=2, Ne=1e+04, vNvS=1e-08, nsynVar=1e-05)
# exvec <- hbInput(invec)
# print( exvec )

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #