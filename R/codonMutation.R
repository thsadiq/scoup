# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Codon-Based HKY85 Mutation Matrix
codonMutation <- function(kappa=4, mrate=0.25){
    codonMutate <- matrix(0, 61, 61)
    hky85 <- mrate * matrix(c(0, 1, kappa, 1, 1, 0,
        1, kappa, kappa, 1, 0, 1, 1, kappa, 1, 0), 4)
    for(a3 in seq(1,61)){ for(a4 in seq(1,61)){
        if(!is.na(codonNuc[a3,a4])){
            diffpost <- eval( parse(text=codonNuc[a3,a4]))
            codonMutate[a3,a4] <- hky85[diffpost[1], diffpost[2]]
    } } }
    return(codonMutate)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #