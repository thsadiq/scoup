# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Frequency-Dependent Simulation                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 23 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Frequency-Dependent Landscape Input
wInput <- function(wList=list()){
    stopifnot("`wInput` arg. must be a list!" = is(wList,"list"))
    compNames <- c("pSize","vNvS","aaPlus","technique","nsynVar")
    runCheck <- warner(wList, compNames)
    pnum <- wList$pSize
    if(is.null(pnum)) pnum <- 1000
    vaRatio <- wList$vNvS
    if(is.null(vaRatio)) vaRatio <- 1
    plusCOEF <- wList$aaPlus
    if(is.null(plusCOEF)) plusCOEF <- seq(1,20)
    approach <- wList$technique
    if(is.null(approach)) approach <- 2
    matched <- paste(approach)
    match.arg(matched, c(1,2))
    nsySigm2 <- wList$nsynVar
    if(is.null(nsySigm2)) nsySigm2 <- 1e-05
    names(pnum) <- NULL
    names(vaRatio) <- NULL
    names(plusCOEF) <- NULL
    names(approach) <- NULL
    names(nsySigm2) <- NULL
    wEntry <- new("omega", nsynVar=nsySigm2, psize=pnum,
        sampler=approach, aaPlus=plusCOEF, vNvS=vaRatio)
    return(wEntry)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #