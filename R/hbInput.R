# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< ====== ><>< # Halpern-Bruno Model Input Generator
hbInput <- function(hbVector=0){
    stopifnot("`hbInput` arg isn't a numeric vector" = is(hbVector,"numeric"))
    runCheck <- warner(hbVector, c("Ne","meth","vNvS","nsynVar"))
    Ne <- hbVector["Ne"]
    popn <- ifelse(is.na(Ne), 1000, Ne)
    mthd <- hbVector["meth"]
    mRule <- is.na(mthd) | (mthd != 1)
    vaRatio <- hbVector["vNvS"]
    vNvS <- ifelse(is.na(vaRatio), 1, vaRatio)
    nsVar <- hbVector["nsynVar"]
    nonVar <- ifelse(is.na(nsVar), 1e-05, nsVar)
    if(mRule){
        methd <- 2
        fname <- "Gamma"
    }else{
        methd <- 1
        fname <- "Gauss"
    }
    cptext <- paste0(fname,"(vNvS=",vNvS,",nsynVar=",nonVar,")")
    mdtxt <- paste0("Model Parameters: population.size=",
                    popn,"  method=", cptext)
    names(popn) <- NULL
    names(vNvS) <- NULL
    names(nonVar) <- NULL
    names(methd) <- NULL
    modelEntry <- new("hbParameters", psize=popn, vNvS=vNvS,
                    nsynVar=nonVar, sampler=methd, words=mdtxt)
    return(modelEntry)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #