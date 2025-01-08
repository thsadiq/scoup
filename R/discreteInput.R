# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Deterministic Landscape Input
discreteInput <- function(defList=list()){
    stopifnot("`discreteInput` arg. must be a list!" = is(defList,"list"))
    everyTag <- c("p02xnodes","technique","psize","nodeIndex","leafModel")
    runCheck <- warner(defList, everyTag)
    aaParams <- defList$p02xnodes[c(1,2),]
    app <- defList$technique
    if(is.null(app)) app <- 2
    matched <- paste(app)
    match.arg(matched, c(1,2))
    nID <- defList$nodeIndex
    if(is.null(nID)) nID <- 1
    xMD <- defList$leafModel
    if(is.null(xMD)) xMD <- NA_character_
    pnum <- defList$pSize
    if(is.null(pnum)) pnum <- 1000
    if(is.null(aaParams)) aaParams <- rbind(rep(1,4),rep(1e-05,4))
    rownames(aaParams) <- c("vNvS","nsynVar")
    names(aaParams) <- NULL
    names(app) <- NULL
    names(nID) <- NULL
    names(pnum) <- NULL
    names(xMD) <- NULL
    defEntry <- new("discrete", lscape=aaParams, sampler=app,
        nodeIndex=nID, psize=pnum, t3mdl=xMD)
    return(defEntry)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #