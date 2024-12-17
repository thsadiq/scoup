# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Deterministic Landscape Coefficient Update
aaSCupdate.discrete <- function(parameters){
    methd <- parameters@sampler
    p02xnodes <- parameters@lscape
    nodeIndex <- parameters@nodeIndex
    newParams <- p02xnodes[,nodeIndex]
    if(methd == 1){ scFunc <- aaGauss }else{ scFunc <- aaGamma }
    newSC <- scFunc(newParams["vNvS"], newParams["nsynVar"])
    return( newSC )
}

# ><>< # Deterministic Landscape Independent Site Simulation
sitesim.discrete <- function(parameters, nodeLength){
    popSize <- parameters@psize
    h85mean <- parameters@mrate
    h85kappa <- parameters@kappa
    s01x22 <- aaSCupdate(parameters)
    parentCodon <- initSeq( s01x22 )
    n_nodes <- round(ncol(parameters@lscape), 0)
    seqVector <- c(parentCodon)
    dndsVec <- c()
    
    for(a4 in seq(1,n_nodes)){
        new_seq_vec <- rep(NA, 2^a4)
        parameters@nodeIndex <- a4
        f01x22 <- aaSCupdate(parameters)
        s01x22 <- f01x22
        new_csc_vec <- codonCoeffs(s01x22)
        qmatrix <- subsMatrix(new_csc_vec, popSize, h85kappa, h85mean)
        dndsVec[a4] <- dndsCalculator( codonFreq(new_csc_vec),
                                    qmatrix, h85kappa, h85mean)
        
        for(a6 in seq(1,length(seqVector))){
            offspringID <- c(a6*2-1, a6*2)
            temp_codon <- c(NA, NA)
            for(a0 in seq(1,2)){
                init_codon <- seqVector[a6]
                new_codon <- branchSimulate(init_codon, nodeLength, qmatrix)
                init_codon <- new_codon["codon"]
                names(init_codon) <- NULL
                temp_codon[a0] <- init_codon
            }
            new_seq_vec[offspringID] <- temp_codon
        }
        seqVector <- new_seq_vec
    }
    return( list(codons=seqVector,dnds=dndsVec,initCodon=parentCodon) )
}

# ><>< # Deterministic Landscape Codon Sequence Alignment Simulation
alignsim.discrete <- function(adaptIn, seqIn, filename=NA){
    theMtrx <- adaptIn@lscape
    nNodes <- round(ncol(theMtrx), 0)
    nsite <- seqIn@sites
    ntaxons <- 2^nNodes
    nodeLen <- seqIn@branchL
    treeData <- biTree(ntaxons, nodeLen, adaptIn@t3mdl)
    dnds_matrix <- array(NA, c(nNodes, nsite))
    alignment <- array(NA, c(ntaxons, nsite))
    xMd <- theMtrx["vNvS",ncol(theMtrx)] > 0
    
    for(b0 in seq(1,nsite)){
        simOut <- sitesim(adaptIn, nodeLen)
        alignment[,b0] <- simOut$codons
        dnds_matrix[,b0] <- simOut$dnds
    }
    aWord <- ifelse(adaptIn@sampler==1, "Gauss", "Gamma")
    mString <- paste( apply(theMtrx, 2,
        function(a) paste(a, collapse=",")), collapse="|")
    notes <- paste0("Discrete Parameters: method=", aWord,
        "Model Parameters: Pop.Size=", adaptIn@psize,", HKY85.Kappa=",
        adaptIn@kappa, ", HKY85.Mu=", adaptIn@mrate,
        " node-wise=(vNvS,nsynVar|", mString, ")")
    commentText <- paste0(seqIn@details, ";   ", notes)
    cdnSEQs <- seqDframe(alignment)
    if(is.null(filename)){
        coloredSQS <- seqColored(cdnSEQs)}else{coloredSQS <- NA}
    if(is(filename, "character")){
        empty <- seqWriter(alignment, treeData, commentText, filename)
    }
    generatedSEQ <- new("scoup", seqs=alignment, cseq=cdnSEQs,
        DNDS=dnds_matrix, aInfo=commentText, seqCOL=coloredSQS)
    return(generatedSEQ)
}

# ><>< # "discrete" Object Methods
setMethod("effpop", "discrete", function(x) x@psize)
setMethod("lscape", "discrete", function(x) x@lscape)
setMethod("sampler", "discrete", function(x) detectApp(x))
setMethod("kappa", signature("discrete"), function(x) x@kappa)
setMethod("hky85mu", signature("discrete"), function(x) x@mrate)
setMethod("aaSCupdate", signature("discrete"), aaSCupdate.discrete)
setMethod("show", "discrete", function(object) {
    cat("\n", is(object)[[1]])
    cat("\nDeterministic simulation model:\n")
    cat("\tPopulation size = ", object@psize,"\n")
    cat("\tNumber of terminal taxa = ",2^ncol(object@lscape),"\n")
})
setMethod("alignsim",
    signature(adaptIn="discrete", seqIn="seqParameters"), alignsim.discrete)
setMethod("sitesim",
    signature(parameters="discrete", nodeLength="numeric"), sitesim.discrete)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #