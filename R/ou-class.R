# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Ornstein-Uhlenbeck Landscape Coefficient Update
aaSCupdate.ou <- function(parameters, oldSC, bLength){
    newSC <- oldSC@coeffs
    for(q0 in seq(1,20)){
        newSC[q0] <- ouEvolve(oldSC@coeffs[q0], bLength,
            parameters@theta, parameters@var, parameters@mu) }
    aacoeff <- new("aminoSC",
        coeffs=abs(newSC), synVar=oldSC@synVar, nsynVar=oldSC@nsynVar)
    return( aacoeff )
}

# ><>< # Ornstein-Uhlenbeck Landscape Independent Site Simulation
sitesim.ou <- function(parameters, nodeLength, popSize, ntaxa, s01x22){
    parentCodon <- initSeq( s01x22 )
    n_nodes <- ceiling( log(ntaxa, 2))
    seqVector <- c(parentCodon)
    dndsVec <- c()
    
    for(a4 in seq(1,n_nodes)){
        new_seq_vec <- rep(NA, 2^a4)
        f01x22 <- aaSCupdate(parameters, s01x22, nodeLength)
        s01x22 <- f01x22
        new_csc_vec <- codonCoeffs(s01x22)
        qmatrix <- subsMatrix(new_csc_vec, popSize)
        dndsVec[a4] <- dndsCalculator( codonFreq(new_csc_vec), qmatrix)
        
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

# ><>< # Ornstein-Uhlenbeck Landscape Codon Sequence Alignment Simulation
alignsim.ou <- function(adaptIn, seqIn, modelIn, filename=NA){
    if(modelIn@sampler == 1){ scFunc <- aaGauss }else{ scFunc <- aaGamma }
    dnds_matrix <- array(NA, c(seqIn@nodes, seqIn@sites))
    alignment <- array(NA, c(seqIn@taxa, seqIn@sites))
    
    for(b0 in seq(1,seqIn@sites)){
        aaSC <- scFunc(modelIn@vNvS, modelIn@nsynVar)
        simOut <- sitesim(adaptIn,seqIn@branchL,modelIn@psize,seqIn@taxa,aaSC)
        alignment[,b0] <- simOut$codons
        dnds_matrix[,b0] <- simOut$dnds
    }
    cdnSEQs <- seqDframe(alignment)
    if(is.null(filename)){
        coloredSQS <- seqColored(cdnSEQs)}else{coloredSQS <- NA}
    commentText <- paste0(adaptIn@words,
        ";   ",seqIn@details,";   ",modelIn@words)
    if(is(filename, "character")){
        empty <- seqWriter(alignment, seqIn@phylogeny, commentText, filename)
    }
    ou_scoup <- new("scoup", seqs=alignment, DNDS=dnds_matrix,
                aInfo=commentText, cseq=cdnSEQs, seqCOL=coloredSQS)
    return(ou_scoup)
}

# ><>< # "ou" Object Methods
setMethod("show", signature("ou"), function(object){
    cat("\n", is(object)[[1]])
    cat("\nOrnstein-Uhlenbeck-Based Genetic Sequence Simulation Settings:")
    cat("\n\tAsymptotic Mean =", object@mu)
    cat("\n\tAsymptotic Variance =", object@var)
    cat("\n\tReversion Parameter =", object@theta, "\n\n")
})
setMethod("asymVar", signature("ou"), function(x) x@var)
setMethod("asymMean", signature("ou"), function(x) x@mu)
setMethod("reversion", signature("ou"), function(x) x@theta)
setMethod("alignsim",
    signature(adaptIn="ou", seqIn="seqParameters"), alignsim.ou)
setMethod("sitesim",
    signature(parameters="ou", nodeLength="numeric"), sitesim.ou)
setMethod("aaSCupdate", signature(parameters="ou"), aaSCupdate.ou)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #