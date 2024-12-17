# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Frequency-Dependent Landscape Coefficient Update
wUpdate <- function(omegaInput){
    aaFull <- seq(1,20)
    dndsBank <- rep(NA, 20)
    editedQmat <- matrix(NA, 61, 61)
    aaZero <- aaFull[!(aaFull %in% omegaInput@aaPlus)]
    if(omegaInput@sampler==1){
        aaCreate <- aaGauss
    }else{
        aaCreate <- aaGamma
    }
    
    for(a3 in seq(1,20)){
        synonymsID <- which(amino2codon == a3)
        if(a3 %in% aaZero){ s01x22 <- aaCreate(0, 0) }else{
            s01x22 <- aaCreate(omegaInput@vNvS, omegaInput@nsynVar) }
        csc01x61 <- codonCoeffs(s01x22)
        subMTX <- subsMatrix(csc01x61, omegaInput@psize,
            omegaInput@kappa, omegaInput@mrate)
        dndsBank[a3] <- dndsCalculator( codonFreq(csc01x61), subMTX,
            omegaInput@kappa, omegaInput@mrate)
        editedQmat[synonymsID,] <- subMTX[synonymsID,]
    }
    cqOutput <- list(qFrame=editedQmat, dnds=dndsBank)
    return(cqOutput)
}

# ><>< # Frequency-Dependent Landscape Independent Site Simulation
sitesim.omega <- function(parameters, nodeLength, ntaxa){
    wOutput <- wUpdate( parameters )
    if(parameters@sampler==1){
        aaCreate <- aaGauss
    }else{
        aaCreate <- aaGamma
    }
    s01x22 <- aaCreate(parameters@vNvS, parameters@nsynVar)
    n_nodes <- ceiling( log(ntaxa, 2))
    parentCodon <- initSeq( s01x22 )
    seqVector <- c(parentCodon)
    
    for(a4 in seq(1,n_nodes)){
        new_seq_vec <- rep(NA, 2^a4)
        
        for(a6 in seq(1,length(seqVector))){
            offspringID <- c(a6*2-1, a6*2)
            temp_codon <- c(NA, NA)
            for(a0 in seq(1,2)){
                init_codon <- seqVector[a6]
                new_codon <- branchSimulate(init_codon,
                                            nodeLength, wOutput$qFrame)
                init_codon <- new_codon["codon"]
                names(init_codon) <- NULL
                temp_codon[a0] <- init_codon
            }
            new_seq_vec[offspringID] <- temp_codon
        }
        seqVector <- new_seq_vec
    }
    return( list(codons=seqVector, dnds=wOutput$dnds, initCodon=parentCodon) )
}

# ><>< # Frequency-Dependent Landscape Codon Sequence Alignment Simulation
alignsim.omega <- function(adaptIn, seqIn, filename=NA){
    alignment <- array(NA, c(seqIn@taxa, seqIn@sites))
    dnds_matrix <- array(NA, c(20, seqIn@sites))
    
    for(b0 in seq(1,seqIn@sites)){
        simOut <- sitesim(adaptIn, seqIn@branchL, seqIn@taxa)
        alignment[,b0] <- simOut$codons
        dnds_matrix[,b0] <- simOut$dnds
    }
    aWord <- ifelse(adaptIn@sampler==1, "Gauss", "Gamma")
    nzaa <- paste0("(", paste0(adaptIn@aaPlus,collapse=","), ")")
    notes <- paste0("Omega-Based Parameters: method=", aWord, " vNvS=",
        adaptIn@vNvS, " nsynVar=", adaptIn@nsynVar, " NonZero-AA=", nzaa,
        " hky85mu=", adaptIn@mrate, " hky85kappa=", adaptIn@kappa)
    commentText <- paste0(seqIn@details, ";   ", notes)
    cdnSEQs <- seqDframe(alignment)
    if(is.null(filename)){
        coloredSQS <- seqColored(cdnSEQs)}else{coloredSQS <- NA}
    if(is(filename, "character")){
        empty <- seqWriter(alignment, seqIn@phylogeny, commentText, filename)
    }
    omgOUT <- new("scoup", seqs=alignment, DNDS=dnds_matrix,
        aInfo=commentText, cseq=cdnSEQs, seqCOL=coloredSQS)
    return(omgOUT)
}

# ><>< # "omega" Object Methods
setMethod("lscape", signature("omega"), function(x){
    message("\nAmino acid(s) with non-zero vN/vS:")
    cmnt <- paste(aacids[x@aaPlus])
    message(cmnt, "\n\n")
})
setMethod("vNvS", signature("omega"), function(x) x@vNvS)
setMethod("kappa", signature("omega"), function(x) x@kappa)
setMethod("effpop", signature("omega"), function(x) x@psize)
setMethod("hky85mu", signature("omega"), function(x) x@mrate)
setMethod("nsynVar", signature("omega"), function(x) x@nsynVar)
setMethod("sampler", signature("omega"), function(x) detectApp(x))
setMethod("alignsim",
    signature(adaptIn="omega",seqIn="seqParameters"), alignsim.omega)
setMethod("sitesim",
    signature(parameters="omega",nodeLength="numeric"), sitesim.omega)
setMethod("show", signature("omega"), function(object){
    cat("\n", is(object)[[1]])
    cat("\nFrequency-dependent simulation model:\n")
    cat("\tvNvS = ", object@vNvS,"\n")
    cat("\tNumber of AA with non-zero vNvS = ",length(object@aaPlus),"\n\n")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #