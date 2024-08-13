# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Frequency-Dependent Simulation                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 23 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

wInput <- function(wList=list()){
    pnum <- wList$pSize; if(is.null(pnum)) pnum <- 1000
    vaRatio <- wList$vNvS; if(is.null(vaRatio)) vaRatio <- 1
    plusCOEF <- wList$aaPlus; if(is.null(plusCOEF)) plusCOEF <- seq(1,20)
    approach <- wList$technique; if(is.null(approach)) approach <- 2
    nsySigm2 <- wList$nsynVar; if(is.null(nsySigm2)) nsySigm2 <- 1e-05
    names(pnum) <- NULL; names(vaRatio) <- NULL; names(plusCOEF) <- NULL
    names(approach) <- NULL;  names(nsySigm2) <- NULL
    wEntry <- list(nsynVar=nsySigm2, technique=approach,
        aaPlus=plusCOEF, vNvS=vaRatio, psize=pnum)
    class(wEntry) <- "omega"
    return(wEntry)
}
# ><>< ====== ><>< # Example 1:
# print( wInput() )

wUpdate <- function(omegaInput){
    aaFull <- seq(1,20)
    dndsBank <- rep(NA, 20)
    editedQmat <- matrix(NA, 61, 61)
    aaZero <- aaFull[!(aaFull %in% omegaInput$aaPlus)]
    if(omegaInput$technique==1){
        aaCreate <- aaGauss
    }else{
        aaCreate <- aaGamma
    }
    
    for(a3 in seq(1,20)){
        synonymsID <- which(amino2codon == a3)
        if(a3 %in% aaZero){ s01x22 <- aaCreate(0, 0) }else{
            s01x22 <- aaCreate(omegaInput$vNvS, omegaInput$nsynVar) }
        csc01x61 <- codonCoeffs(s01x22)
        subMTX <- subsMatrix(csc01x61, omegaInput$psize)
        dndsBank[a3] <- dndsCalculator( codonFreq(csc01x61), subMTX)
        editedQmat[synonymsID,] <- subMTX[synonymsID,]
    }
    cqOutput <- list(qFrame=editedQmat, dnds=dndsBank)
    return(cqOutput)
}
## ><>< ======== ><>< ## Example
# omgUpdated <- wUpdate( wInput() )
# print(omgUpdated)

sitesim.omega <- function(parameterz, nodeLength,
                        popSize=NA, ntaxa, s01x22=NA){
    wOutput <- wUpdate( parameterz )
    if(parameterz$technique==1){
        aaCreate <- aaGauss
    }else{
        aaCreate <- aaGamma
    }
    s01x22 <- aaCreate(parameterz$vNvS, parameterz$nsynVar)
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
# ><>< ============ ><>< # Example
# wTest <- sitesim(wInput(), 0.10, 64)
# print( wTest )

alignsim.omega <- function(adaptIn, seqIn, modelIn=NULL, filename=NA){
    nsite <- seqIn$sites; ntaxons <- seqIn$taxa
    nodeLen <- seqIn$bl; treeData <- seqIn$t3
    alignment <- array(NA, c(ntaxons, nsite))
    dnds_matrix <- array(NA, c(20, nsite))
    
    for(b0 in seq(1,nsite)){
        simOut <- sitesim(adaptIn, nodeLen, NA, ntaxons)
        alignment[,b0] <- simOut$codons; dnds_matrix[,b0] <- simOut$dnds
    }
    aWord <- ifelse(adaptIn$technique==1, "Gauss", "Gamma")
    nzaa <- paste0("(", paste0(adaptIn$aaPlus,collapse=","), ")")
    notes <- paste0("Omega-Based Parameters: method=", aWord, " vNvS=",
        adaptIn$vNvS, " nsynVar=", adaptIn$nsynVar, " NonZero-AA=", nzaa)
    commentText <- paste0(seqIn$w, ";   ", notes)
    cdnSEQs <- seqDframe(alignment)
    if(is.null(filename)){
        coloredSQS <- seqColored(cdnSEQs)}else{coloredSQS <- NA}
    if(!is.na(filename)){
        empty <- seqWriter(alignment, treeData, commentText, filename)
    }
    return(list(seqs=alignment, DNDS=dnds_matrix,
                aInfo=commentText, cseq=cdnSEQs, seqCOL=coloredSQS))
}
# ><>< ======== ><>< # Example:
# omgEntry <- alignsim( wInput(), seqDetails(), NULL, "~/Desktop/seq.nex")
# print(omgEntry)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #