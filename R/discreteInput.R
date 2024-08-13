# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><             Discretely Changing Landscape Simulation             ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 22 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

discreteInput <- function(defList=list()){
    aaParams <- defList$p02xnodes[c(1,2),]
    app <- defList$technique; if(is.null(app)) app <- 2
    nID <- defList$nodeIndex; if(is.null(nID)) nID <- 1
    xMD <- defList$leafModel; if(is.null(xMD)) xMD <- NA
    pnum <- defList$pSize; if(is.null(pnum)) pnum <- 1000
    if(is.null(aaParams)) aaParams <- rbind(rep(1,4),rep(1e-05,4))
    rownames(aaParams) <- c("vNvS","nsynVar")
    names(aaParams) <- NULL
    names(app) <- NULL;  names(nID) <- NULL
    names(pnum) <- NULL; names(xMD) <- NULL
    defEntry <- list(p02xnodes=aaParams, technique=app,
        nodeIndex=nID, psize=pnum, t3mdl=xMD)
    class(defEntry) <- "discrete"
    return(defEntry)
}
# ><>< ====== ><>< # Example 1:
# print( discreteInput() )

aaSCupdate.discrete <- function(parameters, oldSC=NA, bLength=NA){
    methd <- parameters$technique
    p02xnodes <- parameters$p02xnodes
    nodeIndex <- parameters$nodeIndex
    newParams <- p02xnodes[,nodeIndex]
    if(methd == 1){ scFunc <- aaGauss }else{ scFunc <- aaGamma }
    newSC <- scFunc(newParams["vNvS"], newParams["nsynVar"])
    return( abs(newSC) )
}
# ><>< ====== ><>< # Example 1:
# disparams <- rbind(vNvS=rep(0,4),nsynVar=rep(1e-04,4))
# defTest <- list(p02xnodes=disparams, nodeIndex=1, techniq=1)
# defInput <- discreteInput(defTest)
# outs <- aaSCupdate(defInput)
# print(outs)

sitesim.discrete <- function(parameterz, nodeLength,
                            popSize=NA, ntaxa=NA, s01x22=NA){
    popSize <- parameterz$psize
    s01x22 <- aaSCupdate(parameterz)
    parentCodon <- initSeq( s01x22 )
    n_nodes <- round(ncol(parameterz$p02xnodes), 0)
    seqVector <- c(parentCodon); dndsVec <- c()
    
    for(a4 in seq(1,n_nodes)){
        new_seq_vec <- rep(NA, 2^a4)
        parameterz[["nodeIndex"]] <- a4
        f01x22 <- aaSCupdate(parameterz)
        s01x22 <- f01x22; new_csc_vec <- codonCoeffs(s01x22)
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
# ><>< ============ ><>< # Example
# scp <- list(p02xnodes=rbind(rep(1e-05,8),rep(1e-02,8)))
# discreteparams <- discreteInput(scp); nodeL <- .30
# simsite1 <- sitesim(discreteparams, nodeL)
# print( simsite1); print( amino2codon[simsite1$codons])

alignsim.discrete <- function(adaptIn, seqIn, modelIn=NULL, filename=NA){
    theMtrx <- adaptIn$p02xnodes; nNodes <- round(ncol(theMtrx), 0)
    nsite <- seqIn$sites; ntaxons <- 2^nNodes; nodeLen <- seqIn$bl
    treeData <- biTree(ntaxons, nodeLen, adaptIn$t3mdl)
    dnds_matrix <- array(NA, c(nNodes, nsite))
    alignment <- array(NA, c(ntaxons, nsite))
    xMd <- theMtrx["vNvS",ncol(theMtrx)] > 0
    
    for(b0 in seq(1,nsite)){
        simOut <- sitesim(adaptIn, nodeLen)
        alignment[,b0] <- simOut$codons; dnds_matrix[,b0] <- simOut$dnds
    }
    aWord <- ifelse(adaptIn$technique==1, "Gauss", "Gamma")
    mString <- paste( apply(theMtrx, 2,
        function(a) paste(a, collapse=",")), collapse="|")
    notes <- paste0("Discrete Parameters: method=", aWord,
        " node-wise=(vNvS,nsynVar|", mString, ")")
    commentText <- paste0(seqIn$w, ";   ", notes)
    cdnSEQs <- seqDframe(alignment)
    if(!is.na(filename)){
        empty <- seqWriter(alignment, treeData, commentText, filename)
    }
    return(list(seqs=alignment, DNDS=dnds_matrix,
                aInfo=commentText, cseq=cdnSEQs))
}
# ><>< ======== ><>< # Example:
# nodeTheta <- matrix(rep(c(1e-04,1e-05),5), 2, byrow=FALSE)
# adptEntry <- discreteInput( list(p02xnodes=nodeTheta) )
# seqEntry <- seqDetails(); sFile <- "~/Desktop/seq.nex"
# outz <- alignsim(adptEntry, seqEntry, NULL, sFile)
# writeLines("\n"); print(summary(c(outz$D)))
# writeLines("\n"); print(outz$D[,1:6])

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #