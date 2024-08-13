# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                  Ornstein-Uhlenbeck Simulation.                  ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 29 Apr, 2024                         ><>< #
# ><>< ================================================================ ><>< #

ouInput <- function(ouVector=0){
    pull <- ouVector["eMean"]; nPul <- ifelse(is.na(pull), 0, pull)
    vary <- ouVector["eVar"]; nVar <- ifelse(is.na(vary), 0.01, vary)
    oscil <- ouVector["Theta"]; nOsc <- ifelse(is.na(oscil), 0.01, oscil)
    scribe <- paste0("OU Parameters: Mu=",nPul,
                    "  Theta=",nOsc,"  Sigma.Sq=",nVar)
    names(nVar) <- NULL; names(nOsc) <- NULL; names(nPul) <- NULL
    ouEntry <- list(var=nVar, theta=nOsc, mu=nPul, words=scribe)
    class(ouEntry) <- "ou"
    return(ouEntry)
}
## ><>< ====== ><>< # Example 1:
# print( ouInput() )
## ><>< ====== ><>< # Example 2:
# ivec <- c(eVar=1e-02, Theta=10)
# ovec <- ouInput(ivec)
# print( ovec )

aaSCupdate.ou <- function(parameters, oldSC, bLength){
    newSC <- oldSC
    for(q0 in seq(1,20)){
        newSC[q0] <- ouEvolve(oldSC[q0], bLength,
            parameters$theta, parameters$var, parameters$mu) }
    return( abs(newSC) )
}
## ><>< ====== ><>< # Example 1:
# ouparams <- ouInput()
# veuxSC <- aaGauss(0, 1e-05)
# outs <- aaSCupdate(ouparams, veuxSC, 0.10)
# print(summary(veuxSC[1:20]-outs[1:20]))

sitesim.ou <- function(parameterz, nodeLength, popSize, ntaxa, s01x22){
    parentCodon <- initSeq( s01x22 )
    n_nodes <- ceiling( log(ntaxa, 2))
    seqVector <- c(parentCodon); dndsVec <- c()
    
    for(a4 in seq(1,n_nodes)){
        new_seq_vec <- rep(NA, 2^a4)
        f01x22 <- aaSCupdate(parameterz, s01x22, nodeLength)
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
## ><>< ============ ><>< # Example
# ouparams <- ouInput()
# ovec <- ouInput( c(eVar=1e-02, Theta=10))
# ntaxon <- 128; aasc <- aaGauss(1e-05,1e-05); nodeL <- .30
# simsite1 <- sitesim(ovec, nodeL, 1000, ntaxon, aasc)
# simsite0 <- sitesim(ouparams, nodeL, 1000, ntaxon, aasc)
# print( table(simsite1$c))#; print( table(simsite0$c))

alignsim.ou <- function(adaptIn, seqIn, modelIn, filename=NA){
    popsize <- modelIn$psize; scMeth <- modelIn$meth
    if(scMeth == 1){ scFunc <- aaGauss }else{ scFunc <- aaGamma }
    nNodes <- seqIn$nodes; nodeLen <- seqIn$bl; treeData <- seqIn$t3
    nsite <- seqIn$sites; ntaxons <- seqIn$taxa
    dnds_matrix <- array(NA, c(nNodes, nsite))
    alignment <- array(NA, c(ntaxons, nsite))
    
    for(b0 in seq(1,nsite)){
        aaSC <- scFunc(modelIn$vNvS, modelIn$nsynVar)
        simOut <- sitesim(adaptIn, nodeLen, popsize, ntaxons, aaSC)
        alignment[,b0] <- simOut$codons; dnds_matrix[,b0] <- simOut$dnds
    }; cdnSEQs <- seqDframe(alignment)
    if(is.null(filename)){
        coloredSQS <- seqColored(cdnSEQs)}else{coloredSQS <- NA}
    commentText <- paste0(adaptIn$w, ";   ", seqIn$w, ";   ", modelIn$w)
    if(is(filename, "character")){
        empty <- seqWriter(alignment, treeData, commentText, filename)
    }
    return(list(seqs=alignment, DNDS=dnds_matrix,
                aInfo=commentText, cseq=cdnSEQs, seqCOL=coloredSQS))
}
# ><>< ======== ><>< # Example:
# adptEntry <- ouInput()
# modEntry <- hbInput()
# seqEntry <- seqDetails()
# sFile <- "~/Desktop/seq.nex"
# outz <- alignsim(adptEntry, seqEntry, modEntry, sFile)
# writeLines("\n"); print(summary(c(outz$D)))
# writeLines("\n"); print(outz$D[,1:6])

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #