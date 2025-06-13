# ><>< ================================================================= ><>< #
# ><><     scoup:- Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process.          ><>< #
# ><>< ================================================================= ><>< #

# ><>< # Clear memory
# rm(list=ls())

# ><>< # Install package from the Bioconductor platform
    if (!require("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }
    # The following initialises usage of Bioc devel with the most recent
    # updates. It requires R version 4.5 or later. With R version 4.4,
    # use ~ version='devel' ~ instead
    BiocManager::install(version='devel')
    BiocManager::install("scoup")

# ><>< # Make package accessible in R session
library(scoup)

# ><>< # ======
# Combine sequences from >1 simulation executions to facilitate PAML analyses.
# - scoupObj: output of 'alignsim' from each simulation execution.
# - simID: index of the current execution i.e. order of alignment.
# - seqLoc: path to file where the merged data shoud be stored.
# ><>< # ======
mergeseq <- function(scoupObj, simID, seqLoc){
    newData <- cseq(scoupObj)
    nleaf <- nrow(newData)
    ndna <- ncol( seqs(scoupObj)) * 3
    apend <- ifelse(simID == 1, FALSE, TRUE)
    preSEQ <- paste0("  ", nleaf, "  ", ndna)
    write.table(preSEQ, seqLoc, apend, FALSE, row.names=FALSE, col.names=FALSE)
    for(a1 in 1:nleaf){
        sq0 <- paste0(">S", sprintf("%03.0f",a1), "   ", newData[a1,])
        write.table(sq0, seqLoc, TRUE, FALSE, row.names=FALSE, col.names=FALSE)
    }; write.table("\n", seqLoc, TRUE, FALSE, row.names=FALSE, col.names=FALSE)
    dndsVALUE <- mean( dNdS(scoupObj) )
    return(dndsVALUE)
}
# ><>< # ======

# ><>< # Specify data storage paths
# pryPath <- path.expand(file.path("~","suppl"))
# infoPath <- file.path(pryPath, "dataInfo.txt")
# dndsPath <- file.path(pryPath, "dnds.csv")

# ><>< # Number of extant taxa
xtant <- 8

# ><>< # Number of data replications for each parameter combination
simSize <- 20

# ><>< # Number of codon sites
seqLength <- 1000

# ><>< # Variance of the non-synonymous selection coefficients
nsynVary <- c(0, 0.10)

# ><>< # Ratio of the variance of the non-synonymous to synonymous coeff.
vNvSvec <- c(0, 5.00)

# ><>< # OU asymptotic variance (Sigma) value
eVary <- c(0.00, 0.34, 0.67, 1.00)

# ><>< # OU reversion parameter (Theta) value
eThta <- c(0.01, 0.34, 0.67, 1.00)

# ><>< # Sequence alignment size information
seqStat <- seqDetails(c(nsite=seqLength, ntaxa=xtant))

# ><>< # Save the evolutionary tree used for the simulation
t3newick <- biTree(xtant, branchL(seqStat))
# t3Path <- file.path(pryPath, "evolTree.txt")
# write.table(t3newick, t3Path, FALSE, FALSE, row.names=FALSE, col.names=FALSE)

# ><>< # Iterate over all listed coefficient variance ratios
for(a in seq(1,length(vNvSvec))){
    
    # ><>< # Iterate over non-synonymous variance values
    for(e in seq(1,length(nsynVary))){
        
        # ><>< # Iterate over all listed OU variance values
        for(b in seq(1,length(eVary))){
            
            # ><>< # Iterate over all listed OU reversion parameter values
            for(d in seq(1,length(eThta))){
                
                # ><>< # Create secondary path name
                secLoc <- paste0("vN", e, "vS", a, "ouS", b, "ouT", d)
                # seqPath <- file.path(pryPath, paste0(secLoc,".txt"))
                
                # ><>< # Create necessary ("omega") object
                hbrunoStat <- hbInput(c(vNvS=vNvSvec[a], nsynVar=nsynVary[e]))
                
                # ><>< # Create appropriate simulation function ("ou") object
                adaptStat <- ouInput(c(eVar=eVary[b],Theta=eThta[d]))
                
                # ><>< # Initialise dN/dS storage vector
                dndsVec <- numeric(simSize)
                
                # ><>< # Iterate over the specified number of replicates
                for(f in seq(1,simSize)){
                    
                    # ><>< # Execute simulation and save outputs
                    simulateSeq <- alignsim(adaptStat, seqStat, hbrunoStat, NA)
                    # saveseq <- mergeseq(simulateSeq, f, seqPath)
                    # dndsVec[f] <- saveseq
                }
                # ><>< # Print update message
                if(nsynVary[e]==0){
                    svar <- ifelse(vNvSvec[a]==0, 0, vNvSvec[a])
                }else{
                    svar <- ifelse(vNvSvec[a]==0, 0,
                        round(nsynVary[e]/vNvSvec[a],2))
                }
                message( paste0("Done: vS = ", svar, " | vN = ", nsynVary[e],
                    " | OU.Sigma = ", eVary[b], " | OU.Theta = ", eThta[d]) )
                
                # ><>< # Save details and outputs of simulation set
                # write.table( paste0(secLoc, "  ||  ", aInfo(simulateSeq)),
                #     infoPath,a+b+d+e!=4,FALSE,row.names=FALSE,col.names=FALSE)
                
                # dndsXter <- paste0(secLoc,",",paste(dndsVec,collapse=",  "))
                # write.table(dndsXter, dndsPath, a+b+d+e!=4,
                #             FALSE, row.names=FALSE, col.names=FALSE)
            }
        }
    }
}

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #