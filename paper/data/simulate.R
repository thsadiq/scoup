# ><>< ================================================================= ><>< #
# ><><                    Protein Sequence Simulation                    ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><><  R code used to simulate the protein sequences discussed in the   ><>< #
# ><><                          JOSS manuscript                          ><>< #
# ><>< ================================================================= ><>< #

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
    write.table(preSEQ, seqLoc, apend, FALSE,
            row.names=FALSE, col.names=FALSE)
    wrt <- vapply(seq(1,nleaf), function(w){
        sq0 <- paste0(">S", sprintf("%03.0f",w), "   ", newData[w,])
        write.table(sq0,seqLoc,TRUE,FALSE,row.names=FALSE,col.names=FALSE)
        return(0)
        }, FUN.VALUE=0)
    write.table("\n", seqLoc, TRUE, FALSE, row.names=FALSE, col.names=FALSE)
    dndsVALUE <- mean( dNdS(scoupObj) )
    return(dndsVALUE)
}
# ><>< # ======

# Make simulation package accessible in R session
library(scoup)

# Number of tree leaves
leaves <- 8

# Number of codon sites
sitesize <- 1000

# Number of data replications for each tree length
sims <- 10

# OU reversion parameter (Theta) value
eThta <- 0.1

# OU asymptotic variance value
eVary <- 0.1

# Sequence alignment size information
seqStat <- seqDetails(c(nsite=sitesize, ntaxa=leaves))

# Specify data storage paths
pryPath <- tempdir()

# Save the evolutionary tree used for the simulation
t3newick <- biTree(leaves, 0.10)
t3Path <- file.path(pryPath, "simtree.txt")
write.table(t3newick, t3Path, FALSE, FALSE, row.names=FALSE, col.names=FALSE)

# Synonymous variance values
synvary <- c(0.00, 0.01)
stags <- paste0("s",sprintf("%03.0f", synvary*1000))

# Non-synonymous variance value
nsynvary <- c(c(1,5,9)/1000, c(3,7)/100, c(1,5,9)/10)
nbind <- paste0("n",sprintf("%03.0f", nsynvary*1000))

# Iterate over all listed synonymous variance values
for(h in seq(1,length(synvary))){

    # Create storage for dN/dS estimates
    dndsPath <- file.path(pryPath, paste0(stags[h],"dnds.csv"))
    dndsarray <- matrix(NA, ncol=length(nsynvary), nrow=sims)
    colnames(dndsarray) <- nbind

    # Create appropriate simulation function ("ou") object
    adaptStat <- ouInput(c(eVar=eVary,Theta=eThta))

     # Iterate over the different branch lengths
    for(k in seq(1,length(nsynvary))){

        # OU landscape shift parameters
        newvnvs <- ifelse(synvary[h] == 0, 0, nsynvary[k]/synvary[h])
        hbrunoStat <- hbInput(c(vNvS=newvnvs, nsynVar=nsynvary[k]))

        # Create sequence file
        seqPath <- file.path(pryPath, paste0(stags[h],nbind[k],".txt"))

        # Iterate over the specified number of replicates
        for(i in seq(1,sims)){ 

            # Execute simulation
            simData <- alignsim(adaptStat, seqStat, hbrunoStat, NA)
            dndsavg <- mergeseq(simData, i, seqPath)
            dndsarray[i,k] <- dndsavg
            message( paste0("Generated sequences for: ",
                stags[h], " x ", nbind[k], " x sims", i))
        }
    }
    write.table(dndsarray, dndsPath, FALSE, FALSE, ";", row.names=FALSE)
}
message( paste0("\nscoup simulation completed and outputs",
    " are saved in:\n\t", pryPath, "\n"))

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #
