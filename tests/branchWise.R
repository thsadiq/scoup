# ><>< ========================================================================= ><>< #
# ><><                           Branch-wise Analyses.                           ><>< #
# ><>< ========================================================================= ><>< #

# Make package accessible in R session
library(scoup)

# Number of internal nodes on the desired balanced tree
iNode <- 3

# Variance of non-synonymous selection coefficients
nsnV <- 0.01

# Number of data replications for each parameter combination
## Edited count was used for the results presented in article
nsim <- 1 # 50

# Ratio of the variance of the non-synonymous to synonymous coeff.
## Excluded values contributed to results presented in article
vNvSvec <- c(0) # c(0, 1e-06, 1e-03, 0.1, 1, 10, 100)

# Sequence alignment size information
seqsBwise <- seqDetails(c(nsite=1000, blength=0.10))

# Iterate over all listed coefficient variance ratios
for(h in 1:length(vNvSvec)){
  
  # Iterate over the specified number of replicates
  for(i in 1:nsim){
    
    # Create the parameter set applicable at each internal tree node
    scInput <- rbind(vNvS=c(rep(0,iNode-1),vNvSvec[h]),
                     nsynVar=rep(nsnV,iNode))
    
    # Create the applicable ("discrete") object for simulation function
    adaptBranch <- discreteInput(list(p02xnodes=scInput))
    
    # Execute simulation
    genSeq <- alignsim(adaptBranch, seqsBwise, NULL, NA)
  }
}

# ><>< ========================================================================= ><>< #
# ><><                              CODE ENDS HERE.                              ><>< #
# ><>< ========================================================================= ><>< #