# ><>< ================================================================= ><>< #
# ><><                       Branch-wise Analyses.                       ><>< #
# ><><                             ~~~~~~~~~                             ><>< #
# ><>< NOTE:- The edited (commented out) values were deemed necessary to ><>< #
# ><><        avoid unnecessary delays during package installation.      ><>< #
# ><>< ================================================================= ><>< #

# Make package accessible in R session
library(scoup)

# Number of internal nodes on the desired balanced tree
iNode <- 3

# Variance of non-synonymous selection coefficients
nsnV <- 0.01

# Number of data replications for each parameter combination
nsim <- 1 # 50

# Ratio of the variance of the non-synonymous to synonymous coeff.
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
        genSeq <- alignsim(adaptBranch, seqsBwise, NA)
    }
}

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #