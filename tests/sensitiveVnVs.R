# ><>< ================================================================= ><>< #
# ><><                    vN/vS Sensitivity Analyses.                    ><>< #
# ><>< ================================================================= ><>< #

# Make package accessible in R session
library(scoup)

# Number of extant taxa
xtant <- 64

# Number of data replications for each parameter combination
## Edited count was used for the results presented in article
simSize <- 1 # 50

# Variance of the non-synonymous selection coefficients
## Excluded values contributed to results presented in article
nsynVary <- c(0) # c(0, 0.001, 0.1)

# Ratio of the variance of the non-synonymous to synonymous coeff.
## Excluded values contributed to results presented in article
vNvSvec <- c(0) # c(0, 0.001,  1, 10)

# Sequence alignment size information
seqStat <- seqDetails(c(nsite=250, ntaxa=xtant))

# Iterate over all listed coefficient variance ratios
for(a in 1:length(vNvSvec)){

    # Iterate over all listed non-synonymous coefficients variance
    for(b in 1:length(nsynVary)){
      
        # Create appropriate simulation function ("omega") object
        adaptData <- wInput(list(vNvS=vNvSvec[a],nsynVar=nsynVary[b]))
        
        # Iterate over the specified number of replicates
        for(i in 1:simSize){
          
            # Execute simulation
            simulateSeq <- alignsim(adaptData, seqStat, NULL, NA)
        }
    }
}

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #