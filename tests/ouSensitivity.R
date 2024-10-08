# ><>< ================================================================= ><>< #
# ><><                      OU Sensitivity Analyses                      ><>< #
# ><>< ================================================================= ><>< #

# Make package accessible in R session
library(scoup)

# Number of extant taxa
leaves <- 64

# Number of data replications for each parameter combination
## Edited count was used for the results presented in article
sims <- 1 # 50

# OU reversion parameter (Theta) value
## Excluded values contributed to results in published article
eThta <- c(0.01) # c(0.01, 0.1, 1)

# OU asymptotic variance value
## Excluded values contributed to results in published article 
eVary <- c(0.0001) # c(0.0001, 0.01, 1)

# OU landscape shift parameters
hbrunoStat <- hbInput(c(vNvS=1, nsynVar=0.01))

# Sequence alignment size information
seqStat <- seqDetails(c(nsite=250, ntaxa=leaves))

# Iterate over all listed OU variance values
for(g in 1:length(eVary)){
  
    # Iterate over all listed OU reversion parameter values
    for(h in 1:length(eThta)){

        # Create appropriate simulation function ("ou") object
        adaptStat <- ouInput(c(eVar=eVary[g],Theta=eThta[h]))

        # Iterate over the specified number of replicates
        for(i in 1:sims){

            # Execute simulation
            simData <- alignsim(adaptStat, seqStat, hbrunoStat, NA)
        }
    }
}

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #