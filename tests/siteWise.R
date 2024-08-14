# ><>< ================================================================= ><>< #
# ><><                        Site-wise Analyses.                        ><>< #
# ><>< ================================================================= ><>< #

# Make package accessible in R session
library(scoup)

# Number of codon sites
sitesize<- 100

# Variance of non-synonymous selection coefficients
nsynVary <- 0.01

# Number of extant taxa
taxasize <- 1024

# Sequence alignment size information
seqsEntry <- seqDetails(c(nsite=sitesize, ntaxa=taxasize))

# Create the applicable ("ou") object for simulation function
## eVar= OU asymptotic variance, Theta=OU reversion parameter
adaptEntry <- ouInput(c(eVar=0.1,Theta=1))

# Ratio of the variance of the non-synonymous to synonymous coeff.
## Excluded values contributed to results presented in article
sratio <- c(0) # c(0, 1e-06, 1e-03, 0.1, 1, 10, 1000)

# Iterate over all listed coefficient variance ratios
for(a0 in 1:length(sratio)){
  
    # OU landscape shift parameters
    mValues <- hbInput(c(vNvS=sratio[a0], nsynVar=nsynVary))
    
    # Execute simulation
    simSeq <- alignsim(adaptEntry, seqsEntry, mValues, NA)
}

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #