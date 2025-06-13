# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# Clear memory
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

# Make package accessible in R session
library(scoup)

# Create toy function for accurate timing
testFunc <- function(){
    adaptEnrty <- ouInput()
    modelEntry <- hbInput()
    sqEntry <- seqDetails()
    seqData <- alignsim(adaptEnrty, sqEntry, modelEntry)
    return(NULL)
}

# Specify data storage paths
pryPath <- path.expand(file.path("~","suppl"))
timesPath <- file.path(pryPath, "times.csv")

# Implement execution time simultions 
simsize <- 50
tvec <- numeric(simsize)
for(a in seq(1,simsize)) tvec[a] <- system.time(testFunc())[3]
write.table(tvec, timesPath, FALSE, FALSE, row.names=FALSE, col.names=FALSE)

# ><>< ================================================================= ><>< #
# ><><                          CODE ENDS HERE.                          ><>< #
# ><>< ================================================================= ><>< #