# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Simulate Balanced Bifurcating Tree
biTree <- function(ntaxa, bLength, terModel=NA){
    nNode <- ntaxa / 2
    tmNote <- ifelse(is.na(terModel), ":", paste0(terModel,":"))
    taxaTags <- paste0("S", sprintf("%03.0f",seq(1,ntaxa)), tmNote, bLength)
    while(nNode > 1){
        newTags <- c()
        for(i in seq(1,nNode)){ j <- i * 2
            newTags[i] <- paste0("(",taxaTags[j-1],",",
                                taxaTags[j], "):", bLength)
        }
        taxaTags <- newTags
        nNode <- nNode / 2
    }
    tree <- paste0("(", taxaTags[1], ",", taxaTags[2], "):", bLength, ";")
    return(tree)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #