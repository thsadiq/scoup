# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Simulate Balanced Bifurcating Tree
biTree <- function(ntaxa, bLength, terModel=NA){
    leafcheck <- ceiling( log(ntaxa, 2))
    newTaxa <- 2^leafcheck
    nNode <- newTaxa / 2
    if(newTaxa != ntaxa){
        statement <- paste0("Inappropriate extant taxa size specification!",
            " ntaxa = ", newTaxa, " was used instead to ensure that",
            " phylogeny is balanced.")
        warning(statement, call.=FALSE)
    }
    tmNote <- ifelse(is.na(terModel), ":", paste0(terModel,":"))
    taxaTags <- paste0("S", sprintf("%03.0f",seq(1,newTaxa)), tmNote, bLength)
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