# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Create Sequence Details Input
seqDetails <- function(seqVector=0){
    warner(seqVector,c("ntaxa","nsite","blength","terModel"))
    criterion <- is(seqVector, "character")
    if(criterion){
        idModel <- which(names(seqVector)=="terModel")
        newVector <- as.numeric(seqVector[-idModel])
        names(newVector) <- names(seqVector[-idModel])
    }else{ newVector <-  seqVector}
    xModel <- seqVector["terModel"]
    ntaxons <- newVector["ntaxa"]
    names(ntaxons) <- NULL
    if(is.na(ntaxons)){
        taxa <- 64
        nNodes <- 6
    }else{
        nNodes <- ceiling( log(ntaxons, 2))
        taxa <- 2^nNodes
    }
    nsite <- newVector["nsite"]
    sites <- ifelse(is.na(nsite), 250, nsite)
    ndLength <- newVector["blength"]
    BL <- ifelse(is.na(ndLength), 0.10, ndLength)
    t3d <- biTree(taxa, BL, xModel)
    ttext <- paste0("Tree Info: branch.length=", BL)
    names(sites) <- NULL
    names(BL) <- NULL
    seqEntry <- new("seqParameters", sites=sites, taxa=taxa,
        nodes=nNodes, branchL=BL, phylogeny=t3d, details=ttext)
    return(seqEntry)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #