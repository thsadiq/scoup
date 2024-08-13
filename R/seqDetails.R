# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                    <- seqDetails -> Function.                    ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 18 May, 2024                         ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Create Sequence Details Input
seqDetails <- function(seqVector=0){
    if(is(seqVector, "character")){
        idModel <- which(names(seqVector)=="terModel")
        newVector <- as.numeric(seqVector[-idModel])
    }else{ newVector <-  seqVector}
    xModel <- seqVector["terModel"]
    ntaxons <- newVector["ntaxa"]; names(ntaxons) <- NULL
    if(is.na(ntaxons)){ taxa <- 64; nNodes <- 6 }else{
        nNodes <- ceiling( log(ntaxons, 2)); taxa <- 2^nNodes }
    nsite <- newVector["nsite"]; sites <- ifelse(is.na(nsite), 250, nsite)
    ndLength <- newVector["blength"]
    BL <- ifelse(is.na(ndLength), 0.10, ndLength)
    t3d <- biTree(taxa, BL, xModel)
    ttext <- paste0("Tree Info: branch.length=", BL)
    names(sites) <- NULL; names(BL) <- NULL
    seqEntry <- list(sites=sites, taxa=taxa, nodes=nNodes,
                    bl=BL, t3=t3d, words=ttext)
    class(seqEntry) <- "seqParameters"
    return(seqEntry)
}
## ><>< ====== ><>< # Example 1:
# print( seqDetails() )
## ><>< ====== ><>< # Example 2:
# putvec <- c(ntaxa=16, nsite=10, blength=0.20)
# extvec <- seqDetails(putvec)
# print( extvec )

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #