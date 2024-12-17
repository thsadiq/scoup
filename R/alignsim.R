# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                       Background Data Sets                       ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #

# ><>< # Amino Acid and Codon Correspondence.
amino2codon <- c( 9, 12,  9, 12, 17, 17, 17, 17, 15, 16, 15, 16, 
    8,  8, 11,  8, 14,  7, 14,  7, 13, 13, 13, 13, 15, 15, 15, 15,
    10, 10, 10, 10,  4,  3,  4,  3,  1,  1,  1,  1,  6,  6,  6,  6,
    18, 18, 18, 18, 20, 20, 16, 16, 16, 16,  2, 19,  2, 10,  5, 10,  5)

# ><>< # Integer to codon transformer.
aminoacid <- c(AAA="K", AAC="N", AAG="K", AAT="N", ACA="T", ACC="T",
    ACG="T", ACT="T", AGA="R", AGC="S", AGG="R", AGT="S", ATA="I", ATC="I",
    ATG="M", ATT="I", CAA="Q", CAC="H", CAG="Q", CAT="H", CCA="P", CCC="P",
    CCG="P", CCT="P", CGA="R", CGC="R", CGG="R", CGT="R", CTA="L", CTC="L",
    CTG="L", CTT="L", GAA="E", GAC="D", GAG="E", GAT="D", GCA="A", GCC="A",
    GCG="A", GCT="A", GGA="G", GGC="G", GGG="G", GGT="G", GTA="V", GTC="V",
    GTG="V", GTT="V", TAC="Y", TAT="Y", TCA="S", TCC="S", TCG="S", TCT="S",
    TGC="C", TGG="W", TGT="C", TTA="L", TTC="F", TTG="L", TTT="F")

# ><>< # Codon Essentials
codonTriplets <- names(aminoacid)
nucleotides <- c("A", "C", "G", "T")
aacids <- LETTERS[c(1,seq(3,9),seq(11,14),seq(16,20),22,23,25)]

# ><>< # Substitution Matrix Indicator.
qmatID <- matrix(0, 61, 61)
codonNuc <- matrix(NA, 61, 61)

for(i in seq(1,61)){
    for(j in seq(1,61)){
        if(i != j){
            wild <- unlist( strsplit(codonTriplets[i],""))
            target <- unlist( strsplit(codonTriplets[j],""))
            ndiffer <- wild == target
            if(sum(ndiffer) == 2){
                qmatID[i,j] <- 1
                nucpos <- which(!ndiffer)
                p1 <- which(nucleotides %in% wild[nucpos])
                p2 <- which(nucleotides %in% target[nucpos])
                codonNuc[i,j] <- paste0("c(", p1, ",", p2, ")")
} } } }

# ><>< ==== ><>< # Non-Synonymous Codon Indicator Matrix
nonsynonymID <- synonymID <- matrix(0, 61, 61)
for(i in seq(1,61)){ for(j in seq(1,61)){
    synonymID[i, j] <- (amino2codon[i] == amino2codon[j]) * qmatID[i,j]
    nonsynonymID[i, j] <- (amino2codon[i] != amino2codon[j]) * qmatID[i,j]
} }

# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                       Background Functions                       ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #

# ><>< # Simulate Along Internal Node of a Tree
branchSimulate <- function(parentcodon, blength, q61x61e){
    init_time <- 0
    nevolve <- 0
    while(init_time < blength){
        init_codon <- parentcodon
        weights <- q61x61e[init_codon,]
        wait_time <- rexp(1, sum(weights))
        new_time <- init_time + wait_time
        if(new_time <= blength){ nevolve <- nevolve + 1
            parentcodon <- sample(seq(1,61), 1, prob=weights) }
        init_time <- new_time
    }
    return(c(codon=parentcodon, evolved=nevolve))
}

# ><>< # Generate Codon Sequence at Root Node
initSeq <- function(s01x22){
    scoefvalues <- codonCoeffs(s01x22)
    coeffs_weights <- coeffs( codonFreq(scoefvalues))
    codonResidue <- sample(seq(1,61), 1, FALSE, coeffs_weights)
    return(codonResidue)
}

# ><>< # Merge Simulated Sequence
seqMerger <- function(alignmentMatrix, erase, t3model, filePrefix){
    nTaxa <- nrow(alignmentMatrix)
    nSite <- ncol(alignmentMatrix)
    codonspace <- names(aminoacid)
    t3File <- paste0(filePrefix, ".tre")
    seqFile <- paste0(filePrefix, ".txt")
    openText <- paste0("  ", nTaxa, "  ", nSite*3)
    write.table(openText, seqFile, !erase, FALSE, 
                row.names=FALSE, col.names=FALSE)
    taxaNames <- paste0(">S", sprintf("%03.0f",seq(1,nTaxa)))
    for(i in seq(1,nTaxa)){ newxters <- codonspace[ alignmentMatrix[i,] ]
        newSequence <- paste0(taxaNames[i],"   ",paste0(newxters,collapse=""))
        write.table(newSequence, seqFile, TRUE, FALSE,
                    row.names=FALSE, col.names=FALSE)
    }
    write.table("\n", seqFile, TRUE, FALSE, row.names=FALSE, col.names=FALSE)
    treeInfo <- biTree(nTaxa, 0.10, t3model[1])
    if(t3model[2]){ 
        write.table(treeInfo, t3File, FALSE, FALSE,
                    row.names=FALSE, col.names=FALSE)
        write.table("\n", t3File, TRUE, FALSE,
                    row.names=FALSE, col.names=FALSE)
    }
}

# ><>< # Erase redundant names
nameDel <- function(y){
    names(y) <- NULL
    return(y)
}

# ><>< # Warn of Redundant Input Labels
warner <- function(pryIn, fullNames){
    nameless <- any(names(pryIn) == "")
    noName <- length(pryIn) - length(names(pryIn)) > 0
    rule0 <- (length(pryIn) !=1) | all(unlist(pryIn) != 0)
    regulate <- (nameless | noName) & rule0
    xtraTag <- NULL
    if(regulate)  xtraTag <- "*unnamed*"
    inNames <- c(xtraTag, names(pryIn))
    chck1 <- inNames %in% fullNames
    if(!all(chck1)){
        faux <- inNames[which(!chck1)]
        communicate <- paste0("The inputs listed below were ignored ",
            "because they are invalid entries of the discrete ",
            "model.\n\t", paste(faux, collapse="\t"))
        warning(communicate, call. = FALSE)
    }
}

# ><>< # Identify Input Distribution for Initial Coefficient Sampling
detectApp <- function(x){ x1 <- x@sampler
    tch <- ifelse(x1==1, "Normal distribution\n",  "Gamma distribution\n")
    message(tch)
}

# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                        Generic Functions.                        ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #

# ><>< # discrete, ou
setGeneric("aaSCupdate", 
    function(parameters, ...) standardGeneric("aaSCupdate") )
# ><>< # discrete, omega, ou
setGeneric("alignsim", 
    function(adaptIn, seqIn, ...) standardGeneric("alignsim"))
setGeneric("sitesim", 
    function(parameters, nodeLength, ...) standardGeneric("sitesim") )
# ><>< # scoup
setGeneric("cseq", function(x) standardGeneric("cseq"))
setGeneric("dNdS", function(x) standardGeneric("dNdS"))
setGeneric("seqs", function(x) standardGeneric("seqs"))
setGeneric("aInfo", function(x) standardGeneric("aInfo"))
setGeneric("seqCOL", function(x) standardGeneric("seqCOL"))
# ><>< # hbParameters, omega
setGeneric("vNvS", function(x) standardGeneric("vNvS"))
setGeneric("nsynVar", function(x) standardGeneric("nsynVar"))
# ><>< # aminoSC, codonvalues
setGeneric("freqs", function(x) standardGeneric("freqs"))
setGeneric("coeffs", function(x) standardGeneric("coeffs"))
# ><>< # discrete, hbParameters, omega
setGeneric("kappa", function(x) standardGeneric("kappa"))
setGeneric("effpop", function(x) standardGeneric("effpop"))
setGeneric("hky85mu", function(x) standardGeneric("hky85mu"))
setGeneric("sampler", function(x) standardGeneric("sampler"))
# ><>< # discrete, omega
setGeneric("lscape", function(x) standardGeneric("lscape"))
# ><>< # aminoSC
setGeneric("synVar", function(x) standardGeneric("synVar"))
# ><>< # ou
setGeneric("asymVar", function(x) standardGeneric("asymVar"))
setGeneric("asymMean", function(x) standardGeneric("asymMean"))
setGeneric("reversion", function(x) standardGeneric("reversion"))
# ><>< # seqParameters
setGeneric("taxa", function(xo) standardGeneric("taxa"))
setGeneric("sites", function(xo) standardGeneric("sites"))
setGeneric("nodes", function(xo) standardGeneric("nodes"))
setGeneric("branchL", function(xo) standardGeneric("branchL"))
setGeneric("details", function(xo) standardGeneric("details"))
setGeneric("phylogeny", function(xo) standardGeneric("phylogeny"))

# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                        Class Definitions.                        ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
setClass("seqParameters",
    representation(sites="numeric", taxa="numeric", nodes="numeric",
        branchL="numeric", phylogeny="character", details="character"))

setClass("discrete",
    representation(lscape="matrix", sampler="numeric",
        nodeIndex="numeric", psize="numeric",
        t3mdl="character", kappa="numeric", mrate="numeric"))

setClass("omega", representation(nsynVar="numeric", psize="numeric",
    sampler="numeric", aaPlus="numeric", vNvS="numeric",
    kappa="numeric", mrate="numeric"))

setClass("ou", representation(var="numeric",
    theta="numeric", mu="numeric", words="character"))

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #