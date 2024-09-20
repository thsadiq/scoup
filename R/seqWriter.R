# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Write Simulated Sequence to File
seqWriter <- function(alignmentMatrix, treeInfo=NA, addText="", fileTag=NULL){
    nTaxa <- nrow(alignmentMatrix)
    nSite <- ncol(alignmentMatrix)
    codonspace <- names(aminoacid)
    if(is.na(treeInfo)) treeInfo <- biTree(nTaxa, 0.10)
    treetext <- paste0("\n\nbegin trees;\n  tree TREE1 = ", treeInfo, "\nend;")
    openText <- paste0("#NEXUS\n[",addText, "]\nbegin data;\ndimensions ntax=",
        nTaxa, " nchar=", nSite*3,
        ";\nformat datatype=dna missing=? gap=-;\nmatrix\n")
    closeText <- paste0("\n;\nend;")
    
    if(is.null(fileTag)){
        seqFile <- file.path(tempdir(),"cranrSeqs.nex")}else{seqFile<- fileTag}
    write.table(openText, seqFile, FALSE, 
        FALSE, row.names=FALSE, col.names=FALSE)
    taxaNames <- paste0(" S", sprintf("%03.0f",seq(1,nTaxa)))
    
    for(i in seq(1,nTaxa)){
        newxters <- codonspace[ alignmentMatrix[i,] ]
        newSequence <- paste0(taxaNames[i],"   ",paste0(newxters, collapse=""))
        write.table(newSequence, seqFile, TRUE, FALSE,
                    row.names=FALSE, col.names=FALSE)
    }
    write.table(closeText, seqFile, TRUE, FALSE,
                row.names=FALSE, col.names=FALSE)
    write.table(treetext, seqFile, TRUE, FALSE,
                row.names=FALSE, col.names=FALSE)
    statement <- paste0("Simulated sequence was written to file: ",seqFile)
    message(statement)
}

# ><>< # Convert Sequence to Data Frame
seqDframe <- function(alignmentMatrix){
    nTaxa <- nrow(alignmentMatrix)
    codonspace <- names(aminoacid)
    seqBank <- matrix(NA, nTaxa, 1)
    colnames(seqBank) <- c("Sequence")
    rownames(seqBank) <- paste0(" S", sprintf("%03.0f",seq(1,nTaxa)))
    
    for(i in seq(1,nTaxa)){
        newxters <- codonspace[ alignmentMatrix[i,] ]
        seqBank[i,] <- paste0(newxters, collapse="")
    }
    return( as.data.frame(seqBank) )
}

# ><>< # Use Biostrings to Generate Fancy Sequences
seqColored <- function(framedSEQ){
    vectorSEQ <- c(framedSEQ)[[1]]
    sqNames <- paste0("S", sprintf("%03.0f",seq(1,nrow(framedSEQ))))
    coloredOUT <- DNAStringSet(vectorSEQ)
    names(coloredOUT) <- sqNames
    return(coloredOUT)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #