# ><>< ========================================================================= ><>< #
# ><><   scoup: Simulate Codon Sequences with Darwinian Selection Incorporated   ><>< #
# ><><                     as an Ornstein-Uhlenbeck Process.                     ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                           Background Functions.                           ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                               23 Jan, 2024.                               ><>< #
# ><>< ========================================================================= ><>< #

# ><>< # Amino Acid and Codon Correspondence.
amino2codon <- c( 9, 12,  9, 12, 17, 17, 17, 17, 15, 16, 15, 16, 
   8,  8, 11,  8, 14,  7, 14,  7, 13, 13, 13, 13, 15, 15, 15, 15,
  10, 10, 10, 10,  4,  3,  4,  3,  1,  1,  1,  1,  6,  6,  6,  6,
  18, 18, 18, 18, 20, 20, 16, 16, 16, 16,  2, 19,  2, 10,  5, 10,  5)

# ><>< # Integer to codon transformer.
aminoacid <- c(AAA="K", AAC="N", AAG="K", AAT="N", ACA="T", ACC="T", ACG="T",
  ACT="T", AGA="R", AGC="S", AGG="R", AGT="S", ATA="I", ATC="I", ATG="M", ATT="I",
  CAA="Q", CAC="H", CAG="Q", CAT="H", CCA="P", CCC="P", CCG="P", CCT="P", CGA="R",
  CGC="R", CGG="R", CGT="R", CTA="L", CTC="L", CTG="L", CTT="L", GAA="E", GAC="D",
  GAG="E", GAT="D", GCA="A", GCC="A", GCG="A", GCT="A", GGA="G", GGC="G", GGG="G",
  GGT="G", GTA="V", GTC="V", GTG="V", GTT="V", TAC="Y", TAT="Y", TCA="S", TCC="S",
  TCG="S", TCT="S", TGC="C", TGG="W", TGT="C", TTA="L", TTC="F", TTG="L", TTT="F")

# ><>< # Codon Essentials
codonTriplets <- names(aminoacid)
nucleotides <- c("A", "C", "G", "T")
aacids <- LETTERS[c(1,seq(3,9),seq(11,14),seq(16,20),22,23,25)]

# ><>< # Substitution Matrix Indicator.
qmatID <- matrix(0, 61, 61); codonNuc <- matrix(NA, 61, 61)

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

# ><>< # Codon Mutation Matrix
mutationMatrix <- function(kappa=4, mu_ij=1){
  codon_mutate <- matrix(0, 61, 61)
  hky85 <- 0.25 * matrix(c(
    0, mu_ij, mu_ij*kappa, mu_ij, mu_ij, 0, mu_ij, mu_ij*kappa,
    mu_ij*kappa, mu_ij, 0, mu_ij, mu_ij, mu_ij*kappa, mu_ij, 0), 4)
  for(a3 in seq(1,61)){ for(a4 in seq(1,61)){ if(!is.na(codonNuc[a3,a4])){
    diffpost <- eval( parse(text=codonNuc[a3,a4]))
    codon_mutate[a3,a4] <- hky85[diffpost[1], diffpost[2]]
  } } }
  return(codon_mutate)
}
codon_m_matrix <- mutationMatrix()

# ><>< ==== ><>< # Non-Synonymous Codon Indicator Matrix
nonsynonymID <- synonymID <- matrix(0, 61, 61)
for(i in seq(1,61)){ for(j in seq(1,61)){
  synonymID[i, j] <- (amino2codon[i] == amino2codon[j]) * qmatID[i,j]
  nonsynonymID[i, j] <- (amino2codon[i] != amino2codon[j]) * qmatID[i,j]
} }

# ><>< ========================================================================= ><>< #
# ><><   scoup: Simulate Codon Sequences with Darwinian Selection Incorporated   ><>< #
# ><><                     as an Ornstein-Uhlenbeck Process.                     ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                             General Functions                             ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                               18 May, 2024.                               ><>< #
# ><>< ========================================================================= ><>< #

# ><>< # Define print() Method for AA Selection Coefficient
print.aminoSC <- function(aacoefs){
  cat("\nAmino acid selection coefficients:\n")
  cat(paste0(aacids[seq(1,4)],"=",round(aacoefs[seq(1,4)],3),", "), " ...\n")
  cat("\nSynonymous variance:\n",aacoefs[21],"\n")
  cat("\nNonsynonymous variance:\n",aacoefs[22],"\n\n")
}

# ><>< # Define print() Method for Codon Related Values
print.codonvalues <- function(codonSCF){
  cat("\n",paste0(codonTriplets[seq(1,6)],"=",
      round(codonSCF[seq(1,6)],3),", "), " ...\n\n")
}

# ><>< # Simulate Along Internal Node of a Tree
branchSimulate <- function(parentcodon, blength, q61x61e){
  init_time <- 0; nevolve <- 0
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
## ><>< ======== ><>< ## Example
# cvec <- codonCoeffs(aaGamma(1e-10,1e-01))
# smatx <- subsMatrix(cvec, 1000)
# newXter <- branchSimulate(15, 0.2, smatx)
# print(newXter)

# ><>< # Generate Codon Sequence at Root Node
initSeq <- function(s01x22){
  scoefvalues <- s01x22[seq(1,20)]
  coeffs_weights <- codonFreq(scoefvalues)
  aa_resid <- sample(seq(1,20), 1, F, coeffs_weights)
  codonResidue <- sample(which(amino2codon==aa_resid),1)
  return(codonResidue)
}
## ><>< ====== ><>< # Example:
# trial <- aaGauss(0.5, 1e-02)
# trialSeq <- initSeq(trial)
# print(aminoacid[trialSeq])

# ><>< # Merge Simulated Sequence
seqMerger <- function(alignmentMatrix, erase, t3model, filePrefix){
  nTaxa <- nrow(alignmentMatrix)
  nSite <- ncol(alignmentMatrix)
  codonspace <- names(aminoacid)
  t3File <- paste0(filePrefix, ".tre")
  seqFile <- paste0(filePrefix, ".txt")
  openText <- paste0("  ", nTaxa, "  ", nSite*3)
  write.table(openText, seqFile, !erase, F, row.names=F, col.names=F)
  taxaNames <- paste0(">S", sprintf("%03.0f",seq(1,nTaxa)))
  for(i in seq(1,nTaxa)){ newxters <- codonspace[ alignmentMatrix[i,] ]
    newSequence <- paste0(taxaNames[i],"   ",paste0(newxters, collapse=""))
    write.table(newSequence, seqFile, T, F, row.names=F, col.names=F)  }
  write.table("\n", seqFile, T, F, row.names=F, col.names=F)
  treeInfo <- biTree(nTaxa, 0.10, t3model[1])
  if(t3model[2]){ write.table(treeInfo, t3File, F, F, row.names=F, col.names=F)
     write.table("\n", t3File, T, F, row.names=F, col.names=F) }
}

# ><>< # Simulate Sequence Alignment
alignsim <- function(adaptIn, seqIn, modelIn, filename="seq.nex"){
  UseMethod("alignsim")
}

# ><>< # Update Selection Coefficient
aaSCupdate <- function(parameters, oldSC=NA, bLength=NA) UseMethod("aaSCupdate")

# ><>< # Simulate Genetic Sequence at a Site
sitesim <- function(parameterz, nodeLength, popSize=NA, ntaxa=NA, s01x22=NA){
  UseMethod("sitesim")
}

# ><>< ========================================================================= ><>< #
# ><><                              CODE ENDS HERE.                              ><>< #
# ><>< ========================================================================= ><>< #