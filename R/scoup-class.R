# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # "scoup" Object Generics and Methods
setClass("scoup",
    representation(seqs="matrix", DNDS="matrix",
        aInfo="character", cseq="data.frame", seqCOL="ANY"))

setMethod("show", "scoup", function(object) {
    cat("\n", is(object)[[1]])
    cat("\nContains a codon sequence alignment with ", ncol(object@seqs),
        " sites and " , nrow(object@seqs), " extant taxa.\n\n")
})

setMethod("cseq", "scoup", function(x) x@cseq)
setMethod("dNdS", "scoup", function(x) x@DNDS)
setMethod("seqs", "scoup", function(x) x@seqs)
setMethod("aInfo", "scoup", function(x) x@aInfo)
setMethod("seqCOL", "scoup", function(x) x@seqCOL)

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #