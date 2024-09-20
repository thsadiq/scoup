# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # "seqParameters" Object Methods
setMethod("taxa", "seqParameters", function(xo) xo@taxa)
setMethod("sites", "seqParameters", function(xo) xo@sites)
setMethod("nodes", "seqParameters", function(xo) xo@nodes)
setMethod("branchL", "seqParameters", function(xo) xo@branchL)
setMethod("details", "seqParameters", function(xo) xo@details)
setMethod("phylogeny", "seqParameters", function(xo) xo@phylogeny)
setMethod("show", "seqParameters", function(object) {
    cat("\n", is(object)[[1]])
    cat("\nContains details of the structure on the simulated alignment.\n\n")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #