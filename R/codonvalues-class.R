# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # "codonvalues" Object Generics and Methods
setClass("codonvalues", representation(cdnums="numeric"))

setMethod("freqs", "codonvalues", function(x) nameDel(x@cdnums))
setMethod("coeffs", "codonvalues", function(x) nameDel(x@cdnums))
setMethod("show", "codonvalues", function(object) {
    cat("\n", is(object)[[1]])
    cat("\n\n",paste0(codonTriplets[seq(1,6)],"=",
        round(object@cdnums[seq(1,6)],3),", "), " ...\n\n")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #